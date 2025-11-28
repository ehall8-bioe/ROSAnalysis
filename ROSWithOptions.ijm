initializeSciJavaParameters();

#@ String (label="Cell Membrane Channel", description="Please use whatever information is used to identify the images",persist=false) membrane
#@ String (label="Channel for Quantificantion",persist=false) channel
#@ String (label = "What part of the cell should be quantified?",choices={"Whole Cell","Nucleus and Cytosol Seperately"}, style="radioButtonHorizontal",persist=false) area
#@ String (label = "Threshold for Nucleus",choices = {"None","Default","Huang","Intermodes","IsoData","IJ_IsoData","Li","MaxEntropy","Mean", "MinError","Minimum","Moments","Otsu","Percentile","RenyiEntropy","Shanbhag","Triangle","Yen"},style="listBox",persist=false) nuc_thr
#@ String (label = "Threshold for Cell Body",choices = {"Default","Huang","Intermodes","IsoData","IJ_IsoData","Li","MaxEntropy","Mean", "MinError","Minimum","Moments","Otsu","Percentile","RenyiEntropy","Shanbhag","Triangle","Yen"},style="listBox",persist=false) cyto_thr
#@ String (label = "Enhance Contrast",choices={"Nucleus","Cell Body","Both","Neither"}, style="radioButtonHorizontal",persist=false) ECloc
#@ Double (label="Contrast Enhance", style="slider,format:0.000",persist = false, value = 0.35, min=0, max=1, stepSize=0.05) EC
#@ String (label = "Subract Background",choices={"Nucleus","Cell Body","Both","Neither"}, style="radioButtonHorizontal",persist=false) SBloc
#@ Double (label = "Rolling ball radius for Background Subtraction", value = 50, persist=false) SBnum
#@ Double (label = "Minimum Cell Size (pix^2)",value=30) cellsize
#@ File (label ="Input Directory",style="directory") in
#@ File (label ="Output Directory",style="directory") out
close("*");
if (area=="Whole Cell") {
	list = getFileList(in);
	Table.create("Cell Intensity");
	for(i=0; i<list.length; i++){ //Iterate through file list to analyze all images
		name = list[i];
		lcname = toLowerCase(name);
		inputDir = in+"/";
		outputDir = out+"/";
		if(indexOf(lcname, toLowerCase(membrane))>=0 && endsWith(name, ".tiff")){
			run("Set Measurements...", "mean display redirect=None decimal=3"); //makes sure program measures fluorescence intensity
			fileLoc = inputDir+name;
			open(fileLoc);
			selectImage(name);
			name1 = substring(name,0,indexOf(name, membrane)-1);
			coi = name1+"_"+channel+".tiff";
			if(ECloc =="Cell Body" || ECloc == "Both"){
				run("Enhance Contrast...", "saturated="+EC+" normalize");
			}
			if(SBloc =="Cell Body" || SBloc =="Both"){
				run("Subtract Background...", "rolling="+SBnum);
			}
			setAutoThreshold(cyto_thr+" dark");
			run("Analyze Particles...", "size="+cellsize+"-Infinity clear add");
			selectImage(name);
			close();
			open(inputDir+coi);
			roiManager("Measure");
			data = Table.getColumn("Mean","Results");
			Table.setColumn(name1,data,"Cell Intensity");
			selectWindow("Results");
			run("Close");
			close("*");
		}
	}
	selectWindow("Cell Intensity");
	Table.save(outputDir+"Cell Intensity.csv");
	run("Close");
	
}
if (area=="Nucleus and Cytosol Seperately"){
	list = getFileList(in);
	Table.create("Nuclear Intensity");
	Table.create("Cytosolic Intensity");
	for(i=0; i<list.length; i++){ //Iterate through file list to analyze all images
		name = list[i];
		lcname = toLowerCase(name);
		inputDir = in+"/";
		outputDir = out+"/";
		if(indexOf(lcname, "dapi")>=0 && endsWith(name, ".tiff")){
			run("Set Measurements...", "mean display redirect=None decimal=3"); //makes sure program measures fluorescence intensity
			fileLoc = inputDir+name;
			open(fileLoc);
			selectImage(name);
			name1 = substring(name,0,indexOf(name, "dapi")-1);
			cell = name1+"_"+membrane+".tiff";
			coi = name1+"_"+channel+".tiff";
			if(ECloc =="Nucleus" || ECloc == "Both"){
				run("Enhance Contrast...", "saturated="+EC+" normalize");
			}
			if(SBloc =="Nucleus" || SBloc =="Both"){
				run("Subtract Background...", "rolling="+SBnum);
			}
			setAutoThreshold(nuc_thr+" dark"); // Threshold type may need to be adjusted depending on the imaging settings
			run("Analyze Particles...", "size="+cellsize+"-Infinity clear add");
			count = roiManager("count");
			selectImage(name);
			close();
			open(inputDir+coi);
			roiManager("Measure"); 
			nucval = Table.getColumn("Mean","Results");
			Table.setColumn(name1,nucval,"Nuclear Intensity");
			selectWindow("Results");
			run("Close");
			open(inputDir+cell);
			if(ECloc =="Cell Body" || ECloc == "Both"){
				run("Enhance Contrast...", "saturated="+EC+" normalize");
			}
			if(SBloc =="Cell Body" || SBloc =="Both"){
				run("Subtract Background...", "rolling="+SBnum);
			}
			setAutoThreshold(cyto_thr+" dark");
			run("Analyze Particles...", "size="+cellsize+"-Infinity add");
			selectImage(cell);
			close();
			count2 = roiManager("count");
			A = newArray(1);
			C = newArray(1);
			for(k=count; k<count2; k++){
				for(m=0; m<count; m++){
					roiManager("select", newArray(k,m));
					roiManager("AND");
					if(selectionType()>-1){
						A[A.length] = k;
						C[C.length] = m;
					}
				}
			}
			A = Array.deleteIndex(A, 0); //Remove values used to start the new arrays
			C = Array.deleteIndex(C, 0);
			for(n=0; n< A.length; n++){ 
				forxor = newArray(1);
				E = newArray(1);
				y = A[n];
				forxor[0] = A[n];
				forxor[1] = C[n];
				for(p=n+1; p<A.length; p++){ // Group indexes for XOR compensates for multiple nuclei. 
					if(A[p] == y){ // If cell boundry overlaps with multiple nuclei, add both idexes to the function
						E[E.length] = p;
						forxor[forxor.length] = C[p]; 
					}
				}
				E = Array.deleteIndex(E,0);
				if(E.length>0){ // Remove repeats
					for(l=0; l<E.length; l++){
						x = E[l];
						C[x] = -1;
						A[x] = -1;
					}
				}
				C = Array.deleteValue(C, -1);
				A = Array.deleteValue(A, -1);
				if(forxor.length >1){ // Combine ROIs using xor
					roiManager("select", forxor);
					roiManager("XOR");
					roiManager("add");
				}
			}
			count3 = roiManager("count");
			roiManager("Select",Array.getSequence(count2)); // Delete the original Nuclear and Cell Body ROIs
			roiManager("delete");
			roiManager("Measure"); // Measure all the XOR ROIs
			data = Table.getColumn("Mean","Results");
			Table.setColumn(name1,data,"Cytosolic Intensity");
			selectWindow("Results");
			run("Close");
			close("*");
		}
	}
	selectWindow("Nuclear Intensity");
	Table.save(outputDir+"Nuclear Intensity.csv");
	run("Close");
	selectWindow("Cytosolic Intensity");
	Table.save(outputDir+"Cytosolic Intensity.csv");
	run("Close");
}
