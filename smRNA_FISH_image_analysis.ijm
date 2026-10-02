// this code is written by Anna Kogan at University of Oxford
// to analyse single-molecule RNA FISH data for Kogan et al., 2026

// It is partially based on the code from Aleksander Szczurek at University of Oxford
// published as ThunderFISH in Dobrinić, Paula et al. NSMB (2021)
// but uses 3D object-counter instead

// The input data structure is a folder with individual folders for each replicate or condition
// where the 4-channel images are moved into a folder named "in"

// The parameters of the 3D-object counter for both Ex and In data need to be optimised
// based on the labelling efficiency of the specific probes


// ask user to select a folder
dir = getDirectory("Select A folder");
// get the list of files (& folders) in it
fileList = getFileList(dir);


//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (w = 0; w < lengthOf(fileList); w++) {
	// define the "path" 
	// by concatenation of dir and the i element of the array fileList
	
	
	
	input = dir+fileList[w];
	// check that the currentFile is not a directory
	print(input);
IN = input + "in/";

//make directories for all the needed files
File.makeDirectory(input + "OutFISH_E");
File.makeDirectory(input + "OutFISH_I");
File.makeDirectory(input + "OutputMASK");

OutFISH_E = input + "OutFISH_E/"; 
OutFISH_I = input + "OutFISH_I/"; 
OutputMASK = input + "OutputMASK/"; 


/////////////////////////////////////////////////////////////////////////
/////////////////// 1. Split TF files into folders //////////////////////
/////////////////////////////////////////////////////////////////////////

function splitter(input, output1, output2, output3, filename) {
	      	while (nImages>0) {                                                	// close open images 
          	selectImage(nImages); 
          	close();       }

			old_Name = substring(filename,0,lengthOf(filename)-4);// deletes '.tif' from file name
      		print(old_Name);
      		open(input + filename);
      		title = getTitle();
       		number = substring(title,lengthOf(title)-5,lengthOf(title)-4);     //finds the number of the series (0-9 for 10 images)
			run("Split Channels");
			channel_4 = "C4-" + title;
			channel_3 = "C3-" + title;
			channel_2 = "C2-" + title;
			channel_1 = "C1-" + title;

			fileName = substring(title,0,lengthOf(title)-46);
			//print(fileName);
			newName = fileName + "_TF_" + number;
			fileName = newName;
			print(fileName);
       		

//make an mRNA substack
				selectWindow(channel_2);
				run("Grays");
					saveAs("Tiff", output1+fileName+"_FISH_E"+".tif"); 
				
//make a TS substack
				selectWindow(channel_1);
				run("Grays");
					saveAs("Tiff", output2+fileName+"_FISH_I"+".tif"); 

//make an Agglutinin mask
				selectWindow(channel_3);
				run("Grays");
					run("Z Project...", "projection=[Average Intensity]");
					run("Subtract Background...", "rolling=300 sliding");
					saveAs("Tiff", output3+fileName+"_GFP_Ave_MASK"+"_C0002.tif"); // don't change the ending!

}

setBatchMode(true); 
list = getFileList(IN);
print(list.length);
for (i = 0; i < list.length; i++)
        splitter(IN, OutFISH_E, OutFISH_I, OutputMASK, list[i]);
setBatchMode(false);
close();


////////////////////////////////////////////////////////////////////////////
/////////// 2. Put stacks and projections into separate subfolders /////////
////////////////////////////////////////////////////////////////////////////


current_FISH = OutFISH_E; ////////////////choose which FISH files to use

///  Put Ex stacks into separate subfolders

list = getFileList(current_FISH);

setBatchMode(true); 
for (n = 0; n < list.length; n++){  // n corresponds to image number
	if(endsWith(list[n],".tif")){  
		folderName = substring(list[n],0,lengthOf(list[n])-4); //cuts .tif
		File.makeDirectory(current_FISH + folderName); //until here it makes a folder per .tif file correctly
		outputPath = current_FISH + folderName;
		open(current_FISH + list[n]);
		imgName_Final=getTitle();
		print(imgName_Final);
		saveAs("Tiff",  outputPath +"\\" +list[n]);
	}
}

setBatchMode(false);

//delete .tiff files in outputFISH directory:
setBatchMode(true); 
	for (n = 0; n < list.length; n++){  
		if(endsWith(list[n],".tif")){
			File.delete(current_FISH + list[n]);                      // here deletes .tif files from the original dir as they're already in folder
		}}
setBatchMode(false); 
wait(300);

current_FISH = OutFISH_I; //////////////////choose which FISH files to use

///  Put Intr stacks into separate subfolders

list = getFileList(current_FISH);

setBatchMode(true); 
for (n = 0; n < list.length; n++){  // n corresponds to image number
	if(endsWith(list[n],".tif")){  
		folderName = substring(list[n],0,lengthOf(list[n])-4); //cuts .tif
		File.makeDirectory(current_FISH + folderName); //until here it makes a folder per .tif file correctly
		outputPath = current_FISH + folderName;
		open(current_FISH + list[n]);
		imgName_Final=getTitle();
		print(imgName_Final);
		saveAs("Tiff",  outputPath +"\\" +list[n]);
	}
}

setBatchMode(false);

//delete .tif files in outputFISH directory:
setBatchMode(true); 
	for (n = 0; n < list.length; n++){  
		if(endsWith(list[n],".tif")){
			File.delete(current_FISH + list[n]);                      // here deletes .tif files from the original dir as they're already in folder
		}}
setBatchMode(false); 
wait(300);


////////////////////////////////////////////////////////////////////////////
////////////// 3. Produce masks using 2D-background images /////////////////
////////////////////////////////////////////////////////////////////////////

close("*");

function masker(input, output, filename) {
	fileName = substring(filename,0,lengthOf(filename)-4); // deletes '.tif' from file name
        open(input + filename);
        run("Median...", "radius=6");                      //@user
		//run("Scale...", "x=0.5 y=0.5 width=600 height=600 interpolation=Bilinear average create");
		setAutoThreshold("Huang");
			setOption("BlackBackground", false);
			run("Convert to Mask");
			run("Invert");
			//run("Fill Holes");
			run("Watershed"); 
//			run("16-bit"); //optional, perhaps necessary
        	run("Invert LUT"); 
        	saveAs("Tiff", output + fileName + "_MASK");	
}

File.makeDirectory(OutputMASK + "masks");                    // make a subdirectory to save Masks there
outputFINAL = OutputMASK+"/masks/";                          // this is directory that the next step will be operating in!

setBatchMode(true); 
listMASK = getFileList(OutputMASK);   
for (i = 0; i < listMASK.length; i++)
	if(endsWith(listMASK[i],"C0002.tif")){
        masker(OutputMASK, outputFINAL, listMASK[i]);
}
setBatchMode(false);
close();
wait(300);
print("Progress 3/6");


/////////////////////////////////////////////////////////////////////////////////////////////////
/////// 4. Produce single cell masks using 2D-MASK images of full field of view /////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


if (isOpen("ROI Manager")) {
     selectWindow("ROI Manager");
     run("Close");
}

File.makeDirectory(outputFINAL+"singlecellmasks");     // make a subdir for single cell masks
outputFINALsc=outputFINAL+"/singlecellmasks/";         // variable storing name of that subdir

function SCcutter(inputPath, output, filename) {       // this is a definition of a function producing single masks
	open(inputPath + filename);
	setAutoThreshold("Default");
	run("Invert");
	run("Analyze Particles...", "size=6000-Infinity pixel pixel circularity=0.75-1.00 show=Masks exclude clear summarize add in_situ"); // @user - adjust your single cell stringecy!
	setBackgroundColor(0, 0, 0);
	print(roiManager("count"));
	for (i=0; i<roiManager("count"); i++) {
		roiManager("Select", i);
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		run("Clear Outside");
		saveAs("Tiff", output + "\\"+i+"_SingleCellMask.tif");    
	}	
} 

listSCmasks = getFileList(outputFINAL);     
setBatchMode(true); 
for (m = 0; m < listSCmasks.length; m++){  
	if(endsWith(listSCmasks[m],"_MASK.tif")){  
		newDir = outputFINALsc + listSCmasks[m];
		folderName = substring(listSCmasks[m],0,lengthOf(listSCmasks[m])-23); // cuts out the core name of image, i.e. without name extension
		File.makeDirectory(outputFINALsc+folderName);
		outputPath = outputFINALsc+folderName;
	    SCcutter(outputFINAL, outputPath, listSCmasks[m]); 
	    if (isOpen("ROI Manager")) {
        selectWindow("ROI Manager");
        run("Close");
}
	}
}
setBatchMode(false);

	      while (nImages>0) { // close any open images 
          selectImage(nImages); 
          close();       }
print("Progress 4/6");        
print("4) Now review your single cell masks and delete the ones that you don't like !");
print("    Manually remove accidentally merged cells or excessively cut cells !");
//waitForUser("Delete your poor-quality cells from \\\\outMASKS\\masks\\singlecellmasks NOW! Then press OK !"); 
//waitForUser("GREAT! You are ready for the second part  !"); 
wait(300);






//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////    second part   /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


ssMasks = OutputMASK + "masks/" + "singlecellmasks/";
current_FISH = OutFISH_E;                   //////////////////choose which FISH files to use

//////////////////// 1. DEFINE FUNCTION PRODUCING 3D ___EXON___ SINGLE CELL STACKS /////////
//////////////////////////////////////////////////////////////////////////////

end_TRITC=".tif";
end_Masks="_SingleCellMask.tif";

function multiplierEx(inputTRITC, inputMasks) {
setBatchMode(true);
listTRITC = getFileList(inputTRITC);
print(listTRITC.length);
listMaskFiles = getFileList(inputMasks);
print(listMaskFiles.length);
for (i = 0; i < listTRITC.length; i++){     
				if (endsWith(listTRITC[i], ".tif")){
					core = substring(listTRITC[i],0,lengthOf(listTRITC[i])-4);    // name ending: ".tif"
					print(listTRITC.length);
					open(inputTRITC + core + end_TRITC);
							run("Subtract Background...", "rolling=4 stack");
							run("Properties...", "channels=1 slices=60 frames=1 unit=px pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000"); // new 20-09-2020
							imgName_TRITC=getTitle();
							print(imgName_TRITC);
				}
}

for (k = 0; k < listMaskFiles.length; k++){         // k - single cell masks
				if (endsWith(listMaskFiles[k], "_SingleCellMask.tif")){
						core = substring(listMaskFiles[k],0,lengthOf(listMaskFiles[k])-19); 
   						open(inputMasks + core + end_Masks);
   						print(inputMasks + core + end_Masks);
						run("Select All");
						run("Invert");
						run("Multiply...", "value=0.0039"); // why?
						imgName_Mask=getTitle();  					

   					    imageCalculator("Multiply create stack", imgName_TRITC ,imgName_Mask);
   					imgName_Final=getTitle();		

                     //Cropper+Aligner below:
							run("Clear Results");
      					    title=getTitle;selectWindow(title);nz=nSlices;
      					    run("Set Measurements...", "center redirect=None decimal=4");
      					    run("Properties...", "channels=1 slices=1 frames=60 unit=px pixel_width=1 pixel_height=1 voxel_depth=1");
      					    run("Measure");
								x=getResult("XM",i-1); //translation vectors!
              				    y=getResult("YM",i-1);
               				    if(i==1){x0=x;y0=y;}
              				    sx=x0-x;sy=y0-y;
               				    run("Translate...", "x="+sx+" y="+sy+" slice");

							    run("Measure");  //measures center of gravity in aligned image
      								 x1=getResult("XM",i-1); //translation vectors!
      							     y1=getResult("YM",i-1);
         							 print(x1);
         							 print(y1);
        							 setTool("rectangle");
									 run("Specify...", "width=450 height=450 x=x1 y=y1 slice=10 constrain centered scaled"); 
									 run("Crop"); //crops the image to a single cell

               				    	 //selectWindow(title);
               				   		 saveAs("Tiff", inputMasks+core+"_smFISH_Ex");	
               				   		 print(inputMasks+core+"_smFISH_Ex");

		
				}
}
setBatchMode(false);
close();
} //end of function----------------------------------------------------------------------


mainTRITCList = getFileList(current_FISH);
print(mainTRITCList.length);
mainMasksList = getFileList(ssMasks);
print(mainMasksList.length);

setBatchMode(true);
for (l=0; l<mainTRITCList.length; l++) {  // for loop to parse through names in main folder
     if(endsWith(mainTRITCList[l], "/")){   // if the name is a subfolder...

          subDir = current_FISH + mainTRITCList[l]; //directory of l-folder, - one of the inputs of multiplier()
          subDirMasks = ssMasks + mainMasksList[l];
          multiplierEx(subDir, subDirMasks); // the actual function!

     }
}
setBatchMode(false);
wait(300);




ssMasks = OutputMASK + "masks/" + "singlecellmasks/";
current_FISH = OutFISH_E;                   //////////////////choose which FISH files to use

mainDirTRITC = current_FISH;
mainTRITCList = getFileList(current_FISH);
//print(mainTRITCList.length);
mainDirMasks = ssMasks;
mainMasksList = getFileList(ssMasks);
//print(mainMasksList.length);



//3. PERFORM 3D-object counting on single-cell stacks ___EXONS___//
////////////////////////////////////////////
setBatchMode(true);
for (b=0; b<mainMasksList.length; b++) { 
singleDir=mainDirMasks + mainMasksList[b];  //define a folder directory
listInDir = getFileList(singleDir);
print(listInDir.length);
for (d=0; d<listInDir.length; d++) { 
if (endsWith(listInDir[d], "Ex.tif")){
	print(listInDir[d]);
	open(singleDir+listInDir[d]); 
	cell=getTitle;
	getDimensions(width, height, channels, slices, frames);
	run("Properties...", "channels=1 slices="+frames+" frames="+slices+" unit=px pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
	//getStatistics(area, mean, min, max, std, histogram);
	Stack.getStatistics(count, mean, min, max, std);      // get statistics from the ENTIRE stack, not only the first frame!
	if (max < 160) threshold = max;                       // if threshold is the max of an image it wont pick any spots anyways!
  	else threshold = 160;
  	print(threshold);
  	    run("3D Objects Counter", "threshold=" + threshold +" slice=30 min.=10 max.=12150000 exclude_objects_on_edges statistics"); //dont exclude objects on edges!	min 30 better!
		g=b+1; //numerator of FOVs
		//saveAs("txt",input+i+"_frame_"+x+"min"); close("statistics"); //close(cell); // TXT
		saveAs("Results",singleDir+listInDir[d]+g+"_fov_" + d + "_cell_3D_mRNA_Spots_160_10.csv"); close("statistics"); //close(cell); // CSV
		
        
		//close immediately 3D-OC result (slow):
        	listWin = getList("window.titles"); //setBatchMode(true);
        	for (j=0; j<listWin.length; j++){
        	winame = listWin[j];
        	selectWindow(winame);
        	run("Close");
        	} //setBatchMode(false);

}
}
}
setBatchMode(false);

		   //close all 3D-OC results left:
        	listWin = getList("window.titles"); setBatchMode(true);
        	for (j=0; j<listWin.length; j++){
        	winame = listWin[j];
        	selectWindow(winame);
        	run("Close");
        	} setBatchMode(false);


current_FISH = OutFISH_I;                   //////////////////choose which FISH files to use

//////////////////// 1. DEFINE FUNCTION PRODUCING 3D ___INTRON___ SINGLE CELL STACKS /////////
//////////////////////////////////////////////////////////////////////////////

end_TRITC=".tif";
end_Masks="_SingleCellMask.tif";

function multiplierIn(inputTRITC, inputMasks) {
setBatchMode(true);
listTRITC = getFileList(inputTRITC);
print(listTRITC.length);
listMaskFiles = getFileList(inputMasks);
print(listMaskFiles.length);
for (i = 0; i < listTRITC.length; i++){     
				if (endsWith(listTRITC[i], ".tif")){
					core = substring(listTRITC[i],0,lengthOf(listTRITC[i])-4);    // name ending: ".tif"
					print(listTRITC.length);
					open(inputTRITC + core + end_TRITC);
							run("Subtract Background...", "rolling=4 stack");
							getDimensions(width, height, channels, slices, frames);
							run("Properties...", "channels=1 slices="+frames+" frames="+slices+" unit=px pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000"); // new 20-09-2020
							imgName_TRITC=getTitle();
							print(imgName_TRITC);
				}
}

for (k = 0; k < listMaskFiles.length; k++){         // k - single cell masks
				if (endsWith(listMaskFiles[k], "_SingleCellMask.tif")){
						core = substring(listMaskFiles[k],0,lengthOf(listMaskFiles[k])-19); 
   						open(inputMasks + core + end_Masks);
   						print(inputMasks + core + end_Masks);
						run("Select All");
						run("Invert");
						run("Multiply...", "value=0.0039"); // why?
						imgName_Mask=getTitle();  					

   					    imageCalculator("Multiply create stack", imgName_TRITC ,imgName_Mask);
   					imgName_Final=getTitle();		

                     //Cropper+Aligner below:
							run("Clear Results");
      					    title=getTitle;selectWindow(title);nz=nSlices;
      					    run("Set Measurements...", "center redirect=None decimal=4");
      					    getDimensions(width, height, channels, slices, frames);
      					    run("Properties...", "channels=1 slices="+frames+" frames="+slices+" unit=px pixel_width=1 pixel_height=1 voxel_depth=1");
      					    run("Measure");
								x=getResult("XM",i-1); //translation vectors!
              				    y=getResult("YM",i-1);
               				    if(i==1){x0=x;y0=y;}
              				    sx=x0-x;sy=y0-y;
               				    run("Translate...", "x="+sx+" y="+sy+" slice");

							    run("Measure");  //measures center of gravity in aligned image
      								 x1=getResult("XM",i-1); //translation vectors!
      							     y1=getResult("YM",i-1);
         							 print(x1);
         							 print(y1);
        							 setTool("rectangle");
									 run("Specify...", "width=450 height=450 x=x1 y=y1 slice=10 constrain centered scaled"); 
									 run("Crop"); //crops the image to a single cell

               				    	 //selectWindow(title);
               				   		 saveAs("Tiff", inputMasks+core+"_smFISH_In");	
               				   		 print(inputMasks+core+"_smFISH_In");

		
				}
}
setBatchMode(false);
close();
} //end of function----------------------------------------------------------------------



mainTRITCList = getFileList(current_FISH);
print(mainTRITCList.length);
mainMasksList = getFileList(ssMasks);
print(mainMasksList.length);

setBatchMode(true);
for (l=0; l<mainTRITCList.length; l++) {  // for loop to parse through names in main folder
     if(endsWith(mainTRITCList[l], "/")){   // if the name is a subfolder...

          subDir = current_FISH + mainTRITCList[l]; //directory of l-folder, - one of the inputs of multiplier()
          subDirMasks = ssMasks + mainMasksList[l];
          multiplierIn(subDir, subDirMasks); // the actual function!

     }
}
setBatchMode(false);
wait(300);

current_FISH = OutFISH_I;                   //////////////////choose which FISH files to use


mainDirTRITC = current_FISH;
mainTRITCList = getFileList(current_FISH);
//print(mainTRITCList.length);
mainDirMasks = ssMasks;
mainMasksList = getFileList(ssMasks);
//print(mainMasksList.length);



//3. PERFORM 3D-object counting on single-cell stacks ___INTRONS___//
////////////////////////////////////////////
setBatchMode(true);
for (b=0; b<mainMasksList.length; b++) { 
singleDir=mainDirMasks + mainMasksList[b];  //define a folder directory
listInDir = getFileList(singleDir);
print(listInDir.length);
for (d=0; d<listInDir.length; d++) { 
if (endsWith(listInDir[d], "In.tif")){
	print(listInDir[d]);
	open(singleDir+listInDir[d]); 
	cell=getTitle;
	getDimensions(width, height, channels, slices, frames);
	run("Properties...", "channels=1 slices="+frames+" frames="+slices+" unit=px pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
	//getStatistics(area, mean, min, max, std, histogram);
	Stack.getStatistics(count, mean, min, max, std);      // get statistics from the ENTIRE stack, not only the first frame!
	if (max < 500) threshold = max;                       // if threshold is the max of an image it wont pick any spots anyways!
  	else threshold = 500;
  	print(threshold);
  	    // below tested settings, min 20-150, treshold 105 - works robust, max 150px^3 size might not be enough for larger bursts, then use 250!
  		run("3D Objects Counter", "threshold=" + threshold +" slice=30 min.=25 max.=12150000 exclude_objects_on_edges statistics"); //dont exclude objects on edges!	min 30 better!
		g=b+1; //numerator of FOVs
		//saveAs("txt",input+i+"_frame_"+x+"min"); close("statistics"); //close(cell); // TXT
		saveAs("Results",singleDir+listInDir[d]+g+"_fov_" + d + "_cell_3D_TS_500_25.csv"); close("statistics"); //close(cell); // CSV
		
        
		//close immediately 3D-OC result (slow):
        	listWin = getList("window.titles"); //setBatchMode(true);
        	for (j=0; j<listWin.length; j++){
        	winame = listWin[j];
        	selectWindow(winame);
        	run("Close");
        	} //setBatchMode(false);

}
}
}
setBatchMode(false);

		   //close all 3D-OC results left:
        	listWin = getList("window.titles"); setBatchMode(true);
        	for (j=0; j<listWin.length; j++){
        	winame = listWin[j];
        	selectWindow(winame);
        	run("Close");
        	} setBatchMode(false); 
        	
        




}
setBatchMode(false);







