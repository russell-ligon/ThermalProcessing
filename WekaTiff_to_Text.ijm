

// Macro imports a series of tiffs, each with 2 layers, the default output from the Weka probability classification
// it deletes the 2nd layer in each (which corresponds to the probability of a non-cold spot), keeping the first layer
// which corresponds to the probability that a given pixel is a cold pee spot
// Written by Russell Ligon, April 2020
//


  run("Close All");
 
  
  Parentdir = getDirectory("Choose trial directory containing probability TIFs");
  print(Parentdir);
  outdir = getDirectory("Choose outputfolder");
  
  
  fileList=getFileList(Parentdir);
  tifList=newArray();
 
 //This loops through contents of Parentdir
// and keeps only folders (which end with /)

for(i=0; i<fileList.length; i++){ // list only sub-folder names
 	if(endsWith(fileList[i], ".tif")==1)
		tifList = Array.concat(tifList, fileList[i]);		
}

		
		folderDIR2=Parentdir;
		folderDIR3= replace(folderDIR2, "\\", "/");
		folderDIR4=folderDIR3;
		//folderDIR4= replace(folderDIR3, "/", "\\");
		print(folderDIR4);
				
		outviddirOUT=outdir;
		outviddirOUT3= replace(outviddirOUT, "\\", "/");
		outviddirOUT4=outviddirOUT3;
		//outviddirOUT4= replace(outviddirOUT3, "/", "\\");

		 setBatchMode(true); 
 
for(j=0; j<tifList.length; j++){
	
	file = folderDIR4 + tifList[j];
	fileout = outviddirOUT4 + tifList[j];
	fileout = replace(fileout, ".tif", ".txt");
	
	open(file);
	setSlice(2);
	run("Delete Slice");
	saveAs("Text Image", fileout);
	close();
	
}
 
