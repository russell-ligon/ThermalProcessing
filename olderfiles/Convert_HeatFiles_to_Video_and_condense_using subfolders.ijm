




// Macro imports a series of ; or , delimited csv text images as a stack, makes avi video
// and saves composite data where each row = one frame from the stack
// Written by Russell Ligon, February 26, 2018
//


  run("Close All");
  Trialdir = getDirectory("Choose trial directory containing subfolders, each containing semicolon delimited CSVs");
  print(Trialdir);
  outviddir = getDirectory("Choose thermal video parent folder");
  summdir = getDirectory("Choose summary csv parent folder");
  
  fileList=getFileList(Trialdir);
  folderList=newArray();
  
  
  
  summdir1= replace(summdir, "\\", "/");
  summdir2= replace(summdir1, "\\", "/");
  summdir3= replace(summdir2, "\\", "/");
  summdir4= replace(summdir3, "/", '\\');
  
 

//This loops through contents of Trialdir
// and keeps only folders (which end with /)
for(i=0; i<fileList.length; i++){ // list only sub-folder names
	//print(fileList[i]);
	//print(indexOf(fileList[i], "Standard"));
	//print(matches(fileList[i],".*tandard.*"));
	//print(endsWith(fileList[i], ".jpg")==1 & (indexOf(fileList[i], "Standard")<0));
	
	//Finds JPGs and jpgs w/o "tandard" in the name, adds them to tiffList
	/*
	if(endsWith(fileList[i], ".jpg")==1 && (matches(fileList[i],".*tandard.*")<1))
		tiffList = Array.concat(tiffList, fileList[i]);
	if(endsWith(fileList[i], ".JPG")==1 && (matches(fileList[i],".*tandard.*")<1))
		tiffList = Array.concat(tiffList, fileList[i]);
	if((matches(fileList[i],".*tandard.*")==1) && endsWith(fileList[i], ".JPG")==1)
		standardimage = fileList[i];	
	*/
	if(endsWith(fileList[i], "/")==1)
		folderList = Array.concat(folderList, fileList[i]);		
}

print("There are "+folderList.length+" folders to process");

 
	Dialog.create("CSV to AVI conversion settings");
	/*
		Dialog.addNumber("Future frames (how far into future to look for evap signature)?",240);
		Dialog.addNumber("If pixel is warmed in the future, how many frames back to jump back and check",10);
		Dialog.addNumber("Temperature above which is a mouse",79);
		Dialog.addNumber("Pee min temp",71);
		Dialog.addNumber("Pee evap temp threshold",67);
		Dialog.addNumber("Pixel distance for clustering",10);
		Dialog.addNumber("Minimum number of pixels to be considered a cluster",2);
	*/	
		Dialog.addChoice("Which computer are you working on?", newArray("megalaptop","Namaqua","sheehanlab1"));
		//Dialog.addNumber("Chunk size to focus on?",10000);
		Dialog.addChoice("What is the delimiter?", newArray("semicolon!","comma!"));
		Dialog.addNumber("How many composite output files have already been created?",0);
		Dialog.addChoice("Create AVIs?", newArray("true", "false"));
		Dialog.addChoice("Run AVI compiler too?", newArray("false", "true"));
		Dialog.addNumber("Number of frames per short avi?",2000);
		Dialog.addNumber("Frame rate for AVIs",4);
		Dialog.addChoice("Add timestamp",newArray("true","false"));
		Dialog.addChoice("BatchMode?", newArray("true", "false"));
		
		Dialog.addMessage("Do you want to delete original CSV files when program is done?");
		Dialog.addMessage("(original data IS saved in Summary_xxx csvs, if the program works correctly...)");
		Dialog.addChoice("Delete originals?", newArray("noway", "yes"));

	Dialog.show();
/*
	FutureFrame = Dialog.getNumber();
	JumpSize = Dialog.getNumber();
	MouseT = Dialog.getNumber();
	WarmUp = Dialog.getNumber();
	CoolOff = Dialog.getNumber();
	Distance = Dialog.getNumber();
	minpix = Dialog.getNumber();
*/	
	whichcomputer=Dialog.getChoice();
//	chunksize=Dialog.getNumber();
	runas= Dialog.getChoice();
	startwhere = Dialog.getNumber();
	MakeVid=Dialog.getChoice();
	VideoMaker= Dialog.getChoice();
	FrameNum = Dialog.getNumber();	
	FrameRate = Dialog.getNumber();	
	TimeStamp = Dialog.getChoice();
	BM = Dialog.getChoice();
	deleteEM = Dialog.getChoice();
	//scaleVal = Dialog.getNumber();
	//visualSystem = replace(visualSystem, "_", " ");
	//trials=Dialog.getNumber();
  
																								   
																											   
  if(TimeStamp=="true"){
	  Dialog.create("TimeStamp Settings");
		Dialog.addNumber("FontSize", 10);
		Dialog.addNumber("x coordinates", 5);
		Dialog.addNumber("y coordinates",470);
		Dialog.addChoice("FontColor",newArray("green","white","yellow","pink"));
	  Dialog.show();
	
		fontS= Dialog.getNumber();
		x1= Dialog.getNumber();
		y1 = Dialog.getNumber();	
		col = Dialog.getChoice();
		
  }
  

  /*
  csvList=newArray();

	for(i=0; i<list.length; i++) // list only csv files
		if(endsWith(list[i], ".csv")==1)
			csvList = Array.concat(csvList, list[i]);
*/		
  
  
  
  
// chunks = round(csvList.length/chunksize);
  

 
 
 /*
 startlist=newArray();
 endlist=newArray();
 
 for(t=0;t<csvList.length;){
	 startlist=Array.concat(startlist,t);
	 t=t+chunksize;
 }		
 endme=csvList.length;
  for(s=chunksize+1;s<csvList.length;){
	 endlist=Array.concat(endlist,s);
	 s=s+chunksize;
	 if(s>=csvList.length){
		endlist=Array.concat(endlist,endme);
	 }
 }
 
 */
//startlist=seq(1,csvList.length,chunksize);// #seq from 1 to Nfiles, increasing by HowBigMegas
//endlist=seq(chunksize,csvList.length,chunksize);// #seq from HowBigMegas, to Nfiles, increasing by HowBigMegas, and with the last file index added
 // 

  
 /*for(q=0;q<chunks+1;q++ ){ 
		showProgress(q, chunks);
		showStatus("Processing chunk "+q+" out of "+chunks);

		startem= startlist[q];
		endem = endlist[q];
		chunkList = Array.slice(csvList,startem,endem);
*/		

  if(BM=="true"){
	   setBatchMode(true);
  }
 


for(q=0; q<folderList.length; q++){ // cycle through every folder in folderList
  print("Processing folder "+(q+1)+" out of "+folderList.length);
  
		foldernombre=folderList[q];
		folderDIR2=Trialdir+foldernombre;
		folderDIR3= replace(folderDIR2, "\\", "/");
		folderDIR4= replace(folderDIR3, "/", "\\");
				print(folderDIR4);

				
				
		summdirOUT=summdir+foldernombre;
		summdirOUT3= replace(summdirOUT, "\\", "/");
		summdirOUT4= replace(summdirOUT3, "/", "\\");
		
		outviddirOUT=outviddir+foldernombre;
		outviddirOUT3= replace(outviddirOUT, "\\", "/");
		outviddirOUT4= replace(outviddirOUT3, "/", "\\");
		
		File.makeDirectory(summdirOUT4);
		File.makeDirectory(outviddirOUT4);
  
  list = getFileList(folderDIR4);
  csvList = list;


  start=0;
  stop=FrameNum;
//  hundogroups = round(chunkList.length/FrameNum);
  hundogroups = round(csvList.length/FrameNum);
  defaultvalue = 0;
  if(startwhere>0){
	defaultvalue=startwhere;
  }
//  print("Processing chunk "+q+" out of "+chunks);
  
  for(k=defaultvalue;k<hundogroups+1;k++ ) {
	print("Processing "+k+" out of "+hundogroups+"stacks");
		start=(k)*FrameNum;
		stop=(k+1)*FrameNum;
		openimages=0;
		
	//	 print("Processing stack "+k+" out of "+hundogroups+" for this chunk");

		
	  for (i=start; i<stop; i++) {
		if(i<csvList.length){ 
		file = folderDIR4 + csvList[i];
					
					if(runas=="comma!"){
					run("Text Image... ", "open=&file"); //use this for true comma separated csvs
					}
					if(runas=="semicolon!"){
					run("read semicolon", "open=&file");//read_semicolon is a custom java plugin written by R.A. Ligon, February 26, 2018
					}
					if(TimeStamp=="true"){		
						
						endloc=indexOf(file, ".");
						startloc=indexOf(file, "Record_2018");;
						timeinfo=substring(file, startloc,endloc);
						
						fontSize = fontS; 
						x = x1; 
						y = y1; 
						setColor(col); 
						setFont("SansSerif", fontSize); 
						
						//Overlay.remove; 
						Overlay.drawString(timeinfo, x, y); 
						Overlay.show; 
					}
			
		openimages=openimages+1;
		}
	  }
	  
	  if(start>=csvList.length){
		  
	  } else {
		 
		if(openimages>1){ 
		run("Images to Stack", "use");
		wait(1000);
		}
		//make sure odd extra,empty line didn't creep in
		makeRectangle(0, 0, 640, 480);
		run("Crop");
	
			nk=d2s(k, 0);	

			if(k<10){
			newk="0000"+nk;
			}
			if(k>9 && k<100){
			newk="000"+nk;
			}
			if(k>99 && k<1000){
			newk="00"+nk;
			}
			if(k>999 && k<10000){
			newk="0"+nk;
			}
	
	
	
	
	
		gargamel=start;
		//nameout="nameout= "+dir+gargamel+"_Summary_.csv";
		intermediate=summdirOUT4+'Summary_'+newk+'.csv';
		//nameout='nameout="'+summdir+'Summary_'+newk+'.csv"';
		nameout='nameout="'+intermediate+'"';
		
	if(whichcomputer=="Namaqua"){
		script = File.openAsString("C:\\Users\\Rusty\\Desktop\\Fiji.app\\plugins\\ThermalProcessing\\slice_rower.py"); 
	}
	if(whichcomputer=="sheehanlab1"){
		script = File.openAsString("C:\\Users\\Sheehan Lab\\Desktop\\Fiji.app\\plugins\\ThermalProcessing\\slice_rower.py"); 
	}
	if(whichcomputer=="megalaptop"){
		script = File.openAsString("C:\\Users\\Rusty Ligon\\Desktop\\Fiji.app\\plugins\\ThermalProcessing\\slice_rower.py"); 
	}	
			
		script = nameout+script;
		eval("python", script);
       
		
		
		//call("ij.plugin.Macro_Runner.runPython", script, nameout); 
		
		namesarray=newArray();
        for (n=1; n<=FrameNum; n++) {
			if(n<=nSlices){
			  setSlice(n);
			  namesarray = Array.concat(namesarray,getInfo("slice.label"));
			  setResult("Name", n-1, namesarray[n-1]);
			  updateResults();	
			}
		}

		setMinAndMax(65, 85);
		call("ij.ImagePlus.setDefault16bitRange", 8);
		run("Cool");
		
		
		if(MakeVid=="true"){
		run("AVI... ", "compression=JPEG frame=4 save=["+outviddirOUT4+"Stack_"+newk+".avi]");
		}
		selectWindow("Stack");
		close();
		
		saveAs("Results", summdirOUT4+"\\BatchNames_"+newk+".csv");
		
			if (isOpen("Results")) { 
				selectWindow("Results"); 
				run("Close"); 
			} 
	  }
    call("java.lang.System.gc");
	}


	if(deleteEM=="yes"){
		for(w=0;w<csvList.length;w++){
			file = folderDIR4 + csvList[w];
			File.delete(file);
		}
	}
  
	if(VideoMaker=="true"){ 
	list2 = getFileList(outviddirOUT4);
	  
	  for (i=0; i<list2.length; i++) {
		  list3 = getFileList(outviddirOUT4);
		  print(list3[0]);

			if(i<list2.length-1){
			print(list2[i+1]);
			
			file = outviddirOUT4 + list3[0];
				if (endsWith(file , ".avi")) {
						open(file);
						 }
			file2 = outviddirOUT4 + list2[i+1];
				if (endsWith(file2 , ".avi")) {
						open(file2);
						 }
			run("Concatenate...", "  title=[Concatenated Stacks] image1=["+list3[0]+"] image2=["+list2[i+1]+"] image3=[-- None --]");
			run("AVI... ", "compression=JPEG frame="+FrameRate+" save=["+outviddirOUT4+"CONCATENATED_VIDEO.avi]");
			selectWindow("Concatenated Stacks");
			close();
			}
		}
	}  
//  setBatchMode(false);

}
