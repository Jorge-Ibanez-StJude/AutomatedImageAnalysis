///setup
run("Set Measurements...", "area mean min centroid center bounding shape integrated kurtosis redirect=None decimal=3");
run("Colors...", "foreground=white background=black selection=yellow");
setOption("BlackBackground", false);
setOption("ScaleConversions", true);
////
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("Quantification_run_date:"+year+"/"+month+"/"+dayOfMonth+"["+hour+":"+minute+":"+second+"]");
macroname=File.getName(getInfo("macro.filepath"));
print(macroname);
///functions

////Modify functions in order to go 3D



///Mass fluorescence polarity analysis
function MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Channel){
	run("Clear Results");
	selectWindow("Tcellanalysis");
	run("Duplicate...", "title=MC.tif duplicate channels="+Channel);
	CellZ=nSlices;

	run("Z Project...", "projection=[Sum Slices]");
	getDimensions(preCellwidth, preCellheight, channels, slices, frames);///pixels
	run("Measure");
	x=getResult("XM", 0);///microns
	y=getResult("YM", 0);///microns	
	getPixelSize(unit,pw,ph);
	xMC=x/pw;///pixels
	yMC=y/ph;///pixels
	Cellwidth=preCellwidth;///pixels
	Cellheight=preCellheight;///pixels
	CH2=Cellheight*ph;///microns
	CW2=Cellwidth*pw;///microns	
	run("Clear Results");
	//XZ dimensions
	selectWindow("Tcellanalysis");
	run("Select None");
	run("Reslice [/]...", "output=0.280 start=Top avoid");
	run("Z Project...", "projection=[Sum Slices]");
	run("Measure");
	z=getResult("YM", 0);	
	getPixelSize(unit,pw,ph);
	zMC=z*ph;
	run("Clear Results");
	/// maths for calculate Mass center distance from the synapse
	///3Ddistance synap vs MC (x,y,z)
	MC_Synap_distance= sqrt((xMC-Cellwidth)*(xMC-Cellwidth)+(yMC-Cellheight/2)*(yMC-Cellheight/2)+(zMC-centslice*0.28)*(zMC-centslice*0.28));
	CellMass_Synap_distance= sqrt(synap_cell_d*synap_cell_d + (CellZ*0.28-centslice*0.28)*(CellZ*0.28-centslice*0.28));
	polarityIndex=((CellMass_Synap_distance-MC_Synap_distance)/CellMass_Synap_distance);
	////print(xMC+","+yMC+","+z+","+zMC+","+Cellwidth+","+Cellheight+","+CW2+","+CH2+","+MC_Synap_distance+","+CellMass_Synap_distance+","+polarityIndex);
	close("Reslice of Tcellanalysis");
	close("SUM_Reslice of Tcellanalysis");
	close("SUM_MC.tif");
	close("MC.tif");
	return polarityIndex;
}
///////////Heatmap analysis
function Heatmap(centslice,Channel,dir,name1){
	if (!File.exists(dir+File.separator+"Heatmap_"+Channel)) File.makeDirectory(dir+File.separator+"Heatmap_"+Channel);
	startslice=(centslice-5);
	stopslice=(centslice+5);
	selectWindow("Tcellanalysis");
	run("Duplicate...", "title=Heatmap.tif duplicate channels="+Channel);
	run("Z Project...", "projection=[Sum Slices]");
	run("8-bit");
	selectWindow("SUM_Heatmap.tif");
	run("Size...", "width=200 height=200 depth=35 average interpolation=Bilinear");
	///Cerrar y guardar la imagen de interes
	close("Heatmap.tif");
	selectWindow("SUM_Heatmap.tif");
	save(dir+File.separator+"Heatmap_"+Channel+File.separator+name1+"_Heatmap.tif");
	close("SUM_Heatmap.tif");
	return 1
}
/////////// synaptic area fluorescence analysis
function SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Channel){
	startslice=(centslice-5);
	stopslice=(centslice+5);
	selectWindow("Tcellanalysis");
	run("Duplicate...", "title=SAF.tif duplicate channels="+Channel);
	run("Z Project...", " projection=[Sum Slices]");
	selectWindow("SUM_SAF.tif");
	run("Flip Horizontally");
	selectWindow("SUM_SAF.tif");
 	run("Measure");
 	FT=getResult("RawIntDen");
 	run("Clear Results");
	makeRectangle(0, 0, (imagewd/3), imageht);
	run("Measure");
	FS=getResult("RawIntDen");
	run("Clear Results");
	PtjFSFT=(FS/FT)*100;
	close("SAF.tif");
	close("SUM_SAF.tif");
	run("Clear Results");
	return PtjFSFT;
}
////////Max MFI ration analysis
function MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Channel){
	startslice=(centslice-5);
	stopslice=(centslice+5);
	selectWindow(name);
	run("Duplicate...", "title=SAF.tif duplicate channels="+Channel);
	run("Z Project...", "projection=[Sum Slices]");
	selectWindow("SUM_SAF.tif");
 	run("Measure");
 	MaxF=getResult("Max",0);
 	MFI=getResult("Mean",0);
 	run("Clear Results");
	MaxMFIRes=(MaxF/MFI);
	close("SAF.tif");
	close("SUM_SAF.tif");
	run("Clear Results");
	return MaxMFIRes;
}
///////////////////synapsis analysis
function SpreadAnalisis(Channel){
	run("Duplicate...", "title=Spread.tif duplicate channels=Channel");
	run("Clear Results");	
	roiManager("select",0);
	run("Measure");
	PActLFluoTot=getResult("RawIntDen");
	PActLAreaTot=getResult("Area");
	PActLMaxFTot=getResult("Max");
	PActLMfiTot=getResult("Mean");
	AspectRatio=getResult("AR", 0);
	circularity=getResult("Circ.", 0);
	Roundness=getResult("Round", 0);
	run("Clear Results");
	selectWindow("Spread.tif");			
	roiManager("select",1);
	run("Measure");
	PActLFluoCent=getResult("RawIntDen");
	PActLAreaCent=getResult("Area");
	PActLMaxFCent=getResult("Max");
	PActLMfiCent=getResult("Mean");
	PActLFluoPer=PActLFluoTot-PActLFluoCent;
	PActLpArea=PActLAreaCent/PActLAreaTot*100;
	PActLpFluoCent=PActLFluoCent/PActLFluoTot*100;
	PActLpFluoPer=PActLFluoPer/PActLFluoTot*100;
	ActdenFluoCent=(PActLpFluoCent/PActLpArea-1);
	ActdenFluoPer=(PActLpFluoPer/(100-PActLpArea)-1);
	FDisCent=PActLMaxFCent/PActLMfiCent;
	FDisTot=PActLMaxFTot/PActLMfiTot;

	IntdenCenter=PActLFluoCent/PActLAreaCent;
	IntdenTotal=PActLFluoTot/PActLAreaTot;
	
	print("Immune synapse,"+image+","+Channel+","+name1+","+PActLAreaTot+","+AspectRatio+","+circularity+","+Roundness+","+PActLpFluoCent+","+PActLpFluoPer+","+ActdenFluoCent+","+ActdenFluoPer+","+FDisCent+","+FDisTot+","+PActLMfiTot+","+IntdenCenter+","+IntdenTotal);
	run("Clear Results");	
	close("Spread.tif");
	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///get dir of images and open folderdir
setOption("JFileChooser", true);
dir= getDirectory("Choose the folder that contains the images");
imagenames=getFileList(dir); /// image dir
nbimages=lengthOf(imagenames); /// image names

run("Bio-Formats Importer", "open=["+dir+File.separator+imagenames[0]+"] color_mode=Default view=Hyperstack stack_order=XYCZT");///to import any sort of  file
///to increase the size of the images if is too small
getDimensions(imagew,imageh,channels,slices,frames);
if(imagew<2000){
	Stack.setDisplayMode("composite");
	run("In [+]");
	
	getDimensions(width,height,channels,slices,frames);
	for (d = 1; d < channels; d++) {
	Stack.setChannel(d);
	resetMinAndMax();	
	}
	
}
waitForUser("Check the image and its channels");
///user interface to determine what to do next
Dialog.create("What do you want to do?");
Dialog.addCheckbox("manual crop?", false);
Dialog.addCheckbox("automatic IS crop?", false);
Dialog.addCheckbox("select IS interface? (related to IS crop)", false);
Dialog.addCheckbox("resume a previous analysis?",false);
Dialog.addCheckbox("quantify neurospheres(automatic)?", false);
Dialog.addCheckbox("quantify cocultures (automatic)?", false);
Dialog.show();
///getting the answers
Cut = Dialog.getCheckbox();
CutIS = Dialog.getCheckbox();
selectZ = Dialog.getCheckbox();
resumeanalysis = Dialog.getCheckbox();
Neurospheres= Dialog.getCheckbox();
AutoCoculture= Dialog.getCheckbox();
File.makeDirectory(dir+"Results");
close(imagenames[0]);
///loop to crop images and save them into a new folder
if(AutoCoculture== true){
	///channel discrimination
	Dialog.addNumber("Which channel is T cells?",0);
	Dialog.addNumber("Which channel is Tumor cells?",0);
	Dialog.show();
	TcellChannel = Dialog.getNumber();
	TumorChannel = Dialog.getNumber();
	///Tcell segmentation
	setOption("JFileChooser", true);
	dir_classifier= getDirectory("Choose the classifier folder");
	setOption("JFileChooser", false);
	model_name ="";
	Dialog.create("type the T cell model name plus its extension .arff o .model");
	Dialog.addString(model_name, "model name");
	Dialog.show();
	Tcell_model_name = Dialog.getString();
	Tcell_classifier = dir_classifier+Tcell_model_name;
	///Tumorcell segmentation
	Dialog.create("type the Tumot cell model name plus its extension .arff o .model");
	Dialog.addString(model_name, "model name");
	Dialog.show();
	Tumor_model_name = Dialog.getString();
	Tumor_classifier = dir_classifier+Tumor_model_name;
}




if(Cut == true){
	File.makeDirectory(dir+File.separator+"Cropped");			
	for(image=0; image<nbimages; image++) {
		name=imagenames[image];
		totnamelength=lengthOf(name); /// extensión del nombre
		namelength=totnamelength-4;
		name1=substring(name, 0, namelength);
		extension=substring(name, namelength, totnamelength);
			open(dir+File.separator+name);		
			getDimensions(width,height,channels,slices,frames);
			getPixelSize(unit,pw,ph);
			makeRectangle(0, 0, 20/pw, 20/ph);
			waitForUser ("Select the cells to crop (ctrl+t o solo t)");
			selectWindow(name);
			wait(100);
			numROI=roiManager("count"); 
			for (i=0; i<numROI; i++){
				selectWindow(name);
				roiManager("Select", i);
				run("Duplicate...", "title=["+name+"_"+i+"] duplicate channels=1-"+channels+" slices=1-"+slices);

			}
			selectWindow(name);
			run("Close");
			wait(100);
			openimages=nImages;
			
			for (ROIs=0; ROIs<openimages ; ROIs++) {
			    cutname=getTitle();
			    saveAs("Tiff", dir+File.separator+"Cropped"+File.separator+cutname);    
			    run("Close");
			}
				ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}
		}
	}

if(CutIS == true){
	File.makeDirectory(dir+File.separator+"Cropped");
	Dialog.addNumber("Which channel is F-actin (Phalloidin)?",0);
	Dialog.show();
	ActinChannel = Dialog.getNumber();
	
	setOption("JFileChooser", true);
	dir_classifier= getDirectory("Choose the classifier folder");
	setOption("JFileChooser", false);
	model_name ="";
	Dialog.create("type the model name plus its extension .arff o .model");
	Dialog.addString(model_name, "model name");
	Dialog.show();
	model_name = Dialog.getString();
	classifier = dir_classifier+model_name;
		
	for(image=0; image<nbimages; image++) {
		name=imagenames[image];
		totnamelength=lengthOf(name); /// extensión del nombre
		namelength=totnamelength-4;
		name1=substring(name, 0, namelength);
		extension=substring(name, namelength, totnamelength);
			///setOption("BlackBackground", false);
			open(dir+File.separator+name);		
			if(selectZ== true){
				waitForUser("select the Zstack of the IS");
				run("Duplicate...", "title=Synapticstack");
				wait(100);
				close(name);
			}
			selectWindow("Synapticstack");
			getDimensions(width,height,channels,slices,frames);
			getPixelSize(unit,pw,ph);
			run("Duplicate...", "title=WEKA.tif duplicate channels=ActinChannel");
			run("Subtract Background...", "rolling=20");
			run("Advanced Weka Segmentation");
			///call("trainableSegmentation.Weka_Segmentation.setFeature", "Mean=true");
			///call("trainableSegmentation.Weka_Segmentation.setFeature", "Median=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Variance=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Mean=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Minimum=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Maximum=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Median=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Structure=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Entropy=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Neighbors=true");
			call("trainableSegmentation.Weka_Segmentation.setClassBalance", "true");
			wait(500);
			call("trainableSegmentation.Weka_Segmentation.loadClassifier", classifier);
			call("trainableSegmentation.Weka_Segmentation.changeClassName", "0", "Background");
			call("trainableSegmentation.Weka_Segmentation.changeClassName", "1", "Cell");
			wait(500);
			call("trainableSegmentation.Weka_Segmentation.getProbability");
			wait(500);
			selectWindow("Probability maps");
			wait(500);
			run("Duplicate...", "title=Cell duplicate range=2 use");
			close("Probability maps");
			selectWindow("Cell");
			wait(500);
			run("8-bit");
			///run("Measure");
			///Cell_Background_ratio=getResult("Mean");
			///run("Clear Results");
			///if(Cell_Background_ratio>125){
				///setOption("BlackBackground", false);
			///	}
			///else {
			///	setOption("BlackBackground", true);
		///	}		
			run("Make Binary");
			///setOption("BlackBackground", true);
			run("Erode");
			run("Dilate");
			run("Fill Holes");
			run("Watershed");
			wait(200);
			selectWindow("Cell");
			run("Analyze Particles...", "size=30-Infinity exclude add");
			wait(500);
			waitForUser ("check (use paintbrush if its needed to delimit cells)");

			///close("Background");
			///close("Probability maps");
			selectWindow("Trainable Weka Segmentation v3.2.35");
			wait(100);
			run("Close");
			close("WEKA.tif");
			close("Cell");
			selectWindow("Synapticstack");
			wait(100);
			numROI=roiManager("count"); 
			for (i=0; i<numROI; i++){
				selectWindow("Synapticstack");
				roiManager("Select", i);
				run("Duplicate...", "title=["+name+"_"+i+"] duplicate channels=1-"+channels+" slices=1-"+slices);
				run("Clear Outside");
			}
			selectWindow("Synapticstack");
			run("Close");
			wait(100);
			openimages=nImages;
			
			for (ROIs=0; ROIs<openimages ; ROIs++) {
			    cutname=getTitle();
			    saveAs("Tiff", dir+File.separator+"Cropped"+File.separator+cutname);    
			    run("Close");
			}
				ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}
		}
		selectWindow("Log");
		saveAs("text", dir+File.separator+"Results"+File.separator+"Segmentation_Log.csv");
		close("Log");
	}

if (Neurospheres == true) {
	setOption("JFileChooser", true);
	dir_classifier= getDirectory("Choose the classifier folder");
	setOption("JFileChooser", false);
	model_name ="";
	Dialog.create("type the model name plus its extension .arff o .model");
	Dialog.addString(model_name, "model name");
	Dialog.show();
	model_name = Dialog.getString();
	classifier = dir_classifier+model_name;
		
	for(image=0; image<nbimages; image++) {
		name=imagenames[image];
		totnamelength=lengthOf(name); /// extensión del nombre
		namelength=totnamelength-4;
		name1=substring(name, 0, namelength);
		extension=substring(name, namelength, totnamelength);
			open(dir+File.separator+name);	
			run("Split Channels");
			close(name+" (green)");	
			close(name+" (red)");
			selectWindow(name+" (blue)");
			getDimensions(width,height,channels,slices,frames);
			getPixelSize(unit,pw,ph);
			run("Advanced Weka Segmentation");
			///call("trainableSegmentation.Weka_Segmentation.setFeature", "Mean=true");
			///call("trainableSegmentation.Weka_Segmentation.setFeature", "Median=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Entropy=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Neighbors=true");
			call("trainableSegmentation.Weka_Segmentation.setFeature", "Structure=true");
			call("trainableSegmentation.Weka_Segmentation.setClassBalance", "true");
			wait(500);
			call("trainableSegmentation.Weka_Segmentation.loadClassifier", classifier);
			call("trainableSegmentation.Weka_Segmentation.changeClassName", "0", "Background");
			call("trainableSegmentation.Weka_Segmentation.changeClassName", "1", "Neurospheres");
			wait(500);
			call("trainableSegmentation.Weka_Segmentation.getProbability");
			wait(500);
			selectWindow("Probability maps");
			run("Duplicate...", "title=Neurospheres duplicate range=2 use");
			selectWindow("Neurospheres");
			run("8-bit");
			///run("Measure");
			///Cell_Background_ratio=getResult("Mean");
			///run("Clear Results");
			///if(Cell_Background_ratio>125){
				///setOption("BlackBackground", false);
			///	}
			///else {
			///	setOption("BlackBackground", true);
		///	}		
			run("Make Binary");
			///0run("Erode");
			///run("Dilate");
			///run("Fill Holes");
			///run("Watershed");
			wait(200);
			run("Erode");
			selectWindow("Neurospheres");
			run("Analyze Particles...", "size=800-Infinity exclude add");
			wait(500);
			close("Background");
			close("Probability maps");
			selectWindow("Trainable Weka Segmentation v3.2.34");
			wait(100);
			run("Close");
			wait(100);
			numROI=roiManager("count");
			print("analysis_type,name, number of neurospheres, Neurosphere number, Area, Circularity"); 
			for (i=0; i<numROI; i++){
				selectWindow(name+" (blue)");
				roiManager("Select", i);
				run("Measure");
				NeuArea=getResult("Area");
				Circularity=getResult("Circ.");
				print("Neurosphere_quant,"+name+","+numROI+","+i+","+NeuArea+","+Circularity);
				run("Clear Results");
			}
			close("WEKA.tif");
			close("Neurospheres");
			close(name+" (blue)");
			wait(100);
			ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}

	}

	selectWindow("Log");
	saveAs("text", dir+File.separator+"Results"+File.separator+name+".csv");
	close("Log");
}



///Loop to resume a previous analysis, if its a new analysis start will be 0
start = 0;
dir_Cropped= dir+File.separator+"Cropped"+File.separator;
if(resumeanalysis == true){
	Dialog.create("indicate the number of the image to resume");
	Dialog.addNumber("add number", 0);
	Dialog.show();
	restart = Dialog.getNumber();
}

imagenames=getFileList(dir_Cropped); /// image dir
nbimages=lengthOf(imagenames); /// image names
name=imagenames[0];
totnamelength=lengthOf(name); /// extensión del nombre
namelength=totnamelength-4;
name1=substring(name, 0, namelength);
extension=substring(name, namelength, totnamelength);
		
	run("Bio-Formats Importer", "open=["+dir_Cropped+imagenames[0]+"] color_mode=Default view=Hyperstack stack_order=XYCZT");

	if(imagew<2000){
	run("In [+]");
	run("In [+]");
	run("In [+]");
	}
///user interface to check channels ???check
	waitForUser("Check every channel");
	Dialog.create("select the correct color for each channel");
	Dialog.addMessage("select the number corresponding to each color, if the image doesnt have that color just let none");
	Dialog.addChoice("How many channels have yout image?", newArray("1","2","3","4","5"),0);
	Dialog.addChoice("Green:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.addChoice("Red:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.addChoice("Blue:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.addChoice("Widefield:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.addChoice("Cyan:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.addChoice("Magenta:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.addChoice("Yellow:",newArray("*None*","1","2","3","4","5"),"*None*");
	Dialog.show();
	
	NChannels = Dialog.getChoice();
	Green = Dialog.getChoice();
	Red = Dialog.getChoice();
	Blue = Dialog.getChoice();
	Widefield = Dialog.getChoice();	
	Cyan = Dialog.getChoice();
	Magenta = Dialog.getChoice();
	Yellow = Dialog.getChoice();


	/// Analysis selection
	
	EachChannel = newArray(Green,Red,Blue,Cyan,Magenta,Yellow);
	ChannelColor = newArray("Green","Red","Blue","Cyan","Magenta","Yellow");
	/// analysis interface
		rows = 3;
		columns = 6;
		Analisis = newArray("Mass Center Polarity","Heatmap","SinapticFluo","FMax/MFI","Cell Shape","Immune Synapse");
		defaults = newArray(false,false,false,false,false,false);	
		Dialog.create("Analysis selection");
		Dialog.addMessage("Green Channel --------------------------  Green Channel ------------------------");
		Dialog.addCheckboxGroup(rows,columns,Analisis,defaults);
		Dialog.addMessage("Red Channel ---------------------------- Red Channel -------------------------");
		Dialog.addCheckboxGroup(rows,columns,Analisis,defaults);
		Dialog.addMessage("Blue Channel ---------------------------- Blue Channel -------------------------");
		Dialog.addCheckboxGroup(rows,columns,Analisis,defaults);
		Dialog.addMessage("Cyan Channel---------------------------- Cyan Channel-------------------------");
		Dialog.addCheckboxGroup(rows,columns,Analisis,defaults);
		Dialog.addMessage("Magenta Channel------------------------- Magenta Channel----------------------");
		Dialog.addCheckboxGroup(rows,columns,Analisis,defaults);
		Dialog.addMessage("Yellow Channel ------------------------ Yellow Channel ---------------------");
		Dialog.addCheckboxGroup(rows,columns,Analisis,defaults);
		Dialog.show();
		
		//C1
		a1 = Dialog.getCheckbox();///Mass center polarity
		b1 = Dialog.getCheckbox();///Heatmap
		c1 = Dialog.getCheckbox();///Synaptic fluorescence
		d1 = Dialog.getCheckbox();///Max fluo and MFI ratio
		e1 = Dialog.getCheckbox();///Cell shape
		f1 = Dialog.getCheckbox();///Immune synapse
		//C2
		a2 = Dialog.getCheckbox();///Mass center polarity
		b2 = Dialog.getCheckbox();///Heatmap
		c2 = Dialog.getCheckbox();///Synaptic fluorescence
		d2 = Dialog.getCheckbox();///Max fluo and MFI ratio
		e2 = Dialog.getCheckbox();///Cell shape
		f2 = Dialog.getCheckbox();///Immune synapse
		//C3
		a3 = Dialog.getCheckbox();///Mass center polarity
		b3 = Dialog.getCheckbox();///Heatmap
		c3 = Dialog.getCheckbox();///Synaptic fluorescence
		d3 = Dialog.getCheckbox();///Max fluo and MFI ratio
		e3 = Dialog.getCheckbox();///Cell shape
		f3 = Dialog.getCheckbox();///Immune synapse
		//C4
		a4 = Dialog.getCheckbox();///Mass center polarity
		b4 = Dialog.getCheckbox();///Heatmap
		c4 = Dialog.getCheckbox();///Synaptic fluorescence
		d4 = Dialog.getCheckbox();///Max fluo and MFI ratio
		e4 = Dialog.getCheckbox();///Cell shape
		f4 = Dialog.getCheckbox();///Immune synapse
		//C5
		a5 = Dialog.getCheckbox();///Mass center polarity
		b5 = Dialog.getCheckbox();///Heatmap
		c5 = Dialog.getCheckbox();///Synaptic fluorescence
		d5 = Dialog.getCheckbox();///Max fluo and MFI ratio
		e5 = Dialog.getCheckbox();///Cell shape
		f5 = Dialog.getCheckbox();///Immune synapse
		//C6
		a6 = Dialog.getCheckbox();///Mass center polarity
		b6 = Dialog.getCheckbox();///Heatmap
		c6 = Dialog.getCheckbox();///Synaptic fluorescence
		d6 = Dialog.getCheckbox();///Max fluo and MFI ratio
		e6 = Dialog.getCheckbox();///Cell shape
		f6 = Dialog.getCheckbox();///Immune synapse


while (nImages>0) { 
	selectImage(nImages); 
	close(); 
} 

///setBatchMode("show");
		if(f1||f2||f3||f4||f5||f6 == true){
			print("Analysis,Number,Analyzed Channel,Name, Area, Aspect-Ratio, Circularity, Roundness, %Fluo-center/total, %Fluo-Periphery/total, denFluo-Center/total, denFluo-Periphery/Total, Centro Max Fluo/MFI, Total Max Fluo/MFI, total MFI, IntDenCenter, IntDenTotal");
			setTool("polygon");
			ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}
		Dialog.create("For immune synapse analysis");
		Dialog.addNumber("Which channel is F-actin(Phalloidin)?",0);
		Dialog.show();
		ActinChannel = Dialog.getNumber();
		}
		


imagenames=getFileList(dir_Cropped); /// Cell images directory
nbimages=lengthOf(imagenames); /// images name

limit= 10000;


///to define in which cell to begin the analysis (when is aanalysis resume)
if(resumeanalysis ==true){
	start=restart;
}
else{
	start=0;
}


/// Loop for image analysis
for(image=start; image<nbimages; image++) { 
	
	name=imagenames[image];
	totnamelength=lengthOf(name); /// name extension
	namelength=totnamelength-4;
	name1=substring(name, 0, namelength);
	extension=substring(name, namelength, totnamelength);
	dirFoto = dir_Cropped+name;
	
	if(extension==".tif" || extension==".tiff" || extension==".nd2") {
		run("Bio-Formats Importer", "open=["+dirFoto+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
		if(imagew<2000){
			run("In [+]");
			run("In [+]");
			run("In [+]");
		run("Properties...", "pixel_width=0.1100000 pixel_height=0.1100000 voxel_depth=0.28");
		getDimensions(width,height,channels,slices,frames);
		Stack.setDisplayMode("composite");
		imageH=height;
		imageW=width;
		}
		
		if(f1||f2||f3||f4||f5||f6 == true){
			
		selectWindow(name1+".tif");
		run("Duplicate...", "title=Actin.tif duplicate channels=ActinChannel");
		///waitForUser("check");
		run("8-bit");
		///setAutoThreshold("Mean dark");
		///setOption("BlackBackground", false);
		///run("Convert to Mask");
		setAutoThreshold("Huang dark");
		run("Convert to Mask");
		run("Dilate");
		run("Close-");
		run("Fill Holes");	
		///run("Outline");
		///run("Clear Results");
		////manual selection
		/////waitForUser("Delimit the cell");
		////roiManager("Add");
		///automatic selection
		run("Analyze Particles...", "size=0.6-Infinity add");
		///roiManager("Add");
		roiManager("Select",0);
		run("Measure");
		actprex=getResult("X", 0);
		actprey=getResult("Y", 0);
		actprew=getResult("Width",0);
		actpreh=getResult("Height",0);
		getPixelSize(unit,pw,ph);
		actx=actprex/pw;
		acty=actprey/ph;
		actw=actprew/pw;
		acth=actpreh/ph;
		makeOval(actx-actw/4,acty-acth/4,actw/2,acth/2);
		///waitForUser("check");
		roiManager("Add");
		roiManager("Select",1);
		run("Clear Results");
		close("Actin.tif");
		

		}





		if(a1||a2||a3||a4||a5||a6||b1||b2||b3||b4||b5||b6||c1||c2||c3||c4||c5||c6||d1||d2||d3||d4||d5||d6 == true){

			if(AutoCoculture== true){
				ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}
				///T cell cropping
				selectWindow(name1+".tif");
				Tcc=TcellChannel;
				run("Duplicate...", "title=TcellWEKA duplicate channels=Tcc");
				selectWindow("TcellWEKA");
				run("Trainable Weka Segmentation 3D");
				call("trainableSegmentation.Weka_Segmentation.setFeature", "Mean=true");
				call("trainableSegmentation.Weka_Segmentation.setFeature", "Median=true");
				call("trainableSegmentation.Weka_Segmentation.setFeature", "Entropy=true");
				call("trainableSegmentation.Weka_Segmentation.setFeature", "Variance=true");
				call("trainableSegmentation.Weka_Segmentation.setFeature", "Structure=true");
				call("trainableSegmentation.Weka_Segmentation.setFeature", "Minimum=true");
				call("trainableSegmentation.Weka_Segmentation.setClassBalance", "true");
				wait(500);
				call("trainableSegmentation.Weka_Segmentation.loadClassifier",Tcell_classifier);
				call("trainableSegmentation.Weka_Segmentation.changeClassName", "0", "Tcell");
				call("trainableSegmentation.Weka_Segmentation.changeClassName", "1", "Background");
				wait(500);
				call("trainableSegmentation.Weka_Segmentation.getProbability");
				wait(500);
				selectWindow("Probability maps");
				run("Duplicate...", "title=Tcell duplicate channels=1");
				selectWindow("Tcell");
				getDimensions(width,height,channels,slices,frames);
				midslice=slices/2;
				setSlice(1);
				run("8-bit");
				setAutoThreshold("Minimum dark");
				run("Threshold...");
				run("Convert to Mask", "method=Minimum background=Dark");
				run("Fill Holes");
				run("Dilate", "stack");
				run("Dilate", "stack");
				waitForUser("Check the Tcell segmentation, if it is something odd, use paintbrush tool to delete or modify undesired segmentations");
				wait(500);
				run("Analyze Particles...", "size=10-Infinity exclude add stack");			
				numROI=roiManager("count");
				roiManager("Select",0);
				slice0=getSliceNumber();
				finalslice=numROI-1;
				roiManager("Select",finalslice);
				sliceN=getSliceNumber();
				close("C1");
				selectWindow("Trainable Weka Segmentation v3.2.35");
				close();
				close("Probability maps");
				close("TcellWEKA");
				close("Threshold");
				for (z = 0; z < NChannels; z++) {
					channelactive=z+1;
					wait(100);
					selectWindow(name1+".tif");
					run("Duplicate...", "title=C"+z+" duplicate channels=channelactive slices=slice0-sliceN");
					wait(10);
					selectWindow("C"+z);
					roiManager("Select", 0);
					run("Clear Outside", "slice");
					for (i=0; i<numROI; i++){
						wait(10);
						selectWindow("C"+z);
						roiManager("Select", i);
						run("Clear Outside", "slice");
					}		
				}
				///Merge cropped channels
				if(NChannels==4){
					run("Merge Channels...", "c1=[C0] c2=[C1] c3=[C2] c4=[C3] create");
					selectWindow("Composite");
					rename("Tcellanalysis");	
				}
				if(NChannels==3){
					run("Merge Channels...", "c1=[C0] c2=[C1] c3=[C2] create");
					selectWindow("Composite");
					rename("Tcellanalysis");	
				}


			ROI=isOpen("ROI Manager");
			if(ROI == true){
			selectWindow("ROI Manager");
			run("Close");
			}
			selectWindow(name1+".tif");
			run("Duplicate...", "title=SynapCoordinates duplicate slices=slice0-sliceN");
			setTool("point");
			waitForUser("indicate T cell tumor cell interaction point (it has to be at the middle of the interaction (XYZ wise))");
			Stack.getPosition(channel, slice, frame);
			centslice=slice;
			getPixelSize(unit,pw,ph);
			run("Measure");
			precxsynap=getResult("X", 0);
			precysynap=getResult("Y", 0);
			wait(10);
			run("Clear Results");

			falsesynapse=(width*pw)/2;
		
			if(precxsynap == falsesynapse){
			waitForUser("indicate T cell tumor cell interaction point (second shot)");
			Stack.getPosition(channel, slice, frame);
			centslice=slice;
			getPixelSize(unit,pw,ph);
			run("Measure");
			precxsynap=getResult("X", 0);
			precysynap=getResult("Y", 0);
			wait(10);			
			}
			
			run("Clear Results");
			if(precxsynap == falsesynapse){
			waitForUser("indicate T cell tumor cell interaction point (third and last shot)");
			Stack.getPosition(channel, slice, frame);
			centslice=slice;
			getPixelSize(unit,pw,ph);
			run("Measure");
			precxsynap=getResult("X", 0);///micron dimension
			precysynap=getResult("Y", 0);///micron dimension
			wait(10);			
			}
			
			///synapse coordinates transform from microns to pixels
			xsynap = (precxsynap/pw);///pixel dimension
			ysynap = (precysynap/ph);	///pixel dimension
			zsynap = centslice*0.28;///micron dimension
			run("Clear Results");
			
			///T cell coordinates
			selectWindow("Tcellanalysis");
			run("Duplicate...", "title=Tcellcoordinates duplicate channels=Tcc");
			run("Z Project...", "projection=[Average Intensity]");
			///XY coordinates
			selectWindow("AVG_Tcellcoordinates");
			run("8-bit");
			setAutoThreshold("Minimum dark");
			run("Convert to Mask");
			run("Analyze Particles...", "size=5-Infinity add");
			roiManager("Select", 0);
			run("Measure");
			getPixelSize(unit,pw,ph);		
			getSelectionBounds(xCell,yCell,wCell,hCell);
			precxCell=getResult("X", 0); /// x coordinate of Cell mass center micron dimension
			precyCell=getResult("Y", 0);/// y coordinate of cell mass center micron dimension
			wait(10);
			///get cell shape parameters
			if(e1||e2||e3||e4||e5||e6 ==true){
				selectWindow("AVG_Tcellcoordinates");
				AspectRatioXY=getResult("AR", 0);
				circularityXY=getResult("Circ.", 0);
				RoundnessXY=getResult("Round", 0);
				AreaXY=getResult("Area", 0);
			}
			else{
				AspectRatioXY=0;
				circularityXY=0;
				RoundnessXY=0;
				AreaXY=0;
			}
			run("Clear Results");
			ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
			}
			close("AVG_Tcellcoordinates");
			wait(100);
			///Z coordinates
			selectWindow("Tcellcoordinates");
			run("Select None");
			run("Reslice [/]...", "output=0.280 start=Top avoid");
			wait(10);

			selectWindow("Reslice of Tcellcoordinates");
			run("Z Project...", "projection=[Average Intensity]");
			wait(10);

			selectWindow("AVG_Reslice of Tcellcoordinates");
			run("8-bit");
			setAutoThreshold("Minimum dark");
			run("Convert to Mask", "method=Minimum background=Dark");
			wait(10);
			run("Analyze Particles...", "size=10-Infinity add");
			roiManager("Select", 0);
			run("Measure");
			prexResCell=getResult("X", 0);/// in micron dimensions
			preczCell=getResult("Y", 0);/// in micron dimensions
			
			close("Tcellcoordinates");
			close("Reslice of Tcellcoordinates");
			///get cell shape parameters
			if(e1||e2||e3||e4||e5||e6 ==true){	
				run("Measure");
				selectWindow("AVG_Reslice of Tcellcoordinates");
				AspectRatioXZ=getResult("AR", 0);
				circularityXZ=getResult("Circ.", 0);
				RoundnessXZ=getResult("Round", 0);
				AreaXZ=getResult("Area", 0);
			}
			
			else{
				AspectRatioXZ=0;
				circularityXZ=0;
				RoundnessXZ=0;
				AreaXZ=0;
			}
			run("Clear Results");
			ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
			}
			close("AVG_Reslice of Tcellcoordinates");
			
			
			///Cell center synapse distance
			wait(10);
			cxCell = precxCell/pw; /// pixel dimension
			cyCell = precyCell/ph; /// pixel dimension
			
			makeLine(cxCell,cyCell,xsynap,ysynap); ///all must be in pixel dimension
			run("Measure");
			angle=getResult("Angle");
			run("Clear Results");
			///transpose the image to the center

			///image middle coordinates all of them are in pixel dimensions
			Xim=imageW/2;
			Yim=imageH/2;
			XT=(Xim-cxCell);
			YT=(Yim-cyCell);
			///transpose
			run("Select None");
			run("Translate...", "x=XT y=YT interpolation=None");
			///transpose synaptic site all of them in pixel dimension
			sxt=(xsynap+XT);
			syt=(ysynap+YT);
			///rotate in XY
			run("Rotate... ", "angle="+angle+" grid=1 interpolation=None");
			selectWindow("Tcellanalysis");
			run("Select None");
			///rotate in XZ
			run("Reslice [/]...", "output=0.280 start=Top avoid");
			cxResCell=prexResCell/pw;///in pixel dimensions
			czCell=preczCell/0.28;///in pixel dimensions

			selectWindow("Reslice of Tcellanalysis");
			makeLine(cxResCell,czCell,xsynap,zsynap);///all of them are in pixels
			run("Measure");
			angleZ=getResult("Angle");
			///print(angleZ);
			
			run("Select None");
			run("Rotate... ", "angle="+angleZ+" grid=1 interpolation=None enlarge");

			run("Reslice [/]...", "output=0.110 start=Top avoid");

			close("Reslice of Tcellanalysis");
			close("Tcellanalysis");
			close("Tcell");
			wait(20);
			selectWindow("Reslice of Reslice of Tcellanalysis");
			rename("Tcellanalysis");
		
			///synaptic coordinates rotated to the rigth side of the image in pixel dimensions
			rotxsynap=(((sxt-Xim)*cos(angle))-((syt-Yim)*sin(angle))+Xim);
			rotysynap=(((syt-Xim)*sin(angle))+((syt-Yim)*cos(angle))+Xim);
			
			/// distance between the cell center and the synapse
			
			synap_cell_d=sqrt(((xsynap-cxCell)*(xsynap-cxCell))+((ysynap-cyCell)*(ysynap-cyCell)));///without rotation
			////????
			x=(cxCell/pw-10);
			y=(cyCell/ph-10);
			width=synap_cell_d*2;
			height=20/ph;
			//////////
			run("Clear Results");
			selectWindow("Tcellanalysis");
			run("Duplicate...", "title=Autocrop duplicate channels=Tcc");
			run("Z Project...", "projection=[Average Intensity]");
			run("8-bit");
			setAutoThreshold("Minimum dark");
			run("Convert to Mask");
			run("Analyze Particles...", "size=10-Infinity add");
			close("AVG_Autocrop");
			close("Autocrop");
			close("Threshold");
			selectWindow("Tcellanalysis");
			roiManager("Select", 0);
			run("Crop");
			///get the new dimensions of the cropped image
			getDimensions(Cellwidth, Cellheight, channels, slices, frames);
	
			imagewd=Cellwidth;
			imageht=Cellheight;
			xCell=Cellwidth/2;
			yCell=Cellheight/2;
			nsynapX=xCell+synap_cell_d;
			nsynapY=yCell;
			ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}
/////////////In summary: T cell should be cropped in all Z stacks and IS coordinates should be (X=image wide value, Y=Half of image height value, and Z= zsynap = centslice*0.28)
			


			close("SynapCoordinates");	
					}

/////////////////////////////////////////////////////////////
////Manual selection of tumor in 2D

			
			else {
			///tumor cell coordinate
			setTool("point");
			waitForUser("indicate T cell tumor cell interaction point (it has to be at the middle of the interaction)");
			Stack.getPosition(channel, slice, frame);
			centslice=slice;
			getPixelSize(unit,pw,ph);
			run("Measure");
			precxsynap=getResult("X", 0);
			precysynap=getResult("Y", 0);
			wait(10);
			xsynap = (precxsynap/pw);
			ysynap = (precysynap/ph);	
			run("Clear Results");
			///T cell coordinates
			setTool("polygon");
			waitForUser("Delimite the T cell forming the IS");
			roiManager("Add");
			roiManager("Select", 0);
			run("Measure");
			getPixelSize(unit,pw,ph);			
			getSelectionBounds(xCell,yCell,wCell,hCell);
			precxCell=getResult("X", 0); /// x coordinate of Cell mass center
			precyCell=getResult("Y", 0);/// y coordinate of cell mass center		
			///get cell shape parameters
			if(e1||e2||e3||e4||e5||e6 ==true){
			AspectRatio=getResult("AR", 0);
			circularity=getResult("Circ.", 0);
			Roundness=getResult("Round", 0);
			Area=getResult("Area", 0);
			}
			else{
				AspectRatio=0;
				circularity=0;
				Roundness=0;
				Area=0;
			}
			wait(10);
			cxCell = precxCell/pw;
			cyCell = precyCell/ph;
			run("Clear Outside");
			run("Clear Results");
			makeLine(cxCell,cyCell,xsynap,ysynap);
			run("Measure");
			angle=getResult("Angle");
			run("Clear Results");
			///transpose the image to the center

			///image middle coordinates
			Xim=imageW/2;
			Yim=imageH/2;
			XT=(Xim-cxCell);
			YT=(Yim-cyCell);
			///transpose
			run("Translate...", "x=XT y=YT interpolation=None");
			///waitForUser("check");

			///transpose synaptic site
			sxt=(xsynap-XT);
			syt=(ysynap-YT);
			///rotate
			run("Rotate... ", "angle="+angle+" grid=1 interpolation=None");
			///waitForUser("check");
			///synaptic coordinates rotated to the rigth side of the image
			rotxsynap=(((sxt-Xim)*cos(angle))-((syt-Yim)*sin(angle))+Xim);
			rotysynap=(((syt-Xim)*sin(angle))+((syt-Yim)*cos(angle))+Xim);
			
			/// distance between the cell center and the synapse
			synap_cell_d=sqrt(((rotxsynap-Xim)*(rotxsynap-Xim))+((rotysynap-Yim)*(rotysynap-Yim)));
			x=(cxCell/pw-10);
			y=(cyCell/ph-10);
			width=synap_cell_d*2;
			height=20/ph;
			setTool("rectangle");
			makeRectangle(x, y, width/2, height/2);
			waitForUser("delimite the cell");
			roiManager("Add");
			roiManager("Select", 1);	
			run("Crop");
			///get the new doimensions of the cropped image
			getDimensions(Cellwidth, Cellheight, channels, slices, frames);
			///waitForUser("check");		
			imagewd=Cellwidth;
			imageht=Cellheight;
			xCell=Cellwidth/2;
			yCell=Cellheight/2;
			ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}		
			}
		}
		
				
////////analysis functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////Mass Center Polarity
	
		if(a1==true){
			MCP1 = MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Green);
		}
		else{
				MCP1=0;
		}
		if(a2==true){
			MCP2 = MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Red);
		}
		else{
				MCP2=0;
		}
		if(a3==true){
			MCP3 = MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Blue);
		}
		else{
				MCP3=0;
		}
		if(a4==true){
			MCP4 = MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Cyan);
		}
		else{
				MCP4=0;
		}
		if(a5==true){
			MCP5 = MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Magenta);
		}
		else{
				MCP5=0;
		}
		if(a6==true){
			MCP6 = MassPolarity(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,synap_cell_d,Yellow);
		}
		else{
				MCP6=0;
		}
/////////Heatmap
		if(b1==true){
			HM1 = Heatmap(centslice,Green,dir,name1);
		}
		else{
				MCP1=0;
		}
		if(b2==true){
			HM2 = Heatmap(centslice,Red,dir,name1);
		}
		else{
				HM2=0;
		}
		if(b3==true){
			HM3 = Heatmap(centslice,Blue,dir,name1);
		}
		else{
				HM3=0;
		}
		if(b4==true){
			HM4 = Heatmap(centslice,Cyan,dir,name1);
		}
		else{
				HM4=0;
		}
		if(b5==true){
			HM5 = Heatmap(centslice,Magenta,dir,name1);
		}
		else{
				HM5=0;
		}
		if(b6==true){
			HM6 = Heatmap(centslice,Yellow,dir,name1);
		}
		else{
				HM6=0;
		}
/////////Syanpotic area fluorescence
		if(c1==true){
			PtjFSFT1 = SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Green);
		}
		else{
				PtjFSFT1=0;
		}
		if(c2==true){
			PtjFSFT2 = SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Red);
		}
		else{
				PtjFSFT2=0;
		}
		if(c3==true){
			PtjFSFT3 = SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Blue);
		}
		else{
				PtjFSFT3=0;
		}
		if(c4==true){
			PtjFSFT4 = SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Cyan);
		}
		else{
				PtjFSFT4=0;
		}
		if(c5==true){
			PtjFSFT5 = SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Magenta);
		}
		else{
				PtjFSFT5=0;
		}
		if(c6==true){
			PtjFSFT6 = SAF(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Yellow);
		}
		else{
				PtjFSFT6=0;
		}

////////Maxfluorescence/MFI ratio
		if(d1==true){
			MaxMFI1 = MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Green);
		}
		else{
				MaxMFI1=0;
		}
		if(d2==true){
			MaxMFI2 = MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Red);
		}
		else{
				MaxMFI2=0;
		}
		if(d3==true){
			MaxMFI3 = MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Blue);
		}
		else{
				MaxMFI3=0;
		}
		if(d4==true){
			MaxMFI4 = MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Cyan);
		}
		else{
				MaxMFI4=0;
		}
		if(d5==true){
			MaxMFI5 = MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Magenta);
		}
		else{
				MaxMFI5=0;
		}
		if(d6==true){
			MaxMFI6 = MaxMFI(centslice,rotxsynap,rotysynap,xCell,yCell,imagewd,imageht,Yellow);
		}
		else{
				MaxMFI6=0;
		}

		//////Immune synapse analysis
		if(f1==true){
			Spread1 = SpreadAnalisis(Green);
		}
		else{
			Spread1 =0;
		}
		if(f2==true){
			Spread2 = SpreadAnalisis(Red);
		}
		else{
			Spread2 =0;
		}
		if(f3==true){
			Spread3 = SpreadAnalisis(Blue);
		}
		else{
			Spread3 =0;
		}
		if(f4==true){
			Spread4 = SpreadAnalisis(Cyan);
		}
		else{
			Spread4 =0;
		}
		if(f5==true){
			Spread5 = SpreadAnalisis(Magenta);
		}
		else{
			Spread5 =0;
		}
		if(f6==true){
			Spread6 = SpreadAnalisis(Yellow);
		}
		else{
			Spread6 =0;
		}


		if(a1||a2||a3||a4||a5||a6||c1||c2||c3||c4||c5||c6 == true){
			print("Type_of_Analysis,Número,Nombre,Aspect_Ratio_XY,Aspect_Ratio_XZ,Circularity_XY,Circularity_XZ,Roundness_XY,Roundness_XZ,Area_XY,Area_XZ,Mass_Center_Polarity_Green,Mass_Center_Polarity_Red,Mass_Center_Polarity_Blue,Mass_Center_Polarity_Cyan,Mass_Center_Polarity_Magenta,Mass_Center_Polarity_Yellow,Synaptic_Fluorescence_Green,Synaptic_Fluorescence_Red,Synaptic_Fluorescence_Blue,Synaptic_Fluorescence_Cyan,Synaptic_Fluorescence_Magenta,Synaptic_Fluorescence_Yellow,MaxFluo_MFI_Ratio_Green,MaxFluo_MFI_Ratio_Red,MaxFluo_MFI_Ratio_Blue,MaxFluo_MFI_Ratio_Cyan,MaxFluo_MFI_Ratio_Magenta,MaxFluo_MFI_Ratio_Yellow");
			print("Co_culture,"+image+","+name1+","+AspectRatioXY+","+AspectRatioXZ+","+circularityXY+","+circularityXZ+","+RoundnessXY+","+RoundnessXZ+","+AreaXY+","+AreaXZ+","+MCP1+","+MCP2+","+MCP3+","+MCP4+","+MCP5+","+MCP6+","+PtjFSFT1+","+PtjFSFT2+","+PtjFSFT3+","+PtjFSFT4+","+PtjFSFT5+","+PtjFSFT6+","+MaxMFI1+","+MaxMFI2+","+MaxMFI3+","+MaxMFI4+","+MaxMFI5+","+MaxMFI6);
		}

		
	close("Tcellanalysis");
	
	}

ROI=isOpen("ROI Manager");
				if(ROI == true){
				selectWindow("ROI Manager");
				run("Close");
				}	
close(name);	


selectWindow("Log");
saveAs("text", dir+File.separator+"Results"+File.separator+"Results_Log"+year+"-"+month+"-"+dayOfMonth+".csv");
}
selectWindow("Log");
saveAs("text", dir+File.separator+"Results"+File.separator+"Results_Log"+year+"-"+month+"-"+dayOfMonth+".csv");
close("Log");

	
		

