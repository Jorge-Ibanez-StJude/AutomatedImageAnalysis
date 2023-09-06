/////colocalized area macros 
c1=1;
print("Name, Spots detected, Spot number, Spot Area");
///open folder
setOption("JFileChooser", true);
run("Set Measurements...", "area mean standard min centroid center integrated median display redirect=None decimal=3");
dir= getDirectory("Choose cell's folder");
setOption("JFileChooser", false);
imagenames=getFileList(dir); /// directorio de las células a analizar
nbimages=lengthOf(imagenames); /// nombre de las imagenes
for(image=0; image<nbimages; image++) { /// Loop de iteración de las imagenes a anlizar
	
	name=imagenames[image];
	totnamelength=lengthOf(name); /// extención del nombre
	namelength=totnamelength-4;
	name1=substring(name, 0, namelength);
	extension=substring(name, namelength, totnamelength);
	
	if(extension==".tif") {
		
		open(dir+name);
		///Canal 1 y 3 (act y prot)
		run("Duplicate...", "title=cA.tif duplicate channels="+c1);
		selectWindow("cA.tif");
		run("8-bit");
		setThreshold(150, 255);//////to modify depending on the channel
		setOption("BlackBackground", false);
		run("Convert to Mask", "method=Default background=Dark");
		///run("Convert to Mask");
		selectWindow("cA.tif");
		run("Analyze Particles...", "size=2-100 add stack");
		numROI=roiManager("count");
		for (i=0; i<numROI; i++){
			wait(10);
			roiManager("Select", i);
			run("Measure");
			Area=getResult("Area", 0);
			print(name+","+numROI+","+i+","+Area);
			run("Clear Results");
		}	
		run("Select None");	
		selectWindow("ROI Manager");
		run("Close");
		close("cA.tif");
	}
}
