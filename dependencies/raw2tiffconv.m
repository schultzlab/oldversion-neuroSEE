addpath('/Applications/Fiji.app/scripts/');
javaaddpath '/Applications/MATLAB_R2018a.app/java/mij.jar'
javaaddpath '/Applications/MATLAB_R2018a.app/java/ij.jar'
Miji('false');



open("/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.09.03/2P/20180903_15_14_16_2P/20180903_15_14_16_2P_XYT.raw");
run("Make Substack...", "  slices=1-7432-1");
saveAs("Tiff", "/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.09.03/2P/20180903_15_14_16_2P/20180903_15_14_16_2P_XYT_1.tif");

selectWindow("20180903_15_14_16_2P_XYT.raw");
run("Make Substack...", "  slices=7433-14864-1");
saveAs("Tiff", "/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.09.03/2P/20180903_15_14_16_2P/20180903_15_14_16_2P_XYT_1.tif");