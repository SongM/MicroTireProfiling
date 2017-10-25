# MicroTireProfiling
20170727
20171020

I. Download the program folder

II. Download the data folder

1. Pre-processing image file

1.A Install Adobe DNG Converter

Download the installment file from https://helpx.adobe.com/photoshop/using/adobe-dng-converter.html
Install Adobe DNG Converter follow the instruction.

1.B Convert Raw Image (.NEF) to .DNG format

Open Adobe DNG Converter

On Step 1, select the folder where .NEF files are.

On Step 2, select the location to save converted images.

On Step 3, select the name for converted images (suggest:img_(4 Digit Serial Number), and select File extension as .dng.

On Step 4, click 'Change Preferences...', 

  Compatibility: Custom...

  Backward Version: DNG 1.4
  
    Linear (demosaiced): uncheck
    
    Uncompressed: uncheck
    
  No need to change the remaining
  
Click Convert, the program will start to convert all .NEF file to .DNG file (it takes several seconds per image)

2. Put images into corresponding subfolders

2.A put all checkerboard images into 'cam_laser_calib' subfolder in data folder root.

2.B put all the reference checkerboard image and images taken by blackboard into 'ps_ball' subfolder in data folder root.

2.C put all whiteboard images into 'whiteboard' subfolder in data folder root.

2.D put all laser and ps images into 'ps_scan' folder.

III. open the following .m files in program folder, 

1.	‘./calib_cam_laser/calibCam_main.m’
2.	‘./calib_cam_laser/calibLaser_main.m’
3.	‘./calib_light/calibLightPos_main.m’
4.	‘./calib_light/calibLightIllu_main.m’
5.	‘./calib_light/manuallyCorrectLightPos_main.m’
6.	‘./laser_rec/laserReconstruction_main.m’
7.	‘./laser_rec/laserReconstruction_45_main.m’
8.	‘./ps/ps_tire_main_20170925.m’
9.	‘./stitching/calculate_tire_wear_main.m’

And make the following changes for each file:

1.	Line 2: change the function_folder parameter to the program folder’s address.
2.	Line 7: change the data_folder parameter to the data folder’s address.
3.	Line 6: change the case_name parameter to the name of a specific experiment.

IV. Then the programs can be run as the following instruction.

1. Calibration

Run calibration programs will create calibration parameters. Normally I will put the processed data in the data folder. But in case you want to get it by yourself, you can follow the following steps.

1.A. Camera intrinsic calib

Run ‘./calib_cam_laser/calibCam_main.m’

1.B. Laser plane calib

Run ‘./calib_cam_laser/calibLaser_main.m’

1.C. PS light position calib

Run ‘./calib_light/calibLightPos_main.m’.

The light position calibration may need some manual correction steps. If so, run ‘./calib_light/manuallyCorrectLightPos_main.m’, and follow the instruction from the program.

1.D. PS Light illumination calib

Run ‘./calib_light/calibLightIllu_main.m’

2. Experiment data processing

2.A. Laser 3D reconstruct

Run ‘./laser_rec/laserReconstruction_main.m’ for full 360 degree reconstruction
Or Run ‘./laser_rec/laserReconstruction_45_main.m’ for 45 degree reconstruction

2.B. PS 3D reconstruct

Run  ‘./ps/ps_tire_main_20170925.m’

2.C. get tire wear result

Run ‘./stitching/calculate_tire_wear_main.m’

After running this, the final result will be displayed.

