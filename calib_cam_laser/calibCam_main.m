%%
function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));
rmpath(genpath([function_folder,'\','New Folder']));
getStarted;

case_name = '20170724';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);
%%

param.file_name.calibcam = 'cam_laser_calib';
saveMatFile;

% convertImFromDNG_individual_folder(mat_file,param.file_name.calibcam,0.001);
% load(mat_file);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(mat_file);

% 2.	Calibration
% 2.A.	Camera intrinsic calib
task_main = '2.A.calibCam';
displayTask(task_main,2);

% param.calibcam.folder_name_calib = [folder_path,'\dng'];
param.calibcam.folder_name_calib = folder_path;
param.calibcam.checkerboard_size = [16,10];
% param.calibcam.checkerboard_grid_size = 9.26;
param.calibcam.checkerboard_grid_size = 8.68;
param.calibcam.sigma = 10;
param.calibcam.min_l_th = 100;
param.calibcam.max_k_th = 0.4;
param.calibcam.var_k_th = 0.001;
param.calibcam.var_l_th = 1000;
param.calibcam.over_s = 1.1;

param.calibcam.winSize = 20;

saveMatFile;

calibCam(mat_file);
load(mat_file);
load(individual_mat_file.camcalib,'bs_cam_calib');
displayCamCalib(CamCalib,bs_cam_calib,1);
clear bs_cam_calib;
% load(mat_file);
job_done = [job_done; task_main];
clear task_main;
saveMatFile;
