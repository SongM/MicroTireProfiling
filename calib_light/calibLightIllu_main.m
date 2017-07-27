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

param.file_name.caliblightillu = 'whiteboard';
saveMatFile();
% convertImFromDNG_individual_folder(mat_file,param.file_name.caliblightillu);
% load(mat_file);

%% 2.D.	PS Light illumination calib
task_main = '2.D.CalibLightIllu';
displayTask(task_main,2);

% param.tire_region = [1706,1385,4200,5000];
% param.tire_region = [1126,2578,3930,5291];
% param.caliblightillu.whiteboard_region = [498,605,4602,6833];
param.caliblightillu.whiteboard_region = [200,1000,4900,7000];
param.caliblightillu.thickness_whiteboard = 0;
param.caliblightillu.filter_N = 20;
saveMatFile;

calibLightIllu(mat_file);
load(mat_file);

job_done = [job_done;task_main];
clear task_main;
saveMatFile;







