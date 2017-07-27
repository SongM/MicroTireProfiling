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


task_main = '2.B.calibLaser';
displayTask(task_main,2);
param.caliblaser.redness_th = 0.6;
param.caliblaser.line_err_mse_th = 5;

saveMatFile;

calibLaser(mat_file);

load(mat_file);
job_done = [job_done; task_main];
clear task_main;
saveMatFile;
