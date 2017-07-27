%%
function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));
rmpath(genpath([function_folder,'\','New Folder']));
getStarted;

case_name = '20170630_20';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);
%%

param.file_name.caliblightpos = 'ps_ball';
param.N_lights = 16;
saveMatFile;

% convertImFromDNG_individual_folder(mat_file, param.file_name.caliblightpos, 1e-3);
% load(mat_file);

%% 2.C.	PS light position calib
task_main = ['2.C.calibLightPos_',param.file_name.caliblightpos];
displayTask(task_main,2);

param.caliblightpos.d_ball_center_ref_plane = 12.68/2;
param.caliblightpos.ball_radius = 12.68/2;

param.caliblightpos.ball_size_range = [50 150];
param.caliblightpos.resize_scale = 0.1;
param.caliblightpos.ball_range_scale = 1.5;
param.caliblightpos.brightness_th = 0.8;
param.caliblightpos.ratio_radius_brightness_spot = 0.7;

param.caliblightpos.blackboard_region = [1286, 1697, 4100, 5842];


saveMatFile;

calibLightPos(mat_file);
load(mat_file);

load(individual_mat_file.caliblightpos,'balls','ref_board_balls');
displayLightPosAndBalls(LightPos, balls, ref_board_balls, CamCalib , figure_count );
clear balls ref_board_balls

if 0
    manualCorrectLightPos_main;
end

figure
dist = error.caliblightpos.d_ik;
errorbar(1:size(dist,2), mean(dist,1), mean(dist,1)-min(dist), max(dist)-mean(dist,1), 'bo', 'LineWidth', 3);
hold on
plot(median(dist),'ro', 'LineWidth', 3)
plot(mean(dist),'ko', 'LineWidth', 3)
xlabel('LED index')
ylabel('Distance [mm]')
grid on
clear dist;

% 
% temp.(light_position_file_name).LightPos = LightPos;
% temp.(light_position_file_name).error = error;

job_done = [job_done; task_main];
clear task_main;
saveMatFile;

