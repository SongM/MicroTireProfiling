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

param.caliblightpos.ball_range_scale = 1.5;
% param.caliblightpos.brightness_th = 0.8;
% param.caliblightpos.ratio_radius_brightness_spot = 0.7;

% param.caliblightpos.blackboard_region = [1286, 1697, 4100, 5842];


saveMatFile;
%% manually choose ball area
manuallyChooseBallArea(mat_file);

%% get b_ck and b_sik
manuallyChooseBallArea_NEXT_AutomaticallyFindC_k__and__s_ik(mat_file);
load(mat_file);
load(individual_mat_file.caliblightpos,'balls','ref_board_balls');

   [balls, LightPos, error] = getLightPosFromCkAndSik(balls, error);
    fprintf(['LightPos = ']);
    disp(LightPos);
    fprintf(['d_ik = ']);
    disp(error.caliblightpos.d_ik);

