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
if (~isfield(param.file_name,'laserandps'))
    param.file_name.laserandps = 'ps_scan_2s_iso100_and_laser_1';
    [laser_file_list, ps_set_file_list] = seperateLaserAndPS(folder_path, param.file_name.laserandps);
    saveMatFile();
end

%% 3 photometric stereo
task_main = '3.A.PhotometricStereo_GetSurfaceNormal';
displayTask(task_main,2);
param.ps.tire_region = param.tire_region;
param.ps.tire_region(1:2) = max([param.caliblightillu.whiteboard_region(1:2);param.tire_region(1:2)]);
param.ps.tire_region(3:4) = min([param.caliblightillu.whiteboard_region(3:4);param.tire_region(3:4)]);
% param.ps.tire_region = [1121,2343,3950,5015];
param.ps.s_RGB = normalizeColVector([1;1;1]);
saveMatFile();

DROP_RATIO = 5;
load(individual_mat_file.laserrec,'laser_results');
tic;
for ps_i = 2:20
    close all;

    ps_file_list = ps_set_file_list(ps_i);
    ps_res = [];
    [ps_res, figure_count] = getScaledImFromPSFileList(ps_file_list, folder_path, param, CamCalib, LightPos, lightStrength_mat_file, V, C0, DROP_RATIO, figure_count);
    [ps_res, figure_count] = getSurfaceNormalFromScaledIm(ps_res, ...
        param, CamCalib, LightPos, figure_count);

    % save individual ps file.
    [ps_res.Z,ps_res.theta_shift] = getZFromSurfaceNormal(ps_res, C0, V, ps_res.region_due_to_size_drop, CamCalib, laser_results, Laser_C_ref, 10, 1, 1, 0);
    % ps_set_res(ps_i) = ps_res;
    save(['ps_res_',num2str(ps_i)],'ps_res','-v7.3');
    toc;
end
% getSurfaceNormal_PS(mat_file);



