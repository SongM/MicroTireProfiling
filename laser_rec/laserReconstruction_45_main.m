%%
function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));

rmpath(genpath([function_folder,'\','New Folder']));getStarted;

case_name = '20170630_20';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);
%% Seperate Laser Pics and PS Pics



if (~isfield(param.file_name,'laserandps'))
%%
    task = '3.A.1.SeperateLaserAndPSImage';
    displayTask(task,4);

%     param.file_name.laserandps = 'laserandps';
    param.file_name.laserandps = 'ps_scan_before_45';
    [laserrec, ps_set] = seperateLaserAndPS(folder_path, param.file_name.laserandps);    
    job_done = [job_done;task];
    clear task;
    saveMatFile();
    %%
end


% convertImFromDNG_individual_folder(mat_file,param.file_name.laser);
% load(mat_file);
%% 3.A. Laser Reconstruction
task_main = '3.A.LaserReconstruction';
displayTask(task_main,2);
param.tire_region = [700,1800,4300,6000];
param.laserrec.laser_region = [1002,1839,2020,5867];
param.laserrec.redness_th = 0.05;

temp_im = getImFromPath(laserrec(1).im_filepath);
[p_laser,C_Laser] = getOneLaserScanResult(temp_im,LaserCalib,CamCalib,...
    param.laserrec.laser_region,param.laserrec.redness_th,1);
figure();imshow(temp_im/quantile(temp_im(:),1-1e-2));
title(['quantile 1e-2 = ', num2str(quantile(temp_im(:),1-1e-2))]);
rectangle('Position', [param.tire_region(1:2), param.tire_region(3:4)-param.tire_region(1:2)], 'EdgeColor', 'r');
rectangle('Position', [param.laserrec.laser_region(1:2), param.laserrec.laser_region(3:4)-param.laserrec.laser_region(1:2)], 'EdgeColor', 'g');
drawnow();
clear temp_im p_laser C_Laser

%%
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
individual_mat_file.laserrec = [data_folder, '\', case_name, '\processed_data\', 'laserRec_result',timestamp,'.mat'];
clear c timestamp;


%% get strongest Features
task = '3.B.1.i.GetStrongestFeatures';
displayTask(task,4);

param.laserrec.sigma = 3;
param.laserrec.dist_th = 50;
param.laserrec.N_strongest = 1000;
Laser_C_ref = [0;0;0];
%%
sigma = param.laserrec.sigma;
dist_th = param.laserrec.dist_th;
redness_th = param.laserrec.redness_th;
tire_region = param.tire_region;
laser_region = param.laserrec.laser_region;
N_strongest = param.laserrec.N_strongest;

filter = fspecial('Gaussian',4*sigma+1,sigma);
SAT_PERCENT=0.01;
b_display = 0;
N_laser = numel(laserrec);

saveMatFile();
tic;
% for laser_i = 1:N_laser
for laser_i = 1:N_laser
    im = getImFromPath(laserrec(laser_i).im_filepath);

    figure_count = figure_count + b_display;
    [~, laserrec(laser_i).p_laser, laserrec(laser_i).C_Laser,...
        laserrec(laser_i).features, laserrec(laser_i).validPoints,...
        laserrec(laser_i).sat_val]...
        = getStrongestFeatureFromIm(im, filter, N_strongest, dist_th,...
        redness_th, LaserCalib, CamCalib, laser_region, tire_region,...
        b_display, figure_count, SAT_PERCENT);
    figure_count = figure_count + b_display;

    laserrec(laser_i).SAT_PERCENT = SAT_PERCENT;
    time = toc;
    fprintf(['working on ' num2str(laser_i),'/',num2str(N_laser), ...
        ', sat_val = ', num2str(laserrec(laser_i).sat_val),...
        ', time remain = ',num2str((N_laser-laser_i)*time/laser_i),'.\n']);
end
clear laser_i im sigma dist_th tire_region laser_region N_strongest 
clear b_display SAT_PERCENT N_laser filter redness_th time

job_done = [job_done;task];
clear task;



%% get matched pairs
task = '3.B.1.ii.GetMatchedPairs';
displayTask(task,4);

laser_i = 1;
laserrec(laser_i).indexPairs=[];
laserrec(laser_i).matchedPoints1=[];
laserrec(laser_i).matchedPoints2=[];
fprintf(['laser_',num2str(laser_i),': N_pairs = ',...
    num2str(size(laserrec(laser_i).indexPairs,1)),...
    ', sat_val = ', num2str(laserrec(laser_i).sat_val), '.\n']);
N_laser = numel(laserrec);
%%
for laser_i = 2:N_laser
    lr1 = laserrec(laser_i-1);
    lr2 = laserrec(laser_i);
    indexPairs = getPairsInd(lr1, lr2, CamCalib);

    laserrec(laser_i).matchedPoints1 = lr1.validPoints(indexPairs(:,1),:);
    laserrec(laser_i).matchedPoints2 = lr2.validPoints(indexPairs(:,2),:);

    fprintf(['laser_',num2str(laser_i),': N_pairs = ',...
    num2str(size(laserrec(laser_i).matchedPoints1,1)),...
    ', sat_val = ', num2str(laserrec(laser_i).sat_val), '.\n']);
    clear lr1 lr2 indexPairs
end
clear laser_i N_laser
job_done = [job_done;task];
clear task;
%%
if(0)
%%
    laser_i = 108;
    redness_th = param.laserrec.redness_th;
% checkMatchedPairs(mat_file, laser_i);
    laserrec(laser_i) = checkMatchedPairs(laserrec, param, laser_i, figure_count, CamCalib, LaserCalib, redness_th);
%%
end

%% first time estimation of d_theta to get V
task = '3.B.2.i.FirstTimeEstimate_d_theta_toGetC2Obj_rot_axis';
displayTask(task,4);
% param.laserrec.C2ObjRotAxis_est.t = [0;0;1/ref_board_balls.plane_eqn(3)+6+332];
param.laserrec.C2ObjRotAxis_est.t = [0;0;1/ref_board_balls.plane_eqn(3)+6+374];
 % t_V_est=[0;1;0];
t_V_est=[-0.028274;0.99843;0.05222];
param.laserrec.C2ObjRotAxis_est.V = t_V_est/sqrt(sum(t_V_est.^2));
clear t_V_est

param.laserrec.N_loop = 1;

[laserrec, C2ObjRotAxis] = CalculateRotationParams(...
    laserrec, CamCalib, param.laserrec.C2ObjRotAxis_est, param.laserrec.N_loop, figure_count);
fprintf('C0 and V = \n');
disp([param.laserrec.C2ObjRotAxis_est.t,C2ObjRotAxis.t,[0;0;0],param.laserrec.C2ObjRotAxis_est.V,C2ObjRotAxis.V]);
job_done = [job_done;task];
clear task;

save(individual_mat_file.laserrec,'-v7.3');
clear laserrec
saveMatFile();
load(individual_mat_file.laserrec,'laserrec');


%% correct d_theta with match features of 1 full rotation.
% task = '3.B.2.ii.correct_d_theta with match features of 1 full rotation';
% displayTask(task,4);
% if(0)
% getAngelDifferenceBetweenTwoFrames(mat_file, 1, 124, 0);
% getAngelDifferenceBetweenTwoFrames(mat_file, 124, 191, 0);
% getAngelDifferenceBetweenTwoFrames(mat_file, 1, 120, 1);
% getAngelDifferenceBetweenTwoFrames(mat_file, 120, 183, 1);
% getAngelDifferenceBetweenTwoFrames(mat_file, 120, 184, 1);
% 
% end
% 
% 
% manualRotationCorrection = [...
%     1, 120, -2*pi;...
%     120, 184, -2*pi;...
%     ];
% 
% [laserrec, manualRotationCorrection] =...
%     ManuallyCorrectRotationAngle(mat_file, laserrec, manualRotationCorrection,...
%     1);
% 
% param.laserrec.manualRotationCorrection = manualRotationCorrection;
% clear manualRotationCorrection
% job_done = [job_done;task];
% clear task;

%% laserRec
task = '3.C.1.laserRec';
displayTask(task,3);

d_theta_all = [laserrec.d_theta];

for laser_i=1:numel(laserrec)
    C_Laser = laserrec(laser_i).C_Laser;
    M_Laser = convertEuclideanToCylinder(C_Laser, C2ObjRotAxis, Laser_C_ref);
    laserrec(laser_i).theta = sum(d_theta_all(1:laser_i))/pi*180;
    
    M_Laser(2,:) = M_Laser(2,:) + laserrec(laser_i).theta;
    laserrec(laser_i).M_Laser = M_Laser;
    Obj_Laser = convertCylinderToEuclidean(M_Laser,C2ObjRotAxis,Laser_C_ref);
    laserrec(laser_i).Obj_Laser = Obj_Laser;
end

clear d_theta_all laser_i C_Laser M_Laser Obj_Laser
job_done = [job_done;task];
clear task;

figure();pcshow([laserrec.M_Laser]');
view([0,0]);
save(individual_mat_file.laserrec,'-v7.3');
clear laserrec;
saveMatFile();
load(individual_mat_file.laserrec,'laserrec');




%%
if(0)
%%
    current_laser_i = 2;
    t_V = C2ObjRotAxis.V + 0*[-0.01;0;0];
    manuallyCheckC0andV(mat_file, current_laser_i, C2ObjRotAxis.t,t_V)
    clear t_V
%%
end
%% seperate groove And surface
fprintf('seperate groove and surface\n');
param.stitching.ref_groove_region_Z =  [-20,10];
param.stitching.tire_surface_region_Z = [-140,100];

param.stitching.groove_region_half_size_Z = 6;
param.stitching.N_groove = 4;
param.stitching.fitting_err_max_th = 0.5;
laserrec = seperateGrooveAndSurfacePixel(laserrec, param, C2ObjRotAxis, Laser_C_ref, CamCalib);

%% get tire_axis
disp('get tire axis');
groove_i  = 3;
[tire_axis, ~] = findTireAxis_45(groove_i, laserrec, param, C2ObjRotAxis, Laser_C_ref);
T_groove_range_Z = zeros(param.stitching.N_groove, 2);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('other groove for reference');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

for groove_i = 1:param.stitching.N_groove
    [~, T_groove_range_Z(groove_i,:)] = findTireAxis_45(groove_i, laserrec, param, tire_axis, Laser_C_ref);
end

%% get ref surface
fprintf('get ref surface\n');
tic;
[ref_tire_surface, tire_slices] = getRefTireSurf(laserrec, param, T_groove_range_Z, tire_axis, Laser_C_ref);
[ref_tire_surface_with_groove, tire_slices_with_groove] = ...
    getRefTireSurf_with_groove(laserrec, param, T_groove_range_Z, tire_axis, Laser_C_ref);
%%
fprintf('display ref surface');

a=[laserrec.M_Laser];

N_T_ref_tire_surface_theta = ceil(4000/360*(max(a(2,:))-min(a(2,:))));
N_T_ref_tire_surface_Z = 400;
t_T_ref_tire_surface_Z = linspace(...
    max( min(a(3,:)), param.stitching.tire_surface_region_Z(1) ),...
    min( max(a(3,:)), param.stitching.tire_surface_region_Z(2) ),...
    N_T_ref_tire_surface_Z);
t_T_ref_tire_surface_R = ref_tire_surface(t_T_ref_tire_surface_Z);
% t_T_ref_tire_surface_R = ref_tire_surface_with_groove(t_T_ref_tire_surface_Z);
t_T_ref_tire_surface_theta = linspace(min(a(2,:)),max(a(2,:)),N_T_ref_tire_surface_theta+1);
t_T_ref_tire_surface_theta(end) = [];

[t_T_ref_tire_surface_Z_3d,t_T_ref_tire_surface_theta_3d] = meshgrid(t_T_ref_tire_surface_Z,t_T_ref_tire_surface_theta);
t_T_ref_tire_surface_R_3d =  ones(N_T_ref_tire_surface_theta,1)*t_T_ref_tire_surface_R';

figure();pcshow([t_T_ref_tire_surface_R_3d(:),t_T_ref_tire_surface_theta_3d(:),t_T_ref_tire_surface_Z_3d(:)]);
t_Obj_ref_tire_surface = convertCylinderToEuclidean([t_T_ref_tire_surface_R_3d(:),t_T_ref_tire_surface_theta_3d(:),t_T_ref_tire_surface_Z_3d(:)], tire_axis, Laser_C_ref);
figure();pcshow(t_Obj_ref_tire_surface');

groove_R = min(ref_tire_surface_with_groove(linspace(T_groove_range_Z(3,1),T_groove_range_Z(3,2),100)));
figure();subplot(1,2,1);
pcshow([laserrec.Obj_Laser]');
hold on;
N_patch = 300;
t_groove_plane_theta = linspace(min(a(2,:)),max(a(2,:)),N_patch);
for patch_i = 1:(N_patch-1)
groove_plane_patch_T = [...
    groove_R, t_groove_plane_theta(patch_i), T_groove_range_Z(3,1);...
    groove_R, t_groove_plane_theta(patch_i+1), T_groove_range_Z(3,1);...
    groove_R, t_groove_plane_theta(patch_i+1), T_groove_range_Z(3,2);...
    groove_R, t_groove_plane_theta(patch_i), T_groove_range_Z(3,2)];
groove_plane_patch_Obj = convertCylinderToEuclidean(groove_plane_patch_T, tire_axis, Laser_C_ref)';
hold on; fill3(groove_plane_patch_Obj(:,1),groove_plane_patch_Obj(:,2),groove_plane_patch_Obj(:,3),...
    'b','FaceAlpha',0.5,'EdgeAlpha',0);
end

subplot(1,2,2);
pcshow(convertEuclideanToCylinder([laserrec.Obj_Laser], tire_axis, Laser_C_ref)');
groove_plane = [...
    groove_R, min(a(2,:)), T_groove_range_Z(3,1);...
    groove_R, max(a(2,:)), T_groove_range_Z(3,1);...
    groove_R, max(a(2,:)), T_groove_range_Z(3,2);...
    groove_R, min(a(2,:)), T_groove_range_Z(3,2)];
hold on; fill3(groove_plane(:,1),groove_plane(:,2),groove_plane(:,3),...
    'b','FaceAlpha',0.5,'EdgeAlpha',0);


clear N_T_ref_tire_surface_theta t_T_ref_tire_surface_Z t_T_ref_tire_surface_R t_T_ref_tire_surface_theta
clear t_T_ref_tire_surface_Z_3d t_T_ref_tire_surface_R_3d t_T_ref_tire_surface_theta_3d
clear t_Obj_ref_tire_surface
clear a ans groove_i groove_plane groove_plane_patch_Obj groove_plane_patch_T
clear N_patch N_T_ref_tire_surface_Z patch_i t_groove_plane_theta

save(individual_mat_file.laserrec,'-v7.3');
toc;
clear laserrec;
clear tire_slices tire_slices_with groove T_groove_range_Z
saveMatFile();
load(individual_mat_file.laserrec,'laserrec');
