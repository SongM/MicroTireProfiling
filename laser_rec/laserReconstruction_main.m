%%
function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));

rmpath(genpath([function_folder,'\','New Folder']));getStarted;

case_name = '20170724';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);
%% Seperate Laser Pics and PS Pics
if (~isfield(param.file_name,'laserandps'))
    param.file_name.laserandps = 'ps_scan';
    [laser_file_list, ps_set_file_list] = seperateLaserAndPS(folder_path, param.file_name.laserandps);
    
    saveMatFile();
end
% convertImFromDNG_individual_folder(mat_file,param.file_name.laser);
% load(mat_file);
%% 3.A. Laser Reconstruction
task_main = '3.A.LaserReconstruction';
displayTask(task_main,2);
param.tire_region = [400,1800,4800,6700];


param.laserrec.laser_region = [1728,404,3119,6865];
temp_im = getImFromPath([folder_path,'/',param.file_name.laserandps,'/',laser_file_list(1).name]);
figure();imshow(temp_im/quantile(temp_im(:),1-1e-2));
title(['quantile 1e-2 = ', num2str(quantile(temp_im(:),1-1e-2))]);
rectangle('Position', [param.tire_region(1:2), param.tire_region(3:4)-param.tire_region(1:2)], 'EdgeColor', 'r');
rectangle('Position', [param.laserrec.laser_region(1:2), param.laserrec.laser_region(3:4)-param.laserrec.laser_region(1:2)], 'EdgeColor', 'g');
drawnow();
clear temp_im

%%
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
individual_mat_file.laserrec = [data_folder, '\', case_name, '\processed_data\', 'laserRec_result',timestamp,'.mat'];
clear c timestamp;


%% get strongest Features
param.laserrec.sigma = 3;
param.laserrec.dist_th = 50;

param.laserrec.N_strongest = 1000;

saveMatFile();

getStrongestFeatures(mat_file);
getMatchedPairs(mat_file);
Laser_C_ref = laser_results(1).Laser_X_c(:,1);

load(mat_file);
% checkMatchedPairs(mat_file, laser_i);

%% first time estimation of d_theta to get V
param.laserrec.C0_est = [0;0;1065];
t_V_est=[0;1;0];
param.laserrec.V_est = t_V_est/sqrt(sum(t_V_est.^2));
clear t_V_est
C0 = param.laserrec.C0_est;
V = param.laserrec.V_est;

param.laserrec.N_loop = 3;
saveMatFile();


    %% get C0 and V and d_theta
load(individual_mat_file.laserrec,'laser_results');
[laser_results, C0, V] = CalculateRotationParams(...
    laser_results, CamCalib, param.laserrec.C0_est, param.laserrec.V_est, param.laserrec.N_loop, figure_count);
fprintf('C0 and V = \n');
disp([param.laserrec.C0_est,C0,[0;0;0],param.laserrec.V_est,V]);

manualRotationCorrection = ...
    [1, 271, -2*pi;];

[laser_results, manualRotationCorrection] =...
    ManuallyCorrectRotationAngle(mat_file, laser_results, manualRotationCorrection,...
    1);

param.laserrec.manualRotationCorrection = manualRotationCorrection;
clear manualRotationCorrection

% clear manualRotationCorrection

% param.laserrec.tire_surface_region_Z = [-100,60];
% param.laserrec.N_groove = 4;
% t_laser_results = findLaserGroove...
% (t_laser_results, C0, V, CamCalib, param.laserrec.tire_surface_region_Z, param.laserrec.N_groove, Laser_C_ref);


% % V = updateVWithLaserGroove(t_laser_results, param.stitching.N_groove, C0, V, Laser_C_ref);
% 
% %% second time estimation of d_theta
% laser_resulst = t_laser_results;
% for i=1:10
% [laser_results, C0_, V_] = CalculateRotationParams(...
%     t_laser_results, CamCalib, C0, V, 1, figure_count);
% fprintf('C0 and V = \n');
% disp([C0,C0_,[0;0;0],V,V_]);
% 
% [laser_results, ~] =...
%     ManuallyCorrectRotationAngle(mat_file, laser_results, manualRotationCorrection,...
%     1);
% clear manualRotationCorrection
% laser_results = findLaserGroove...
% (laser_results, C0, V, CamCalib, param.laserrec.tire_surface_region_Z, param.laserrec.N_groove, Laser_C_ref);
% V = updateVWithLaserGroove(laser_results, param.stitching.N_groove, C0, V, Laser_C_ref);
% end
% end



% 
% manualRotationCorrection = ...
%     [1, 120, -2*pi;...
%     1, 184, -4*pi;...
%     30, 133, -2*pi;...
%     29, 199, -4*pi;...
%     38, 144, -2*pi;...
%     84, 161, -2*pi];


% 
% manualRotationCorrection = ...
%     [27, 90, -2*pi;...
%     90, 152, -2*pi;...
%     31, 94, -2*pi;...
%     94, 156, -2*pi;...
%     42, 102, -2*pi;...
%     102, 165, -2*pi;...
%     49, 110, -2*pi;...
%     110, 175, -2*pi;...
%     71, 132, -2*pi;...
%     132, 199, -2*pi];


% 
% manualRotationCorrection = ...
%     [1, 120, -2*pi;...
%     120, 185, -2*pi];

% manualRotationCorrection = ...
%     [27, 90, -2*pi;...
%     90, 152, -2*pi;];
%     31, 94, -2*pi;...
%     94, 156, -2*pi;...
%     42, 102, -2*pi;...
%     102, 165, -2*pi;...
%     49, 110, -2*pi;...
%     110, 175, -2*pi;...
%     71, 132, -2*pi;...
%     132, 199, -2*pi];




save(individual_mat_file.laserrec,'-v7.3');
clear laser_results_forC0andV laser_results_beforeCorrection laser_results
clear C0_ V_ d_thetas_beforeCorrection thetas_beforeCorrection
saveMatFile();

%% laserRec
load(individual_mat_file.laserrec,'laser_results');

Laser_C_ref = laser_results(1).Laser_X_c(:,1);
d_theta_all = [laser_results.d_theta];
b_display = 0;

figure_count = figure_count + 4;
for laser_i=1:numel(laser_results)
    Laser_x_C = laser_results(laser_i).Laser_X_c;
    Laser_Cylinder = convertEuclideanToCylinder(Laser_x_C, C0, V, Laser_C_ref);
    laser_results(laser_i).theta = sum(d_theta_all(1:laser_i))/pi*180;
    
    Laser_Cylinder(2,:) = Laser_Cylinder(2,:) + laser_results(laser_i).theta;
    laser_results(laser_i).Laser_Cylinder = Laser_Cylinder;
    Laser_E = convertCylinderToEuclidean(Laser_Cylinder,C0,V,Laser_C_ref);
    laser_results(laser_i).Laser_E = Laser_E;
    
    if (b_display)

        %     Laser_Cylinder_all = [Laser_Cylinder_all, Laser_Cylinder];
        %     figure();scatterPoints(Laser_Cylinder);
            figure(figure_count-3);hold on;
            scatter3(Laser_Cylinder(3,:),Laser_Cylinder(1,:),Laser_Cylinder(2,:),3,'r','filled');pbaspect([2,1,1]);
            title('Cylinder Coordinate System');
        %     figure(figure_count-2);hold on;
        %     scatter3(Laser_Cylinder(3,:),Laser_Cylinder(1,:),Laser_Cylinder(2,:)*pi/180*340,3,'r','filled');
            figure(figure_count-1);hold on;
            scatter3(Laser_E(2,:),Laser_E(3,:),Laser_E(1,:),3,'r','filled');
            axis equal;
            view(0,-70);
        %     xlim([-120,120]);
            title('Camera Coordinate System');
    end
end

clear d_theta_all laser_i Laser_x_C Laser_Cylinder Laser_E


save(individual_mat_file.laserrec,'-v7.3');
clear laser_results
saveMatFile();





%%
manuallyCheckC0andV(mat_file, current_laser_i, C0_guess, V_guess)



