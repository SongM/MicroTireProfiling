%%
function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));

rmpath(genpath([function_folder,'\','New Folder']));getStarted;

case_name = '20170630_20';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);

%%

c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
individual_mat_file.stitching = [data_folder, '\', case_name, '\processed_data\', 'tire_wear',timestamp,'.mat'];
clear c timestamp;
saveMatFile();


clear ps_set_after ps_set_before;
for i=1:3
    load(['C:\Users\s\Google Drive\projects\honda_20170118\laser_photo\20170630_20\processed_data\ps_after',num2str(i),'.mat']);
    ps_set_after{i} = ps_small;
    load(['C:\Users\s\Google Drive\projects\honda_20170118\laser_photo\20170630_20\processed_data\ps_before',num2str(i),'.mat']);
    ps_set_before{i} = ps_small;
end
fprintf('loading done.\n');

load('laserRec_result_before_45.mat', 'tire_slices_with_groove', 'groove_R');
tire_slices_with_groove_before_45 = tire_slices_with_groove;
groove_R_before_45 = groove_R;
load('laserRec_result_after_45.mat', 'tire_slices_with_groove', 'groove_R');
tire_slices_with_groove_after_45 = tire_slices_with_groove;
groove_R_after_45 = groove_R;

figure(); plot(tire_slices_with_groove_before_45(1,:),tire_slices_with_groove_before_45(2,:)-groove_R_before_45)
hold on; plot(tire_slices_with_groove_after_45(1,:),tire_slices_with_groove_after_45(2,:)-groove_R_after_45)

%%
tire_slices = tire_slices_with_groove_after_45;

N_T_ref_tire_surface_theta = 7000;
N_T_ref_tire_surface_Z = 1200;
t_T_ref_tire_surface_Z = linspace(min(tire_slices(1,:)),max(tire_slices(1,:)),N_T_ref_tire_surface_Z);

t_T_ref_tire_surface_theta = linspace(-180,180,N_T_ref_tire_surface_theta+1);
t_T_ref_tire_surface_theta(end) = [];


grid_depth_map_after = getGridDepthMap(tire_slices_with_groove_after_45, ps_set_after);
grid_depth_map_before = getGridDepthMap(tire_slices_with_groove_before_45, ps_set_before);

t_grid_depth_map_after = grid_depth_map_after;
for i=size(t_grid_depth_map_after,2):-1:1
    if (sum(~isnan(t_grid_depth_map_after(:,i)))==0)
        t_grid_depth_map_after(:,i)=[];
    end
end
t_grid_depth_map_before = grid_depth_map_before;
for i=size(t_grid_depth_map_before,2):-1:1
    if (sum(~isnan(t_grid_depth_map_before(:,i)))==0)
        t_grid_depth_map_before(:,i)=[];
    end
end
depthMap_before = t_grid_depth_map_before-groove_R_before_45;
depthMap_after = t_grid_depth_map_after-groove_R_after_45;
figure();subplot(2,1,1);imagesc_jet(depthMap_before);caxis([0,10]);
subplot(2,1,2);imagesc_jet(depthMap_after);caxis([0,10]);

figure();imagesc_jet([depthMap_before,depthMap_after]);


t1 = [depthMap_before, depthMap_before];
t3 = []; t4 = [];
for i=1:(size(t1,2)-size(depthMap_after,2)-1);
    t2 = t1(:, i:(i+size(depthMap_after,2)-1));
    depth_diff = t2 - depthMap_after;
    t3 = [t3,mean(depth_diff(~isnan(depth_diff)))];
    t4 = [t4,sum(depth_diff(:)<0)];
end

[~,ind] = min(t3);

i=806;
t2 = t1(:, i:(i+size(depthMap_after,2)-1));
depth_diff = t2 - depthMap_after;

figure();imagesc_jet(depth_diff);    
