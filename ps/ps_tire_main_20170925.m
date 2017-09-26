function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));
rmpath(genpath([function_folder,'\','New Folder']));
getStarted;

case_name = '20170630_20';
% case_name = '20170912';
% case_name = '20170809';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);
load(mat_file);
load(individual_mat_file.laserrec,'laserrec','tire_slices','tire_slices_with_groove','T_groove_range_Z');
% load('C:\Users\s\Google Drive\projects\honda_20170118\laser_photo\20170809\processed_data\laserRec_result_2017_8_21_17_5.mat', 'laserrec');
% load('C:\Users\s\Google Drive\projects\honda_20170118\laser_photo\20170809\processed_data\20170809_2017_8_28_7_27.mat', 'ref_tire_surface')
%%
% param.caliblightillu.whiteboard_region = [200,200,4900,7300];
param.caliblightillu.whiteboard_region = [1300,1700,4100,5800];
% param.caliblightillu.whiteboard_region = [1422,1709,3354,6377];

param.caliblightillu.thickness_whiteboard = 3;
% whiteboard_region = param.caliblightillu.whiteboard_region;

param.tire_region = param.caliblightillu.whiteboard_region;
param.ps.tire_region(1:2) = max([param.caliblightillu.whiteboard_region(1:2);param.tire_region(1:2)]);
param.ps.tire_region(3:4) = min([param.caliblightillu.whiteboard_region(3:4);param.tire_region(3:4)]);
param.ps.s_RGB = normalizeColVector([1;1;1]);
tire_region = param.ps.tire_region;
DROP_RATIO = 8;
sRGB = [1;1;1];
MIN_NUM_IM = 5; 

N_lights = param.N_lights;
%% whiteboard
param.file_name.caliblightillu = 'whiteboard';

whiteboard_small_info = getWhiteboardSmallInformation(folder_path, param, CamCalib, LightPos, DROP_RATIO);
whiteboard_small = whiteboard_small_info.whiteboard_small;
cosThetai_whiteboard = whiteboard_small_info.cosThetai_whiteboard;
xyzC_whiteboard = whiteboard_small_info.xyzC_whiteboard;
LightUnitVector = whiteboard_small_info.LightUnitVector;

NUM_X = size(xyzC_whiteboard,2);
NUM_Y = size(xyzC_whiteboard,1);

mx = NUM_X/( max(max(xyzC_whiteboard(:,:,1))) - min(min(xyzC_whiteboard(:,:,1))) );
my = NUM_Y/( max(max(xyzC_whiteboard(:,:,2))) - min(min(xyzC_whiteboard(:,:,2))) );
dx = double((1/mx + 1/my)/2);
clear mx my;
ViewingUnitVector = -xyzC_whiteboard./repmat(sqrt(sum(xyzC_whiteboard.^2,3)),1,1,3);

%% laser
% ind = tire_slices_with_groove(1,:)>T_groove_range_Z(3,1) &...
%     tire_slices_with_groove(1,:)<T_groove_range_Z(3,2);
% ref_groove_R = min(tire_slices_with_groove(2,ind));
ref_groove_R = groove_R;
% a = [laserrec.M_Laser];

% N_T_ref_tire_surface_theta = ceil(7000/360*(max(a(2,:))-min(a(2,:))));
N_T_ref_tire_surface_theta = 7000;
N_T_ref_tire_surface_Z = 1200;
t_T_ref_tire_surface_Z = linspace(min(tire_slices(1,:)),max(tire_slices(1,:)),N_T_ref_tire_surface_Z);
t_T_ref_tire_surface_R = smooth(ref_tire_surface(t_T_ref_tire_surface_Z));
% t_T_ref_tire_surface_theta = linspace(min(a(2,:)),max(a(2,:)),N_T_ref_tire_surface_theta+1);
t_T_ref_tire_surface_theta = linspace(-180,180,N_T_ref_tire_surface_theta+1);
t_T_ref_tire_surface_theta(end) = [];

[t_T_ref_tire_surface_Z_3d,t_T_ref_tire_surface_theta_3d] = meshgrid(t_T_ref_tire_surface_Z,t_T_ref_tire_surface_theta);
t_T_ref_tire_surface_R_3d =  ones(N_T_ref_tire_surface_theta,1)*t_T_ref_tire_surface_R';

% figure();pcshow([t_T_ref_tire_surface_R_3d(:),t_T_ref_tire_surface_theta_3d(:),t_T_ref_tire_surface_Z_3d(:)]);
t_T_ref_tire_surface = [t_T_ref_tire_surface_R_3d(:),t_T_ref_tire_surface_theta_3d(:),t_T_ref_tire_surface_Z_3d(:)];

t_Obj_ref_tire_surface = convertCylinderToEuclidean(t_T_ref_tire_surface, tire_axis, Laser_C_ref);
X = xyzC_whiteboard(:,:,1);
Y = xyzC_whiteboard(:,:,2);

ind = t_Obj_ref_tire_surface(1,:) > max(X(:))|...
    t_Obj_ref_tire_surface(1,:) < min(X(:))|...
    t_Obj_ref_tire_surface(2,:) > max(Y(:))|...
    t_Obj_ref_tire_surface(2,:) < min(Y(:))|...
    t_Obj_ref_tire_surface(3,:) > tire_axis.t(3);
t_Obj_ref_tire_surface(:,ind) = [];
figure();pcshow(t_Obj_ref_tire_surface');title('t_Obj_ref_tire_surface');
fprintf('laser_done.\n');
%%
% ps_small.region = whiteboard_region;
t_xn_ref_tire_surface = t_Obj_ref_tire_surface(1:2,:) ./ (ones(2,1)*t_Obj_ref_tire_surface(3,:));
t_p_ref_tire_surface = inv_normalize_pixel(t_xn_ref_tire_surface,CamCalib);

t_p_ref_tire_surface_after_drop_resolution = (t_p_ref_tire_surface-1)/DROP_RATIO + 1;
    
    temp_p = t_p_ref_tire_surface_after_drop_resolution -...
        floor( (tire_region(1:2)-1) / DROP_RATIO + 1)'*ones(1,size(t_p_ref_tire_surface_after_drop_resolution,2)) + 1;
    round_temp_p = round(temp_p);
    round_error = sum( (round_temp_p - temp_p).^2);
    
    Z_measured_ind = zeros(NUM_Y,NUM_X);
    Z_measured = nan(NUM_Y,NUM_X);
    error_map = 2*ones(NUM_Y,NUM_X);
    for laser_p_i = 1:size(t_Obj_ref_tire_surface,2)
        jj = round_temp_p(1,laser_p_i);
        ii = round_temp_p(2,laser_p_i);
        if (error_map(ii,jj)>round_error(laser_p_i))
            Z_measured_ind(ii,jj) = 1;
            Z_measured(ii,jj) = t_Obj_ref_tire_surface(3,laser_p_i);
            error_map(ii,jj) = round_error(laser_p_i);
        end
    end
    figure();imagesc(Z_measured);axis equal;title('Z measured');
    q_ref = ( Z_measured(2:end,:) - Z_measured(1:(end-1),:) ) ./ ( Y(2:end,:) - Y(1:(end-1),:) );
    q_ref(end+1,:) = nan;
    p_ref = ( Z_measured(:,2:end) - Z_measured(:,1:(end-1)) ) ./ ( X(:,2:end) - X(:,1:(end-1)) );
    p_ref(:,end+1) = nan;
    t_l_ref = sqrt(p_ref.^2 + q_ref.^2+1);
    nMap_ref(:,:,1) = p_ref./t_l_ref;
    nMap_ref(:,:,2) = q_ref./t_l_ref;
    nMap_ref(:,:,3) = -1./t_l_ref;
    p_ref( (abs(p_ref)>tand(85))|(abs(q_ref)>tand(85)) ) = nan;
    q_ref( (abs(p_ref)>tand(85))|(abs(q_ref)>tand(85)) ) = nan;
        
    figure();imshow((1-nMap_ref)/2);title('nMap ref');
    figure();
    subplot(1,2,1);imagesc(p_ref);axis equal;colorbar;colormap jet;title('p ref');
    subplot(1,2,2);imagesc(q_ref);axis equal;colorbar;colormap jet;title('q ref');
    figure();pcshow([X(:),Y(:),Z_measured(:)]);
fprintf('Z_measured done.\n');

%% ps scan
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
disp('ps_tire')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
ps_i = 1;

for ps_i = 2:numel(ps_set)
ps = ps_set(ps_i);
ps.region = param.tire_region;

ps_scan_ims = [];
fprintf('load_ims_');
for light_i = 1:N_lights
    fprintf([num2str(light_i),'..']);
    ps_scan_ims(:,:,:,light_i) = getImFromPath(ps.im_filepath_list{light_i});
    if(light_i==1);
        figure();imshow(ps_scan_ims(:,:,:,1)./quantile(ps_scan_ims(:),0.99));hold on;
        rectangle('Position', [tire_region(1:2), tire_region(3:4)-tire_region(1:2)], 'EdgeColor', 'r');
        drawnow();
    end

end
fprintf('done.\n');
ims_small = ps_scan_ims(...
param.tire_region(2):DROP_RATIO:param.tire_region(4),...
param.tire_region(1):DROP_RATIO:param.tire_region(3), :, :);

%%
ps_small = ps;
ps_small.DROP_RATIO = DROP_RATIO;
ps_small.ims_small = ims_small;

%% scaled_im
filter_size = 4;
filter = fspecial('gaussian', filter_size, filter_size/2);
scaled_im = zeros(NUM_Y,NUM_X,3,N_lights);
for color_i=1:3
    scaled_im(:,:,color_i,:) = imfilter(ims_small(:,:,color_i,:),filter) ./ ...
        imfilter(whiteboard_small(:,:,color_i,:),filter) .* cosThetai_whiteboard;
end
fprintf('scaled_ims done.\n');
%% invalidMapEst
% SHADOW_THRES = 0.03;
SHADOW_THRES = 6e-4;
% SHADOW_THRES = 1/254;
SATUATE_THRES = 0.99;
invalidMapEst = estimateInvalidMap( whiteboard_small_info.whiteboard_small, ims_small, SHADOW_THRES, SATUATE_THRES, 1 );
figure();   imagesc(sum(invalidMapEst,3)); colormap jet;   colorbar;   axis equal; title('invalidMapEst');
fprintf('invalidMapEst.\n');
%%
% quantile_th = 0.01;
color_i = 1;
scaled_im_color_i = squeeze(scaled_im(:,:,color_i,:));
dark_th = min(quantile(scaled_im_color_i(~invalidMapEst), 0.01), 0.03);
% bright_th = max(quantile(scaled_im(:), 1-0.001), 1);
bright_th = max(quantile(scaled_im_color_i(~invalidMapEst), 1-0.001));

scaled_im_darkMap = squeeze( ...
    (scaled_im(:,:,1,:)<dark_th)...
    | (scaled_im(:,:,2,:)<dark_th)...
    | (scaled_im(:,:,3,:)<dark_th) );

scaled_im_brightMap = squeeze( ...
    (scaled_im(:,:,1,:)>bright_th)...
    | (scaled_im(:,:,2,:)>bright_th)...
    | (scaled_im(:,:,3,:)>bright_th) );

invalidMap = invalidMapEst | scaled_im_darkMap | scaled_im_brightMap;
figure();   imagesc(sum(invalidMap,3)); colormap jet;   colorbar;   axis equal; title('invalidMap');
fprintf('invalidMap.\n');

%%
% method_mixed_simple = 'mixed_simple';
%%
method_diffuse = 'diffuse';
tic;[nMap_d, kdMap_d, ksMap_d, alphaMap_d, errorMap_d, rerenderIm_d] = ...
    getSurfaceNormal_Map(scaled_im, LightUnitVector, ViewingUnitVector, invalidMap, method_diffuse);toc;
[connected_region_d, p_d, q_d, ~, ~, ~, ~, ~, calculatedZ_d]...
    = displayPSResults(nMap_d, kdMap_d, ksMap_d, alphaMap_d, errorMap_d, method_diffuse, xyzC_whiteboard, dx);
figure();imagesc_jet(acosd(sum(nMap_d.*nMap_ref,3)));caxis([0,10]);title('angular error, d');
%%

method_paper_it = 'paper_it';
tic;[nMap_pi, kdMap_pi, ksMap_pi, alphaMap_pi, errorMap_pi, rerenderIm_pi] =...
    getSurfaceNormal_Map(scaled_im, LightUnitVector, ViewingUnitVector, invalidMap, method_paper_it);toc;
[connected_region_pi, p_pi, q_pi, ~, ~, ~, ~, ~, calculatedZ_pi]...
    = displayPSResults(nMap_pi, kdMap_pi, ksMap_pi, alphaMap_pi, errorMap_pi, method_paper_it, xyzC_whiteboard, dx);
figure();imagesc_jet(acosd(sum(nMap_pi.*nMap_ref,3)));caxis([0,10]);title('angular error, pi');


%%


a=histogram(connected_region_pi,[1:max(connected_region_pi(:))]);
surface_region_ind = find(a.Values>mean(a.Values));
surface_region = zeros(size(connected_region_pi));
for region_ind = surface_region_ind
    surface_region = surface_region | (connected_region_pi == region_ind);
end
    
    
stretched_factor = max(tire_slices(2,:))/180*pi;

calculatedZ_tire = getZfromPQAndLaser_region(p_pi,q_pi,surface_region,dx, Z_measured, Z_measured_ind, 0.01,1);
P_Obj = [X(~isnan(calculatedZ_tire)), Y(~isnan(calculatedZ_tire)), calculatedZ_tire(~isnan(calculatedZ_tire))];
P_T = convertEuclideanToCylinder(P_Obj, tire_axis, Laser_C_ref);
theta_shift = laserrec(find([laserrec.ind]==ps_small.laser_im_ind)).theta;

P_T(2,:) = (P_T(2,:)+theta_shift)*stretched_factor;
ind = zeros(1,size(P_T,2));
for i=1:4
    ind = ind | ( (P_T(3,:)>T_groove_range_Z(i,1)) & (P_T(3,:)<T_groove_range_Z(i,2)));

end
P_T(:,ind) = [];
P_Obj(ind,:) = [];
figure();pcshow(P_T');
figure();pcshow(P_Obj);

ps_small.nMap_pi = nMap_pi;
ps_small.kdMap_pi = kdMap_pi;
ps_small.ksMap_pi = ksMap_pi;
ps_small.alphaMap_pi = alphaMap_pi;
ps_small.errorMap_pi = errorMap_pi;
ps_small.theta = theta_shift;
ps_small.P_Obj_pi = P_Obj;
ps_small.P_T_pi = P_T;
ps_small.SHADOW_THRES = SHADOW_THRES;
ps_small.SATUATE_THRES = SATUATE_THRES;
ps_small.filter_size = filter_size;
save(['C:\Users\s\Google Drive\projects\honda_20170118\laser_photo\20170630_20\processed_data\ps_before',num2str(ps_i),'.mat'],'ps_small','-v7.3');
close all;

end