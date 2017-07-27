function [Z, Laser_E_all] = getZFromSurfaceNormal(ps_res, C0, V, tire_region, CamCalib, laser_results, Laser_C_ref, h1, h2, h3, b_display)

    if (nargin==7) 
        h1 = 10;
        h2 = 1;
        h3 = 1;
        b_display = 0;
    end
    if(nargin==10)
        b_display = 0;
    end
    ps = ps_res;
    clear ps_res;
    n = ps.n;

    rerenderIm = ps.rerenderIm;
    
    p = -n(:,:,1)./n(:,:,3);
    q = -n(:,:,2)./n(:,:,3);
    p(abs(p)>tand(88)) = nan;
    q(abs(q)>tand(88)) = nan;
    p = inpaint_nans(p);
    q = inpaint_nans(q);


    [NUM_Y, NUM_X, ~] = size(n);

    mx = NUM_X/( max(ps.xC(:)) - min(ps.xC(:)) );
    my = NUM_Y/( max(ps.yC(:)) - min(ps.yC(:)) );
    dx = double((1/mx + 1/my)/2);

%     ps.laser_im_filepath=ps.im_filepath{1};
% 
%     ps.laser_im_filepath(end-6:end-4) = num2str(str2num(ps.laser_im_filepath(end-6:end-4)) - 2);
    laser_im_filepath_list = {laser_results.im_path};
    current_laser_im_filepath = ps.laser_im_filepath;
    current_laser_im_idx = str2num(current_laser_im_filepath(end-6:end-4));
    laser_im_idx_list = [];
    for laser_i=1:numel(laser_im_filepath_list)
        laser_im_idx_list(laser_i) = str2num(laser_im_filepath_list{laser_i}(end-6:end-4));
    end
    
    
    current_laser_im_ind = find(laser_im_idx_list == current_laser_im_idx);
    theta_shift = laser_results(current_laser_im_ind).theta;
    Laser_Cylinder_all = [laser_results.Laser_Cylinder];
    Laser_Cylinder_all(2,:) = Laser_Cylinder_all(2,:) - theta_shift;
    Laser_E_all = convertCylinderToEuclidean(Laser_Cylinder_all,C0,V,Laser_C_ref);
    
%     
%     min_xC = min(ps.xC(:));
%     max_xC = max(ps.xC(:));
%     min_yC = min(ps.yC(:));
%     max_yC = max(ps.yC(:));
%     ind = Laser_E_all(1,:)<min_xC | Laser_E_all(1,:)>max_xC | ...
%         Laser_E_all(2,:)<min_yC | Laser_E_all(2,:)>max_yC |...
%         0; %zC
%     xyC_Laser_all_Mdl = KDTreeSearcher(laser_p_all_after_drop_resolution(1:2,:)');
%     xyC_ps = [ps.xC(:),ps.yC(:)];
% 
%     [xyC_ps_IdxNN,dist] = knnsearch(xyC_Laser_all_Mdl,xyC_ps,'K',1);
%     
%     dist_th = ( mean(mean(ps.xC(:,2:end)-ps.xC(:,1:end-1))) + mean(mean(ps.yC(2:end,:)-ps.yC(1:end-1,:)))) / sqrt(2);
%     
%     Z_measured_ind = zeros(NUM_Y,NUM_X);
%     Z_measured_ind(dist<dist_th) = 1;
%     Z_measured = zeros(NUM_Y,NUM_X);
%     Z_measured(dist<dist_th) = Laser_E_all(3,xyC_ps_IdxNN(dist<dist_th));


    laser_xn_all = Laser_E_all(1:2,:)./ (ones(2,1)*Laser_E_all(3,:));
    laser_p_all = inv_normalize_pixel(laser_xn_all,CamCalib);
    ind = laser_p_all(1,:)<tire_region(1) | laser_p_all(1,:)>tire_region(3) | ...
        laser_p_all(2,:)<tire_region(2) | laser_p_all(2,:)>tire_region(4) |...
        Laser_E_all(3,:)>C0(3);


    Laser_E_all(:,ind) = [];
    laser_p_all(:,ind) = [];
    laser_p_all_after_drop_resolution = (laser_p_all-1)/ps.DROP_RATIO + 1;
    
    temp_p = laser_p_all_after_drop_resolution - ps.region(1:2)'*ones(1,size(laser_p_all_after_drop_resolution,2)) + 1;
    round_temp_p = round(temp_p);
    round_error = sum( (round_temp_p - temp_p).^2);
    
    Z_measured_ind = zeros(NUM_Y,NUM_X);
    Z_measured = zeros(NUM_Y,NUM_X);
    error_map = 2*ones(NUM_Y,NUM_X);
    for laser_p_i = 1:size(laser_p_all,2)
        jj = round_temp_p(1,laser_p_i);
        ii = round_temp_p(2,laser_p_i);
        if (error_map(ii,jj)>round_error(laser_p_i))
            Z_measured_ind(ii,jj) = 1;
            Z_measured(ii,jj) = Laser_E_all(3,laser_p_i);
            error_map(ii,jj) = round_error(laser_p_i);
        end
    end
    
    
    Z = getZfromPQ(p,q,NUM_X,NUM_Y,dx,h1,h2,h3,Z_measured,Z_measured_ind);
    
    if (b_display)
        
        
        laser_im = getImFromPath(ps.laser_im_filepath);
        laser_im = laser_im/quantile(laser_im(:),0.99);
        laser_im(laser_im>1) = 1;
        laser_im_ps_region = laser_im(ps.region_due_to_size_drop(2):ps.DROP_RATIO:ps.region_due_to_size_drop(4),...
            ps.region_due_to_size_drop(1):ps.DROP_RATIO:ps.region_due_to_size_drop(3),:);
        
        
        Laser_X_c = laser_results(current_laser_im_ind).Laser_X_c;
        
        current_xn(1,:) = Laser_X_c(1,:)./Laser_X_c(3,:);
        current_xn(2,:) = Laser_X_c(2,:)./Laser_X_c(3,:);
        current_laser_x_pixel = inv_normalize_pixel(current_xn, CamCalib);
        current_laser_x_pixel_ps_region(1,:) = (current_laser_x_pixel(1,:) - ps.region_due_to_size_drop(1))/ps.DROP_RATIO + 1;
        current_laser_x_pixel_ps_region(2,:) = (current_laser_x_pixel(2,:) - ps.region_due_to_size_drop(2))/ps.DROP_RATIO + 1;
        
        xn(1,:) = Laser_E_all(1,:)./Laser_E_all(3,:);
        xn(2,:) = Laser_E_all(2,:)./Laser_E_all(3,:);
        laser_x_pixel_all = inv_normalize_pixel(xn, CamCalib);
        laser_x_pixel_ps_region_all(1,:) = (laser_x_pixel_all(1,:) - ps.region_due_to_size_drop(1))/ps.DROP_RATIO + 1;
        laser_x_pixel_ps_region_all(2,:) = (laser_x_pixel_all(2,:) - ps.region_due_to_size_drop(2))/ps.DROP_RATIO + 1;

        
        figure();imshow((1-n)/2);
        hold on; scatterPoints(laser_x_pixel_ps_region_all,2,'r');

        figure();
        subplot(1,2,1); imagesc(p); axis equal;
        colorbar;   colormap jet; title('p');caxis([-1,1]);
        hold on; scatterPoints(laser_x_pixel_ps_region_all,2,'r');
        subplot(1,2,2); imagesc(q); axis equal;
        colorbar;   colormap jet; title('q');caxis([-2,2]);
        hold on; scatterPoints(laser_x_pixel_ps_region_all,2,'r');
        
        
        figure();imshow(laser_im_ps_region);
        hold on; scatterPoints(laser_x_pixel_ps_region_all,2,'r');
        hold on; scatterPoints(current_laser_x_pixel_ps_region,1,'g');
        
        figure();imshow(rerenderIm/quantile(rerenderIm(:),0.99));        

        
%         figure();
%         set(gcf,'Name','normalmap with laser measurement');
%         subplot(1,2,1);imshow((1-n)/2);axis on;
% 
%         subplot(1,2,2);imagesc(Z_measured);
%         axis equal;caxis([min(Z_measured(Z_measured>0)),max(Z_measured(:))]);colorbar;
        
        figure();surf(ps.xC,ps.yC,Z,rerenderIm);view(90,0);shading interp;axis equal;
        hold on;scatter3(Laser_E_all(1,:),Laser_E_all(2,:),Laser_E_all(3,:),5,'r');

        
        laser_points_cylinder = convertEuclideanToCylinder...
            (Laser_E_all, C0, V, Laser_C_ref);
        
        
        points_cylinder = convertEuclideanToCylinder...
            ([ps.xC(:), ps.yC(:), Z(:)], C0, V, Laser_C_ref);
        pc_tire = pointCloud(points_cylinder');
        pc_tire_downsize = pcdownsample(pc_tire,'gridAverage',0.1);
        figure();pcshow([pc_tire_downsize.Location(:,2)/180*pi*340,...
            pc_tire_downsize.Location(:,3), pc_tire_downsize.Location(:,1)]);
        hold on; scatter3(laser_points_cylinder(2,:)/180*pi*340,...
            laser_points_cylinder(3,:),laser_points_cylinder(1,:),5,'r');
    end

    
    if (b_display)
        temp_Z_measured_ind = zeros(NUM_Y,NUM_X);
        temp_ind = find(Z_measured_ind==1);
        temp_Z_measured_ind(temp_ind(1)) = 1;
        temp_Z_measured = temp_Z_measured_ind.*Z_measured;
        calculatedZ = getZfromPQ(p,q,NUM_X,NUM_Y,dx,1,1,1,temp_Z_measured,temp_Z_measured_ind);
%         calculatedZ_2 = getZfromPQ(p,q,NUM_X,NUM_Y,dx,h1/10,h2,h3,Z_measured,Z_measured_ind);
%         calculatedZ_3 = getZfromPQ(p,q,NUM_X,NUM_Y,dx,h1*10,h2,h3,Z_measured,Z_measured_ind);
        figure();surf(ps.xC,ps.yC,calculatedZ,rerenderIm);view(90,0);shading interp;axis equal;
%         hold on;scatter3(Laser_E_all(1,:),Laser_E_all(2,:),Laser_E_all(3,:),5,'r');
% 
%         figure();surf(ps.xC,ps.yC,calculatedZ_2,rerenderIm);view(90,0);shading interp;axis equal;
%         hold on;scatter3(Laser_E_all(1,:),Laser_E_all(2,:),Laser_E_all(3,:),5,'r');
% 
%         figure();surf(ps.xC,ps.yC,calculatedZ_3,rerenderIm);view(90,0);shading interp;axis equal;
%         hold on;scatter3(Laser_E_all(1,:),Laser_E_all(2,:),Laser_E_all(3,:),5,'r');% 

        %% convert xyzC to cylinder coordinate system.
        laser_points_cylinder = convertEuclideanToCylinder...
            (Laser_E_all, C0, V, Laser_C_ref);

        points_cylinder = convertEuclideanToCylinder...
            ([ps.xC(:), ps.yC(:), calculatedZ(:)], C0, V, Laser_C_ref);
        pc_tire = pointCloud(points_cylinder');
        pc_tire_downsize = pcdownsample(pc_tire,'gridAverage',0.1);
        figure();pcshow([pc_tire_downsize.Location(:,2)/180*pi*340,...
            pc_tire_downsize.Location(:,3), pc_tire_downsize.Location(:,1)]);
%         hold on; scatter3(laser_points_cylinder(2,:)/180*pi*340,...
%             laser_points_cylinder(3,:),laser_points_cylinder(1,:),5,'r');

%         points_cylinder = convertEuclideanToCylinder...
%             ([ps.xC(:), ps.yC(:), calculatedZ_2(:)], C0, V, Laser_C_ref);
%         pc_tire = pointCloud(points_cylinder');
%         pc_tire_downsize = pcdownsample(pc_tire,'gridAverage',0.1);
%         figure();pcshow([pc_tire_downsize.Location(:,2)/180*pi*340,...
%             pc_tire_downsize.Location(:,3), pc_tire_downsize.Location(:,1)]);
%         hold on; scatter3(laser_points_cylinder(2,:)/180*pi*340,...
%             laser_points_cylinder(3,:),laser_points_cylinder(1,:),5,'r');
% 
%         points_cylinder = convertEuclideanToCylinder...
%             ([ps.xC(:), ps.yC(:), calculatedZ_3(:)], C0, V, Laser_C_ref);
%         pc_tire = pointCloud(points_cylinder');
%         pc_tire_downsize = pcdownsample(pc_tire,'gridAverage',0.1);
%         figure();pcshow([pc_tire_downsize.Location(:,2)/180*pi*340,...
%             pc_tire_downsize.Location(:,3), pc_tire_downsize.Location(:,1)]);
%         hold on; scatter3(laser_points_cylinder(2,:)/180*pi*340,...
%             laser_points_cylinder(3,:),laser_points_cylinder(1,:),5,'r');
% 
% 
    end











end