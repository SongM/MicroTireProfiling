function [ps_res, figure_count] = getScaledImFromPSFileList(ps_file_list, folder_path, ...
    param, CamCalib, LightPos, lightStrength_mat_file, V, C0, DROP_RATIO, figure_count)

    tire_region = param.ps.tire_region;
    whiteboard_region = param.caliblightillu.whiteboard_region;

    N_lights = param.N_lights;
%     s_RGB = param.ps.s_RGB;
%     NUM_IMAGES = param.N_lights;

%     map16 = jet(16);


    %% 3.A.1.Loading ps image and light strength
    task = '3.A.1.Loading ps image and light strength';
    displayTask(task,3);
    
    
    % 3.A.1.i load ps image
    ps_raw_temp_mat_file = [folder_path,'\ps_raw_temp.mat'];
    if(1)
        [ps_raw, figure_count] = loadPSIm...
            ([folder_path,'\',param.file_name.laserandps], ps_file_list.im_file_list, N_lights, tire_region, figure_count);
        ps_raw.laser_im_filepath = [folder_path,'\',param.file_name.laserandps,'\',ps_file_list.laser_im_file.name]; 
        save(ps_raw_temp_mat_file,'ps_raw','-v7.3');
    else
        load(ps_raw_temp_mat_file,'ps_raw');
    end
        
 
    % 3.A.1.ii load light strength
    load(lightStrength_mat_file);
    ps_raw.region = tire_region;

    R = 332;
    
    xyzC_Cylinder = getCylinderxyzC(R,V,C0,CamCalib,tire_region);
    ps_raw.xC = xyzC_Cylinder(:,:,1);
    ps_raw.yC = xyzC_Cylinder(:,:,2);
    ps_raw.zC_tire_original_est = xyzC_Cylinder(:,:,3);

    clear xyzC_Cylinder
    
    ps_raw.whiteboard_region = whiteboard_region;

    disp('load light strength done.');
    ps_raw.original_region = tire_region;
    ps_raw.region_due_to_size_drop = tire_region;
    
    % 3.A.1.iii Drop process size for ps
    ps = dropProcessingSizeForPS(ps_raw, DROP_RATIO);
    disp('Drop process size for ps done.');
    clear ps_raw
    
    figure_count = displaySurfaceAndLightPos(ps.xC, ps.yC, ps.zC_tire_original_est, LightPos, CamCalib, figure_count);

    [NUM_Y, NUM_X] = size(ps.xC);
    
    % 3.A.1.v get scaled intensity for ps
    xyzC_Cylinder(:,:,1) = ps.xC;
    xyzC_Cylinder(:,:,2) = ps.yC;
    xyzC_Cylinder(:,:,3) = ps.zC_tire_original_est;
    
    [~, RiC_Cylinder] = getUnitLight(LightPos, xyzC_Cylinder);
    RiC_Cylinder = reshape(RiC_Cylinder,NUM_Y,NUM_X,N_lights);         
    
    xyzC_whiteboard_reshaped = reshape(xyzC_whiteboard,numel(xyzC_whiteboard)/3,3);
    xyzC_Cylinder_reshaped = reshape(xyzC_Cylinder,numel(xyzC_Cylinder)/3,3);
    xyC_whiteboard_reshaped = xyzC_whiteboard_reshaped(:,1:2);
    xyC_Cylinder_reshaped = xyzC_Cylinder_reshaped(:,1:2);
    clear xyzC_whiteboard_reshaped xyzC_Cylinder_reshaped;

    % lightStrength_adjusted_due_to_the_difference_between xyzC_tire and
    % xyzC_whiteboard;

    xyC_whiteboard_Mdl = KDTreeSearcher(xyC_whiteboard_reshaped);
    xyC_Cylinder_IdxNN = knnsearch(xyC_whiteboard_Mdl,xyC_Cylinder_reshaped,'K',1);
    lightStrength_reshaped = reshape(lightStrength,numel(lightStrength)/3/N_lights,3,N_lights);
    lightStrength_adjusted_reshaped = lightStrength_reshaped(xyC_Cylinder_IdxNN,:,:);
    lightStrength_adjusted = reshape(lightStrength_adjusted_reshaped, NUM_Y, NUM_X, 3, N_lights);
    
    ps.RiC_Cylinder = RiC_Cylinder;
    ps.lightStrength_adjusted = lightStrength_adjusted;

    
    RiC_Cylinder_sqr_repeated_for_3_color(:,:,1,:) = RiC_Cylinder.^2;
    RiC_Cylinder_sqr_repeated_for_3_color(:,:,2,:) = RiC_Cylinder.^2;
    RiC_Cylinder_sqr_repeated_for_3_color(:,:,3,:) = RiC_Cylinder.^2;
%     
    temp_ps_im = ps.im_raw.* RiC_Cylinder_sqr_repeated_for_3_color ./ lightStrength_adjusted;
%     temp_ps_im = ps.im_raw.* RiC_Cylinder_sqr_repeated_for_3_color ./ repmat(lightStrength_adjusted(:,:,2,:),1,1,3,1);
    

    % 3.A.1.vi rescale image irradiance to [0,1]
%     SAT_PERCENT = 0.0001;    % proportion of saturated pixels (disregard too bright pixels)
    ps.SAT_PERCENT = 0.0001;
    ps.scaled_im = temp_ps_im/quantile(temp_ps_im(:),1-ps.SAT_PERCENT);
    fprintf(['SAT_PERCENT = ', num2str(ps.SAT_PERCENT),...
        ', capped_intensity = ', num2str( quantile(temp_ps_im(:),1-ps.SAT_PERCENT) ),...
        ', rescale image irradiance down.\n']);
    clear SAT_PERCENT;
    disp('get scaled im done');
    ps_res = ps;
% %     job_done = [job_done;task];
%     clear task
    
end

function [ps_raw, figure_count] = loadPSIm...
    (folder_name, ps_im_file_list, N_lights, region, figure_count)

%     ps_raw.laser_im_filepath = [folder_name,'\',ps_im_file_list(1).name]; 
    ims = [];
    fprintf('loading raw ps_im:16/');
    
    for ps_i = 1:N_lights
        fprintf(['...',num2str(ps_i)]);
        ps_raw.im_filepath{ps_i} = [folder_name,'\',ps_im_file_list(ps_i).name];
        [ps_im, max_val] = getImFromPath(ps_raw.im_filepath{ps_i});
        ps_im = double(ps_im).*max_val;
        ps_im = ps_im(region(2):region(4),region(1):region(3),:);
        % just to check if the brightest point is on the surface of tire
        fprintf(['.(',num2str(max_val),',',num2str(max(ps_im(:))),')']);
        ims(:,:,:,ps_i) = ps_im;
    end
    ps_raw.region = region;

    ps_raw.im_raw = ims;
    
    % disp intensity
    figure_count = figure_count + 1;
    temp_scaled_factor = quantile(ims(:),0.99);
    for ps_i = 1:N_lights
        figure(figure_count);
        subplot(ceil(N_lights/ceil(sqrt(N_lights))),ceil(sqrt(N_lights)),ps_i);
        imshow( ims(:,:,:,ps_i) / temp_scaled_factor );
    end
    disp('...ps_scan loading done.');
end

% function Z = estimateTireOriginalZ(param, X, Y)
% 
%     R = param.tire.original_est.R;
%     C = param.tire.original_est.C;
%     V = param.tire.original_est.V;
%     
%     
%     % |(P-C) - V((P-C)'V)| = R
%     % P-C = Pr
%     % (P-C)'V = Xr*u + Yr*v + Zr*w 
%     %         = Xr*u + Yr*v + (Z-C3)*w 
%     %         = w*Z + Xru_Yrv_minus_wC3
%     %         = w*Z + tmp
%     
%     % |Xr - u*(w*Z + tmp)  | = |Xr - u*tmp + (-u*w)Z| = |a1 + b1*Z|
%     % |Yr - v*(w*Z + tmp)  | = |Yr - v*tmp + (-v*w)Z| = |a2 + b2*Z| 
%     % |Z-C3 - w*(w*Z + tmp)| = |-C3-w*tmp + (1-w^2)Z| = |a3 + b3*Z|
% 
%     % (a1+b1*Z)^2 + (a2+b2*Z)^2 + (a3+b3*Z)^2 - R^2 = 0;
%     % (b1^2+b2^2+b3^2)*Z^2 + 2(a1b1+a2b2+a3b3)*Z +(a1^2+a2^2+a3^2-R^2) = 0;
%     % c1*Z^2 + 2c2*Z + c3 = 0;
%     % (c1*Z + c2)^2 = c2^2-c1*c3
%     % Z = (-c2 - sqrt(c2^2-c1*c3))/c1
% 
%     X_remain = X - C(1);
%     Y_remain = Y - C(2);
% 
%     tmp = X_remain*V(1) + Y_remain*V(2) - V(3)*C(3);
%     
%     a1 = X_remain - V(1)*tmp;
%     a2 = Y_remain - V(2)*tmp;
%     a3 = -C(3) - V(3)*tmp;
%     
%     b1 = -V(1)*V(3);
%     b2 = -V(2)*V(3);
%     b3 = 1-V(3)^2;
%     
%     c1 = b1.^2+b2.^2+b3.^2;
%     c2 = a1.*b1+a2.*b2+a3.*b3;
%     c3 = a1.^2+a2.^2+a3.^2-R.^2;
%     Z = (-c2 - sqrt(c2.^2-c1.*c3)) ./ c1;
%     
% end

function xyzC = getCylinderxyzC(R,V,C0,CamCalib,region)
%     
% P = l*Z - C0;
% (P-(P'*V)*V)' * (P-(P'*V)*V) - R^2 = 0;
% 
% P'*P - 2*(P'*V)*P'*V + [(P'*V)^2]*V'*V - R^2 = 0;
% P'*P - (P'*V)^2 - R^2 = 0;
% 
% (l'*l)*Z^2 - 2*(l'*C0)*Z + C0'*C0
%     - ( (l'*V)*Z - (C0'*V) )^2 - R^2 = 0;
% 
% [ (l'*l) - (l'*V)^2 ] *Z^2 
%     + 2*[-(l'*C0) + (l'*V)*(C0'*V)] *Z
%     + [C0'*C0 - (C0'*V)^2 - R^2] = 0;
%     
% Z = ( b - sqrt(b^2-4ac) ) / 2a
    

    x = region(1):region(3);
    y = region(2):region(4);
    NUM_Y = numel(y);
    NUM_X = numel(x);
    [xx,yy] = meshgrid(x,y);
    l = normalize_pixel([xx(:),yy(:)]',CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
    l(3,:) = 1;
    a = (dot(l,l)' - (l'*V).^2);
    b = 2*(-l'*C0 + (l'*V) * (C0'*V));
    c = C0'*C0 - (C0'*V)^2 - R^2;
    
    Z = (-b-sqrt(b.^2-4*a*c))./(2*a);
    
    xyzC(:,:,1) = reshape(Z'.*l(1,:),NUM_Y,NUM_X);
    xyzC(:,:,2) = reshape(Z'.*l(2,:),NUM_Y,NUM_X);
    xyzC(:,:,3) = reshape(Z',NUM_Y,NUM_X);
end

