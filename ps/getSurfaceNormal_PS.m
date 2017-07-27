function getSurfaceNormal_PS(mat_file)
    load(mat_file);
    tire_region = param.ps.tire_region;
    
    N_lights = param.N_lights;
    s_RGB = param.ps.s_RGB;
    NUM_IMAGES = param.N_lights;
    DROP_RATIO = 5;

    %     map16 = flipud(hot(16))*3;
    %     map16(13:16,:) = 0;
    %     map16(map16>1) = 1;
    map16 = jet(16);


    %% 3.A.1.Loading ps image and light strength
    task = '3.A.1.Loading ps image and light strength';
    displayTask(task,3);
    
    
    % 3.A.1.i load ps image
    folder_name = [folder_path,'\dng'];
    ps_file_list = dir([folder_name,'\*',param.file_name.ps,'*mat']);
    [ps_raw, figure_count] = loadPSIm...
        (folder_name, ps_file_list, N_lights, tire_region, figure_count);

    % 3.A.1.ii load light strength
%     fprintf(2,'param.caliblightillu.whiteboard_region is assigned, remove this line later');
%     param.caliblightillu.whiteboard_region = [1121,2343,3950,5015];
    load(lightStrength_mat_file);
    ps_raw.lightStrength = lightStrength;
    ps_raw.xyzC_whiteboard = xyzC_whiteboard;
    ps_raw.whiteboard_region = param.caliblightillu.whiteboard_region;

    % notice this xC should comes from tire_region, instaed of
    % whiteboard_region.
    
    ps_raw.xC = xyzC_whiteboard(:,:,1);
    ps_raw.yC = xyzC_whiteboard(:,:,2);

    clear lightStrength xyzC_whiteboard
    disp('load light strength done.');
    ps_raw.original_region = tire_region;
    ps_raw.region_due_to_size_drop = tire_region;
    ps_raw.region = tire_region;
    
    % 3.A.1.iii Drop process size for ps
%     load('temp.mat');
    ps = dropProcessingSizeForPS(ps_raw, DROP_RATIO);
    disp('Drop process size for ps done.');
    clear ps_raw
    
    
    
    % 3.A.1.iv get original tire surface guess
    param.tire.original_est.R = 332;
    % param.tire.original_est.C = RotateCalib.cis(:,1);
    % param.tire.original_est.V = RotateCalib.V;
%     param.tire.original_est.C = [0; 0; param.tire.original_est.R+mean(mean(ps.xyzC_whiteboard(:,:,3)))+10];
%     param.tire.original_est.V = [0;1;0];
    param.tire.original_est.C = [0.0064;-0.1473;1]*(param.tire.original_est.R + mean(mean(ps.xyzC_whiteboard(:,:,3))) + 10);
    param.tire.original_est.C = [8.3;0;1198.1];
    param.tire.original_est.V = [0.0133;0.9998;0.0159];

%     param.tire.original_est.C = C0;
%     param.tire.original_est.V = V;
    

    
    
    ps.zC_tire_original_est = estimateTireOriginalZ(param, ps.xC, ps.yC);
    
%     xyzC_tire_original_est = initTireSurface(param,CamCalib);
%     figure_count = displaySurfaceAndLightPos(ps.xC, ps.yC, ps.zC_tire_original_est, ps.xyzC_whiteboard, LightPos, CamCalib, figure_count);
    figure_count = displaySurfaceAndLightPos(ps.xC, ps.yC, ps.zC_tire_original_est, [], LightPos, CamCalib, figure_count);

    
    
    % 3.A.1.v get scaled intensity for ps
%     ps_raw = getScaledIntensityForPS;
    R_tire = sqrt(ps.xC.^2 + ps.yC.^2 + ps.zC_tire_original_est.^2);
    % lightStrength_adjusted_due_to_the_difference_between xyzC_tire and
    % xyzC_whiteboard;
    start_point = ps.region_due_to_size_drop - ps.original_region+1;
    DROP_RATIO = ps.DROP_RATIO;
    lightStrenth_adjusted = ps.lightStrength( start_point(1):DROP_RATIO:size(ps.lightStrength,1), start_point(2):DROP_RATIO:size(ps.lightStrength,2), :, :); 
    ps.im = ps.im_raw.*repmat(R_tire,1,1,3,N_lights)./repmat(lightStrenth_adjusted(:,:,2,:),1,1,3,1);
    % ps.im = ims.*repmat(R_tire,1,1,3,N_lights) ./ lightStrength;
    disp('get scaled im done');
    

    NUM_X = size(ps.im,2);
    NUM_Y = size(ps.im,1);

    % 3.A.1.vi rescale image irradiance to [0,1]
    SAT_PERCENT = 0.0001;    % proportion of saturated pixels (disregard too bright pixels)
    ps.scaled_im = ps.im/quantile(ps.im(:),1-SAT_PERCENT);
    fprintf(['SAT_PERCENT = ', num2str(SAT_PERCENT),...
        ', capped_intensity = ', num2str( quantile(ps.im(:),1-SAT_PERCENT) ),...
        ', rescale image irradiance down.\n']);
    clear SAT_PERCENT;
    
    job_done = [job_done;task];
    clear task
    
        %% 3.A.2. PS
    task = '3.A.2.Loading ps image and light strength';

    %%%%%%%%
    % Shadow map
    %%%%%%%%
    % Here use simple thresholding (optional)
    scaled_im = ps.scaled_im;
    SHADOW_THRES = 0.03; % lower limit of valid image irradiance
    shadowMapEst = scaled_im(:,:,1,:) <= SHADOW_THRES & scaled_im(:,:,2,:) <= SHADOW_THRES & scaled_im(:,:,3,:) <= SHADOW_THRES;
    shadowMapEst = squeeze(shadowMapEst);
    figure();imagesc(sum(shadowMapEst,3));axis equal;
    colorbar;caxis([0,16]);colormap(map16);
    title('shadowMapEst, 0->no shadow influence, larger->high shadow influence');
    % % Use proposed shadow detection

    %%%%%%%%%
    % Saturation map
    %%%%%%%%%
    satMapEst = scaled_im(:,:,1,:) >= 1 | scaled_im(:,:,2,:) >= 1 | scaled_im(:,:,3,:) >= 1;  % If one channel is saturated, then not usable
    satMapEst = squeeze(satMapEst);
    figure();imagesc(sum(satMapEst,3));axis equal;
    colorbar;caxis([0,16]);colormap(map16);
    title('satMapEst, 0->no saturation influence, larger->high saturation influence');
    invalidMapEst = shadowMapEst | satMapEst;

    %% RPCA to find unit diffuse color
    Td = 0.035; % Threshold for robust PCA in detecting outliers
    % Td = 0.05;
    tic
    [dRGBEst, outlierMap, mean_color_error] = pcaDiffColor(scaled_im, invalidMapEst, Td);
    figure();imshow(dRGBEst);title('dRGBEst');
    figure();imagesc(sum(outlierMap,3));axis equal;
    colorbar;caxis([0,16]);colormap(map16);
    title('outlierMap, num of outliers');
    figure();imagesc(mean_color_error);colorbar;
    colormap jet;axis equal;title('Mean Color Error');
    % [dRGBEst, ~] = pcaDiffColor(im, maskRegion, invalidMapEst, Td);
    toc


    %% Transform RGB to SUV
    Tc = 1;    % chromatic angle in degrees, > Tc, then SUV transformable, otherwise, not
    rgb2suv;
    figure_count = figure_count+1;
    figure(figure_count);set(figure_count,'Name','fsEst');
    for i=1:NUM_IMAGES
        subplot(4,4,i);
        imagesc(fsEst(:,:,i));colorbar;axis equal;colormap colorcube;
    end

    %% Estimate unit light directions (assume surface deviation is small compared with light-surface distance in z)
    xyzRefC(:,:,1) = ps.xC;
    xyzRefC(:,:,2) = ps.xC;
    xyzRefC(:,:,3) = ps.zC_tire_original_est;
    xyzRefC_rep = repmat(xyzRefC,[1,1,1,NUM_IMAGES]);
    temp_LightPos = [];
    temp_LightPos(1,1,:,:) = LightPos(:,1:NUM_IMAGES);
    LightPos_rep = repmat(temp_LightPos, [size(xyzRefC_rep,1),size(xyzRefC_rep,2),1,1]);
    LightVector = LightPos_rep - xyzRefC_rep;
    LightVector_length = sqrt(sum(LightVector.^2,3));
    LightUnitVector = LightVector./repmat(LightVector_length,[1,1,3,1]);

    
    %% PS in UV
    IMNOISE_STD = 0.1; % image noise standard deviation estimate, assume white additive gaussian
    MIN_NUM_IM = 4;
    outlierParaPsUV.To = 2.5;    % tolerance of studentized residual std, the bigger, the more tolerant, default value is 2.5
    outlierParaPsUV.Tm = 9*(IMNOISE_STD^2); % mean square error tolerance, smaller than this one, define tolerable
    outlierParaPsUV.MIN_NUM_IM = MIN_NUM_IM; % minimum number of image for PS required
    lC=permute(reshape(LightUnitVector,NUM_X*NUM_Y,3,NUM_IMAGES),[2,1,3]);
    psUV;
    figure();subplot(1,2,1);
    imshow((-nCRecon+1)/2);title('nCRecon');
    subplot(1,2,2);display_surface_normal_map;title('surface normal color indicator');


    %% Clear bad fittings
    % BAD_FIT_PERCENT = 0.01;    % bad fitting proportion (residual)
    % MAX_SLOPE = 89; % in degrees (reject sharp changes)
    BAD_FIT_PERCENT = 0.03;    % bad fitting proportion (residual)
    MAX_SLOPE = 85; % in degrees (reject sharp changes)
    clearBadFit;
    figure();
    imshow((-nCReconClear+1)/2);title('nCReconClear');
    figure();
    subplot(1,2,1); imagesc(nCRecon(:,:,1)./nCRecon(:,:,3)); colorbar; axis equal; caxis([-4,4]); colormap jet;title('p');
    subplot(1,2,2); imagesc(nCRecon(:,:,2)./nCRecon(:,:,3)); colorbar; axis equal; caxis([-4,4]); colormap jet;title('q');
    set(gcf,'Name', 'nCRecon');

    figure();
    subplot(1,2,1); imagesc(nCReconClear(:,:,1)./nCReconClear(:,:,3)); colorbar; axis equal; caxis([-4,4]); colormap jet; title('p');
    subplot(1,2,2); imagesc(nCReconClear(:,:,2)./nCReconClear(:,:,3)); colorbar; axis equal; caxis([-4,4]); colormap jet; title('q');
    set(gcf,'Name', 'nCReconClear');
    figure();
    imagesc(kdClear);colorbar;axis equal;
    title('kdClear');

    fprintf('Relighting the scene...');
    vC = [0;0;-1]*ones(1,NUM_Y*NUM_X);    % Used as light direction
    rerenderIm = max(0,dot(vC, reshape(nCReconClear, [NUM_Y*NUM_X, 3])', 1));
    rerenderIm = reshape(dRGBEst, [NUM_Y*NUM_X, 3])' .* repmat(kdClear(:)',[3,1]) .* repmat(rerenderIm,[3,1]);
    rerenderIm = reshape(rerenderIm', [NUM_Y, NUM_X, 3]);
    fprintf('done!\n');
    figure();imshow(rerenderIm/max(rerenderIm(:)));
    title('rerenderIm');



    
    job_done = [job_done;task];
    clear task

    SNI_SONG;
    temp_show_ps_result;
end

function [ps_raw, figure_count] = loadPSIm...
    (folder_name, ps_file_list, N_lights, region, figure_count)

    ps_raw.laser_im_filepath = [folder_name,'\',ps_file_list(1).name]; 
    ims = [];
    fprintf('loading raw ps_im:');
    
    for ps_i = 1:N_lights
        fprintf(['...',num2str(ps_i)]);
        ps_raw.im_filepath{ps_i} = [folder_name,'\',ps_file_list(ps_i+1).name];
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

function Z = estimateTireOriginalZ(param, X, Y)

    R = param.tire.original_est.R;
    C = param.tire.original_est.C;
    V = param.tire.original_est.V;
    
    
    % |(P-C) - V((P-C)'V)| = R
    % P-C = Pr
    % (P-C)'V = Xr*u + Yr*v + Zr*w 
    %         = Xr*u + Yr*v + (Z-C3)*w 
    %         = w*Z + Xru_Yrv_minus_wC3
    %         = w*Z + tmp
    
    % |Xr - u*(w*Z + tmp)  | = |Xr - u*tmp + (-u*w)Z| = |a1 + b1*Z|
    % |Yr - v*(w*Z + tmp)  | = |Yr - v*tmp + (-v*w)Z| = |a2 + b2*Z| 
    % |Z-C3 - w*(w*Z + tmp)| = |-C3-w*tmp + (1-w^2)Z| = |a3 + b3*Z|

    % (a1+b1*Z)^2 + (a2+b2*Z)^2 + (a3+b3*Z)^2 - R^2 = 0;
    % (b1^2+b2^2+b3^2)*Z^2 + 2(a1b1+a2b2+a3b3)*Z +(a1^2+a2^2+a3^2-R^2) = 0;
    % c1*Z^2 + 2c2*Z + c3 = 0;
    % (c1*Z + c2)^2 = c2^2-c1*c3
    % Z = (-c2 - sqrt(c2^2-c1*c3))/c1

    X_remain = X - C(1);
    Y_remain = Y - C(2);

    tmp = X_remain*V(1) + Y_remain*V(2) - V(3)*C(3);
    
    a1 = X_remain - V(1)*tmp;
    a2 = Y_remain - V(2)*tmp;
    a3 = -C(3) - V(3)*tmp;
    
    b1 = -V(1)*V(3);
    b2 = -V(2)*V(3);
    b3 = 1-V(3)^2;
    
    c1 = b1.^2+b2.^2+b3.^2;
    c2 = a1.*b1+a2.*b2+a3.*b3;
    c3 = a1.^2+a2.^2+a3.^2-R.^2;
    Z = (-c2 - sqrt(c2.^2-c1.*c3)) ./ c1;
    
end

