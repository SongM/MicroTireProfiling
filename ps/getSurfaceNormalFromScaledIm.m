function [ps_res, figure_count] = getSurfaceNormalFromScaledIm(ps_res, ...
    param, CamCalib, LightPos, figure_count)

    %% 3.A.2. PS
    task = '3.A.2.Loading ps image and light strength';
    ps = ps_res;
    scaled_im = ps.scaled_im;

    clear ps_res;
    
    % param
    color_i = 1;

    NUM_Y = size(scaled_im,1);
    NUM_X = size(scaled_im,2);
    NUM_IMAGES = param.N_lights;
    N_lights = param.N_lights;
    % colormap
    map16 = jet(16);

    %%%%%%%%
    % Shadow map
    %%%%%%%%
    % Here use simple thresholding (optional)
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
    
    color_vec_map = [];
    rerenderIm = [];
    for pixel_i = 1:size(scaled_im,1)
        for pixel_j = 1:size(scaled_im,2)
            invalidMapEst_pixel = squeeze(invalidMapEst(pixel_i, pixel_j, :));
            imIrrad_pixel = squeeze(scaled_im(pixel_i, pixel_j, :, :));
            imIrrad_pixel_valid = imIrrad_pixel(:,~invalidMapEst_pixel);
            [color_vec_pixel, int_pixel]=normalizeColVector(imIrrad_pixel_valid);
            color_vec_map(pixel_i, pixel_j,:) = mean(color_vec_pixel,2)/norm(mean(color_vec_pixel,2));
            rerenderIm(pixel_i,pixel_j,:) = color_vec_map(pixel_i, pixel_j,:) * mean(int_pixel);
        end
    end
    
    figure();imshow(color_vec_map);title('color_vec_map');
    figure();imshow(rerenderIm);title('rerenderIm');
    
%% get surface normal (PS)
    est_xyzC(:,:,1) = ps.xC;
    est_xyzC(:,:,2) = ps.xC;
    est_xyzC(:,:,3) = ps.zC_tire_original_est;
    est_xyzC = repmat(est_xyzC,[1,1,1,NUM_IMAGES]);

    
    temp_LightPos = [];
    temp_LightPos(1,1,:,:) = LightPos(:,1:NUM_IMAGES);
    LightPos_rep = repmat(temp_LightPos, [size(est_xyzC,1),size(est_xyzC,2),1,1]);
    LightVector = LightPos_rep - est_xyzC;
    LightVector_length = sqrt(sum(LightVector.^2,3));
    LightUnitVector = LightVector./repmat(LightVector_length,[1,1,3,1]);

    % params
    MIN_NUM_IM = 4;
    IMNOISE_STD = 0.1;
    % To = 2.5;
    Tm = (N_lights+1)*(IMNOISE_STD^2);

    % input:    
    %     LightUnitVector: NUM_Y x NUM_X x 3 x NUM_IMAGES
    %     invalidMapEst: NUM_Y x NUM_X x NUM_IMAGES
    %     scaled_im: NUM_Y x NUM_X x 3 x NUM_IMAGES

    % output:
    %   n:  NUM_Y x NUM_X x 3
    %   rho: NUM_Y x NUM_X

    n = nan(NUM_Y,NUM_X,3);
    rho = nan(NUM_Y,NUM_X);
    residualMap = nan(NUM_Y,NUM_X,NUM_IMAGES); % NUM_Y x NUM_X x NUM_IMAGES
    outlierMap = nan(NUM_Y,NUM_X,NUM_IMAGES); % NUM_Y x NUM_X x NUM_IMAGES
    
    for pixel_i = 1:size(scaled_im,1)
        for pixel_j = 1:size(scaled_im,2)
            invalidMap_ind = squeeze(invalidMapEst(pixel_i, pixel_j, :));
            imIrrad_ind = squeeze(scaled_im(pixel_i, pixel_j, :, :));
            lC_ind = squeeze(LightUnitVector(pixel_i, pixel_j, :, :));
            
            t_flag = 0;
            
            while (t_flag==0)
                if (sum(~invalidMap_ind)<MIN_NUM_IM)
                    outlierMap(pixel_i,pixel_j,:) = 1;
%                     if(b_debug)
%                         fprintf(2,['ii=',num2str(ii),', N_valid=',num2str(sum(~invalidMap_ind)),'.\n']);
%                     end
                    t_flag = 1;
                else
                    t_imIrrad_ind = imIrrad_ind(color_i,~invalidMap_ind)';
                    %             [~,max_ind] = max(t_imIrrad_ind);
                    %             valid_inds = find(invalidMap_ind==0);
                    %             t_ind = valid_inds(max_ind);
                    %             invalidMap_ind(t_ind) = 1;
                    %             t_imIrrad_ind = imIrrad_ind(color_i,~invalidMap_ind)';            
                    t_lC_ind = lC_ind(:,~invalidMap_ind)';
                    t_lC_ind_inv = pinv(t_lC_ind);
                    scaledNormal = t_lC_ind_inv*t_imIrrad_ind;
                    t_h = diag(t_lC_ind*t_lC_ind_inv);
                    t_e = t_imIrrad_ind - t_lC_ind*scaledNormal;
                    t_mse = mean(t_e.^2);
                    t_mse_ratio = t_mse/(Tm*mean(abs(t_imIrrad_ind))^2);
                    if (t_mse_ratio > 1)
%                         if(b_debug)
%                             fprintf(['ii=',num2str(ii),', N_valid=',num2str(sum(~invalidMap_ind)),', t_mse_ratio=',num2str(t_mse_ratio),', large.\n']);
%                         end
                        [~,tt_ind] = max(abs(t_e));
                        valid_inds = find(invalidMap_ind==0);
                        t_ind = valid_inds(tt_ind);
                        invalidMap_ind(t_ind) = 1;
                    else
                        [n(pixel_i,pixel_j,:),rho(pixel_i,pixel_j)] = normalizeColVector(scaledNormal);
                        t_residualMap = nan(1,NUM_IMAGES);
                        t_residualMap(~invalidMap_ind) = t_e';
                        residualMap(pixel_i,pixel_j,:) = t_residualMap;
                        outlierMap(pixel_i,pixel_j,:) = invalidMap_ind;
%                         if(b_debug)
%                             fprintf(['ii=',num2str(ii),', N_valid=',num2str(sum(~invalidMap_ind)),...
%                                 ', t_mse_ratio=',num2str(t_mse/(Tm*mean(t_imIrrad_ind)^2)),'.\n']);
%                         end


                        t_flag = 1;
                    end
                end
            end
        end
    end
    ps.color_i = color_i;
    ps.color_vec_map = color_vec_map;
    ps.rerenderIm = rerenderIm;
    ps.n = n;
    ps.rho = rho;
    ps.residualMap = residualMap;
    ps.outlierMap = outlierMap;
    ps_res = ps;
    
%% display result    
    map100 = jet(100);
    map100(1,:) = 1;
    map100(end,:) = 0;

    x=-1:0.01:1;
    y=-1:0.01:1;
    [xx,yy] = meshgrid(x,y);
    zz = 1 - xx.^2 - yy.^2;
    zz(zz<0) = 1;
    zz = sqrt(zz);
    figure();
    subplot(1,2,1); imagesc(acosd(-n(:,:,3))); axis equal; caxis([-2,90]);
    title('direction1, indicate the slope of the surface change');
    subplot(1,2,2); imagesc( acosd(zz) ); colorbar; axis equal; axis off; caxis([-2,90]);
    colormap(map100); 
    title({'direction1 indicator, direction of a semisphere', 'indicates the angle between the surface normal and the vector which perpendicular to the image'});

    figure();
    subplot(1,2,1); imagesc( atan2(-n(:,:,1),-n(:,:,2))*180/pi ); 
    axis equal; caxis([-190,180]);
    title('direction2, indicate the direction of the surface change');
    subplot(1,2,2); imagesc( atan2(-xx,-yy)*180/pi  .* (zz<0.9) );
    colorbar; axis off; axis equal; caxis([-190,180]);
    colormap(map100); 

    figure();imagesc(sum(outlierMap,3));axis equal;
    colorbar;caxis([0,16]);colormap(map16);
    title('outlierMap');
    figure();set(gcf, 'Name', 'residualMap');
    for i=1:16
        subplot(4,4,i);imagesc(residualMap(:,:,i));colorbar,axis equal;colormap(map100),caxis([-0.01,0.1])
    end
    
    figure();set(gcf, 'Name', 'scaled_im,color_i');
    for i=1:16
        subplot(4,4,i);imagesc(scaled_im(:,:,color_i,i));colorbar,axis equal;colormap(map100),caxis([-0.1,1]);
    end
end
