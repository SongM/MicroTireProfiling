function calibLightIllu(mat_file)
    load(mat_file);
        
        

    whiteboard_file_name = param.file_name.caliblightillu;
    whiteboard_region = param.caliblightillu.whiteboard_region;
    filter_N = param.caliblightillu.filter_N;

    b_display = 0;
    
    
    folder_name = [folder_path];
    whiteboard_file_list = dir([folder_name,'\',whiteboard_file_name,'\*JPG']);
%     whiteboard_file_list = dir([folder_name,'\',whiteboard_file_name,'\*DNG']);
        
    t_im = getImFromPath([folder_name,'\',whiteboard_file_name,'\',whiteboard_file_list(2).name]);
    figure();imshow(t_im/quantile(t_im(:),1-1e-2));
    hold on; rectangle('Position',[whiteboard_region(1:2),whiteboard_region(3:4)-whiteboard_region(1:2)],'EdgeColor','r');
    drawnow();
    clear t_im;
    
    [whiteboard_plane_eqn, error, ref_board_calib_light_illu] =...
        getWhiteBoardPlane([folder_name,'\',whiteboard_file_name,'\',whiteboard_file_list(1).name],...
        CamCalib, param, error);
    
    whiteboard.whiteboard_plane_eqn = whiteboard_plane_eqn;
    whiteboard.ref_board_calib_light_illu = ref_board_calib_light_illu;
         
    xyzC_whiteboard = getWhiteBoardxyzC(whiteboard_plane_eqn,CamCalib,whiteboard_region);
    NUM_X = whiteboard_region(3)-whiteboard_region(1)+1;
    NUM_Y = whiteboard_region(4)-whiteboard_region(2)+1;

    nC_whiteboard = getUnitSurfNormC(xyzC_whiteboard);
    [lightDirectionC_whiteboard, RiC_whiteboard] = getUnitLight(LightPos, xyzC_whiteboard);
    [cosThetai_whiteboard, ~] = getSurfIrradSimple(reshape(nC_whiteboard, [NUM_Y*NUM_X, 3])', lightDirectionC_whiteboard, 1);

    N_lights = size(LightPos,2);
    
    
    cosThetai_whiteboard = reshape(cosThetai_whiteboard,NUM_Y,NUM_X,N_lights);
    RiC_whiteboard = reshape(RiC_whiteboard,NUM_Y,NUM_X,N_lights);         

%%
% the light illuminate calibration may get better by considering modelling
% the illumination function of the led.
% the distance will be take into consideration, and the theta and phi angle
% from the led to the lighting surface will form a look-up table, and a 2d
% interperation function can be used to get better lighting condition
% estimation.
%%
% k = (LightPos(:,light_i)'*whiteboard_plane_eqn(1:3)+whiteboard_plane_eqn(4)) ./ (lC(:,:,light_i)'*whiteboard_plane_eqn(1:3));
%         xyzC_on_whiteboard_xC = reshape(LightPos(1,light_i)-k'.*lC(1,:,light_i),NUM_Y,NUM_X)
%         xyzC_on_whiteboard_yC = reshape(LightPos(2,light_i)-k'.*lC(2,:,light_i),NUM_Y,NUM_X)
%     [theta1,theta2] = getAngleFromVector(lc);
%     figure_count = figure_count+1;
%     figure(figure_count);
        
%%        

    fprintf('getLightIllu...');
    for light_i = 1:N_lights
        fprintf([num2str(light_i),'...']);
        whiteboard.im_filepath{light_i} = [folder_name,'\',whiteboard_file_name,'\',whiteboard_file_list(light_i+2).name];
        [wb_im_ind, max_val_ind] = getImFromPath(whiteboard.im_filepath{light_i});
        wb_im_trimed = double( wb_im_ind(whiteboard_region(2):whiteboard_region(4), ...
                whiteboard_region(1):whiteboard_region(3), :) ) * max_val_ind;

            
        cosThetai_whiteboard_ind = cosThetai_whiteboard(:,:,light_i);
        RiC_whiteboard_ind = RiC_whiteboard(:,:,light_i);
        
        for color_i = 1:3
            wb_Intensity_filtered_ind(:,:,color_i) = ...
                imfilter(wb_im_trimed(:,:,color_i),...
                ones(filter_N)/filter_N^2);
            
            lightStrength_ind(:,:,color_i) =...
                wb_Intensity_filtered_ind(:,:,color_i) ./ cosThetai_whiteboard_ind .* (RiC_whiteboard_ind.^2);

            rerendingIm(:,:,color_i) = lightStrength_ind(:,:,color_i) .*cosThetai_whiteboard_ind ./ (RiC_whiteboard_ind.^2);
        end
        rerendering_error = wb_im_trimed - rerendingIm;
        rerendering_error_gray = sqrt(sum(rerendering_error.^2,3));
        error.caliblightillu.whiteboard_rerending_error(light_i) = ...
            mean(rerendering_error_gray(:));

        
        lightStrength(:,:,:,light_i) = lightStrength_ind;
        
        if(b_display)
            figure_count = figure_count+1;
            displayResult(figure_count, wb_im_ind, wb_im_trimed, ...
            wb_Intensity_filtered_ind, lightStrength_ind, ...
            cosThetai_whiteboard_ind, rerendering_error_gray, whiteboard_region)
            set(figure_count,'Name',['light_illu_',num2str(light_i),'_max_val_',num2str(max_val_ind)]);
            drawnow();
        end
    end
    fprintf('\n');
    
    clear whiteboard_file_name tire_region filter_N b_display;
    clear folder_name whiteboard_file_list
    clear NUM_X NUM_Y nC_whiteboard
    clear lightDirectionC_whiteboard cosThetai_whiteboard
    clear N_lights light_i color_i
    clear wb_im_ind max_val_ind wb_im_trimed cosThetai_whiteboard_ind
    clear wb_Intensity_filtered_ind lightStrength_ind rerendingIm 
    clear RiC_whiteboard RiC_whiteboard_ind
    clear rerendering_error rerendering_error_gray
    
    
    lightStrength_mat_file = [folder_path,'\lightStrength_',case_name,'.mat'];
    save(lightStrength_mat_file,'lightStrength','xyzC_whiteboard','-v7.3');
    fprintf('    lightStrength saved to: \n');
    disp(lightStrength_mat_file);
    clear lightStrength xyzC_whiteboard
    c = clock;
    disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
    timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];

    individual_mat_file.caliblightillu = [data_folder, '\', case_name, '\processed_data\', 'lightIllu_procedure',timestamp,'.mat'];
    save(individual_mat_file.caliblightillu,'-v7.3');
        clear c timestamp;

    clear ref_board_calib_light_illu whiteboard whiteboard_plane_eqn 

    saveMatFile;
    
end    


function [whiteboard_plane_eqn, error, ref_board_calib_light_illu] =...
    getWhiteBoardPlane(ref_board_file_name, CamCalib, param, error)
    
    ref_board_calib_light_illu =...
        getBoardFromIm(ref_board_file_name,...
        param.calibcam.checkerboard_size, param.calibcam.checkerboard_grid_size,...
        param.calibcam.sigma, param.calibcam.min_l_th, param.calibcam.max_k_th,...
        param.calibcam.var_k_th, param.calibcam.var_l_th,...
        1, 1, param.calibcam.over_s, 0,...
        1);
    [ref_board_calib_light_illu, error_record] =...
        getBoardExtParam(ref_board_calib_light_illu, CamCalib, 5, 20);
    ref_board_calib_light_illu.plane_eqn =...
        linearFit(ref_board_calib_light_illu.c_X);

    thickness_whiteboard = param.caliblightillu.thickness_whiteboard;

    ref_plane_eqn = [-ref_board_calib_light_illu.plane_eqn;1] /...
        ref_board_calib_light_illu.plane_eqn(3) ;

    whiteboard_plane_eqn = [ref_plane_eqn(1:3);...
        ref_plane_eqn(4)-thickness_whiteboard/norm(ref_plane_eqn(1:3))];

    error.caliblightillu.error_whiteplane_x_proj = error_record(:,1:2);
end

function xyzC = getWhiteBoardxyzC(plane_eqn,CamCalib,region)
    x = region(1):region(3);
    y = region(2):region(4);
    NUM_Y = numel(y);
    NUM_X = numel(x);
    [xx,yy] = meshgrid(x,y);
    [xn] = normalize_pixel([xx(:),yy(:)]',CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
    Z = -plane_eqn(4)./(plane_eqn(1)*xn(1,:)+plane_eqn(2)*xn(2,:)-1);
    xyzC(:,:,1) = reshape(Z.*xn(1,:),NUM_Y,NUM_X);
    xyzC(:,:,2) = reshape(Z.*xn(2,:),NUM_Y,NUM_X);
    xyzC(:,:,3) = reshape(Z,NUM_Y,NUM_X);
end

function displayResult(figure_count, wb_im_ind, wb_im_trimed, ...
    wb_Intensity_filtered_ind, lightStrength_ind, ...
    cosThetai_whiteboard_ind, rerendering_error_gray, tire_region)
    figure(figure_count);
    for color_i = 1:3
        subplot(3,4,color_i);
        imagesc(wb_im_trimed(:,:,color_i));
        axis equal; colorbar;   colormap colorcube;
        title(['raw',num2str(color_i)]);

        subplot(3,4,color_i+4);
        imagesc(wb_Intensity_filtered_ind(:,:,color_i));
        axis equal; colorbar;   colormap colorcube;
        title(['filtered',num2str(color_i)]);

        subplot(3,4,color_i+4*2);
        imagesc(lightStrength_ind(:,:,color_i));
        axis equal; colorbar;   colormap colorcube;
        title(['lightStrength',num2str(color_i)]);
    end

    subplot(3,4,4);
    imagesc(cosThetai_whiteboard_ind);
    axis equal; colorbar;   colormap colorcube;
    title('cosThetai due to light position');

    max_lightStrength=max(wb_Intensity_filtered_ind(:));
    p=[tire_region(1),tire_region(2);...
        tire_region(1),tire_region(4);...
        tire_region(3),tire_region(4);...
        tire_region(3),tire_region(2);...
        tire_region(1),tire_region(2)];
    subplot(3,4,8);
    imshow(wb_im_ind/max_lightStrength);
    hold on; plot(p(:,1),p(:,2),'r');
    title(['rawIm,maxLS=',num2str(max_lightStrength)]);

    subplot(3,4,12);
    imagesc(rerendering_error_gray);
    axis equal; colorbar;   colormap colorcube;
    title(['Rerending error, meanError=',num2str(mean(rerendering_error_gray(:)))]);
    drawnow();
end