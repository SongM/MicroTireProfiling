function [balls, figure_count, balls_im] = Find__C_k__and__s_ik...
    (param, CamCalib, balls, balls_center_plane_eqn, blackboard_region, ref_board_balls,...
    b_display, figure_count)

    tic;
    balls_center_plane_eqn_3_variable = -balls_center_plane_eqn(1:3)/balls_center_plane_eqn(4);
    ball_size_range = param.caliblightpos.ball_size_range;
    resize_scale = param.caliblightpos.resize_scale;
    ball_range_scale = param.caliblightpos.ball_range_scale;
    brightness_th = param.caliblightpos.brightness_th;
    ratio_radius_brightness_spot = param.caliblightpos.ratio_radius_brightness_spot;
    
    N_lights = param.N_lights;

    if (b_display)
        figure_count = figure_count + 1;
        cam_board_idx = figure_count;
        displayCamCalib(CamCalib,ref_board_balls,1,cam_board_idx);
        drawnow();
    end

    N_balls = numel(balls);
    for ball_k = 1:N_balls
        fprintf(['working on ball_', num2str(ball_k),'/', num2str(N_balls), '...']);
        temp_im = 0;
        ims = [];
      
        fprintf([':',num2str(N_lights)]);
        ball = balls(ball_k);
        for im_i = 1:N_lights
            fprintf(['.',num2str(im_i)]);
            t_im = getImFromPath( ball.bright_spot_im_filepath{im_i} );
            t_im = t_im/quantile(t_im(:), 1-(1e-2) );
            t_im(t_im>1) = 1;
            ims{im_i} = t_im;
            clear t_im;
            temp_im = temp_im + ims{im_i};
        end
        fprintf('\n');
        t_temp_im = temp_im/N_lights;
        temp_im = t_temp_im*0;
%         temp_im(blackboard_region(1):blackboard_region(3),blackboard_region(2):blackboard_region(4),:) = ...
%             t_temp_im(blackboard_region(1):blackboard_region(3),blackboard_region(2):blackboard_region(4),:);

        temp_im(blackboard_region(2):blackboard_region(4),blackboard_region(1):blackboard_region(3),:) = ...
            t_temp_im(blackboard_region(2):blackboard_region(4),blackboard_region(1):blackboard_region(3),:);

        
        
        %% 2.C.3.i.	Find c_k
        [b_ck, ball_area, center_trimmed, radius, figure_count] =...
            Find__C_k(temp_im, resize_scale, ball_size_range, ball_range_scale,...
            b_display, figure_count);
        %%  2.C.3.ii.	Find s_ik 
        [b_sik, figure_count, balls_im{ball_k}] = Find__S_ik(ims, temp_im,...
            ball_area, center_trimmed, b_ck, radius,...
            ratio_radius_brightness_spot, brightness_th,...
            b_display, figure_count);

        C_ck = get3DCoorFromPlaneEqnAnd2DCoor...
            (balls_center_plane_eqn_3_variable, b_ck', CamCalib);
        C_sik = get3DCoorOfBrightestSpot...
            (b_sik', C_ck, param.caliblightpos.ball_radius, CamCalib);

        RadiusFromImage = radius/mean(CamCalib.fc)*C_ck(3);
        fprintf(['b_ck = [', num2str(b_ck), '], C_ck = [', num2str(C_ck'),...
            '], radius = ', num2str(RadiusFromImage),...
            ' (default_radius = ',num2str(param.caliblightpos.ball_radius),')',...
            ', time_remain = ', num2str( (N_balls - ball_k)*(toc)/ball_k ),'.\n']);


        if(b_display)
            figure(cam_board_idx);
            drawSphere([C_ck(1),C_ck(3),-C_ck(2)],RadiusFromImage);        

            Z_min = 300;

            len = (Z_min-C_ck(3))./(C_sik(3,:)-C_ck(3));
            line_endpoint_at_Z_min = (C_sik - repmat(C_ck,1,numel(ims))).*repmat(len,3,1) + repmat(C_ck,1,numel(ims)    );

            for light_i = 1:numel(ims)
                line_light = [line_endpoint_at_Z_min(:,light_i),C_ck];
                hold on;
                plot3(line_light(1,:),line_light(3,:),-line_light(2,:),'r');
            end
            drawnow();
        end


        balls(ball_k).ball_center_plane_eqn = balls_center_plane_eqn_3_variable;
        balls(ball_k).C_Radius_measured = param.caliblightpos.ball_radius;

        balls(ball_k).b_ck = b_ck;
        balls(ball_k).C_ck = C_ck;

        balls(ball_k).b_radius = radius;
        balls(ball_k).C_Radius_FromImage = RadiusFromImage;

        balls(ball_k).b_ball_area = ball_area;

        balls(ball_k).b_sik = b_sik;
        balls(ball_k).C_sik = C_sik;
    end
    
    
    
    
    
    
end

function [b_sik, figure_count, ball_im] = Find__S_ik(ims, temp_im,...
    ball_area, center_trimmed, b_ck, radius,...
    ratio_radius_brightness_spot, brightness_th,...
    b_display, figure_count)
    % to find the brightest spot, we only convern the bright pixels on the
    % ball surface, we use a circle mask to enable that only the bright
    % pixels on the ball surface will be count.
    
    circle_mask = zeros(ball_area(3)+1,ball_area(4)+1);
    [ny,nx] = size(circle_mask);
    [X,Y] = meshgrid(1:nx,1:ny);
    circle_ind = ( (X-center_trimmed(1)).^2+(Y-center_trimmed(2)).^2 ) <...
        (ratio_radius_brightness_spot*radius)^2;
    circle_mask(circle_ind) = 1;
%     circle_mask = repmat(circle_mask,1,1,3);
    
    
    if(b_display)
        figure_count = figure_count+1;
        figure(figure_count); 
        set(figure_count,'Name','b_sik');
    end
    
    b_sik = [];
%     ball_im.temp_im = temp_im;
    ball_im.brgiht_spot_im = zeros(ball_area(3)+1,ball_area(4)+1,numel(ims));
    ball_im.brgiht_spot_im_whole = zeros(ball_area(3)+1,ball_area(4)+1,numel(ims));
    ball_im.bright_spot_size = zeros(1,numel(ims));
    for im_i = 1:numel(ims)
        
        im_trimmed_gray = ims{im_i}( ball_area(2):(ball_area(2)+ball_area(4)), ...
        ball_area(1):(ball_area(1)+ball_area(3)), 1 ).*circle_mask;

        ball_im.brgiht_spot_im(:,:,im_i) = im_trimmed_gray;
        ball_im.brgiht_spot_im_whole(:,:,im_i) = ims{im_i}( ball_area(2):(ball_area(2)+ball_area(4)), ...
        ball_area(1):(ball_area(1)+ball_area(3)), 1 );
        
%         im_trimmed_gray = rgb2gray(imcrop(ims{im_i},ball_area).*circle_mask);
        
        spot_brightness_th = mean( im_trimmed_gray(im_trimmed_gray>brightness_th) );
        [bright_spot_x_trimmed, bright_spot_y_trimmed] =...
            find( im_trimmed_gray>spot_brightness_th );
        center_spot_brightness =...
            [mean(bright_spot_y_trimmed), mean(bright_spot_x_trimmed)];
        sb_minx = min(bright_spot_x_trimmed);
        sb_miny = min(bright_spot_y_trimmed);
        sb_maxx = max(bright_spot_x_trimmed);
        sb_maxy = max(bright_spot_y_trimmed);
        
        b_sik(im_i,:) = center_spot_brightness + ball_area(1:2) - 1;
        
        ball_im.bright_spot_size(1,im_i) = (sb_maxx + sb_maxy - sb_minx - sb_miny)/2;

        if(b_display)
            subplot(4,4,im_i); imshow(im_trimmed_gray);
            hold on;         
            scatterPoints(b_ck, 3, 'r');
            scatterPoints(center_spot_brightness, 10, 'g');
            rectangle('Position', [sb_miny, sb_minx, sb_maxy-sb_miny, sb_maxx-sb_minx], 'EdgeColor', 'g');
        end
    end
    drawnow();
%     if(b_display)
%         figure(figure_count);
%         imshow(temp_im);viscircles(b_ck,radius);
%         hold on; scatterPoints(b_sik, 3, 'r');
%         scatterPoints(b_ck, 3, 'r');
%         figure_count = figure_count+1;
%         drawnow();
%     end
    

end

function C_sik = get3DCoorOfBrightestSpot(b_sik, C_ck, R, CamCalib)

    fc = CamCalib.fc;
    cc = CamCalib.cc;
    kc = CamCalib.kc;
    alpha_c = CamCalib.alpha_c;
    %     kc = zeros(5,1);

    i_x_undistort = normalize_pixel(b_sik, fc, cc, kc, alpha_c);
    l_s = [i_x_undistort;ones(1,size(i_x_undistort,2))];
    alpha = sum(l_s.^2);
    beta = C_ck'*l_s;
    gamma = sum(C_ck.^2)-R^2;
    C_sik_Z = (beta - sqrt(beta.^2-gamma*alpha))./alpha;
    C_sik = l_s.*repmat(C_sik_Z,3,1);
end


