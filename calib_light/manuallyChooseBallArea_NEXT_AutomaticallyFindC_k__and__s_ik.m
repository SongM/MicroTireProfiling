function manuallyChooseBallArea_NEXT_AutomaticallyFindC_k__and__s_ik(mat_file)
    load(mat_file);
    load(individual_mat_file.caliblightpos,'balls','balls_im','ref_board_balls','balls_center_plane_eqn_3_variable');
    tic;
    N_balls = numel(balls);
    N_lights = param.N_lights;
    ball_size_range_margin = 0.3;
    brightness_th = 0.8;
    ratio_radius_brightness_spot = 0.8;
    b_display = 1;

    %% b_ck & b_sik
    
    if (b_display)
        figure_count = figure_count + 1;
        cam_board_idx = figure_count;
        displayCamCalib(CamCalib,ref_board_balls,1,cam_board_idx);
        drawnow();
    end

    
    for ball_k = 1:N_balls
        ball = balls(ball_k);
        ball_im = balls_im{ball_k};
        

        [b_ck, radius, b_sik, ball_im] = get_b_ck_and_b_sik(ball,ball_im,N_lights,...
            ball_size_range_margin, ratio_radius_brightness_spot, brightness_th, b_display);

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
            line_endpoint_at_Z_min = (C_sik - repmat(C_ck,1,N_lights)).*repmat(len,3,1) + repmat(C_ck,1,N_lights);

            for light_i = 1:N_lights
                line_light = [line_endpoint_at_Z_min(:,light_i),C_ck];
                hold on;
                plot3(line_light(1,:),line_light(3,:),-line_light(2,:),'r');
            end
            drawnow();
            clear Z_min len line_endpoint_at_Z_min light_i line_light
        end


        balls(ball_k).ball_center_plane_eqn = balls_center_plane_eqn_3_variable;
        balls(ball_k).C_Radius_measured = param.caliblightpos.ball_radius;

        balls(ball_k).b_ck = b_ck;
        balls(ball_k).C_ck = C_ck;

        balls(ball_k).b_radius = radius;
        balls(ball_k).C_Radius_FromImage = RadiusFromImage;

%         balls(ball_k).b_ball_area = ball_area;

        balls(ball_k).b_sik = b_sik;
        balls(ball_k).C_sik = C_sik;

        
        balls_im{ball_k} = ball_im;
    end
    clear ans ball ball_im ball_k N_balls N_lights b_display cam_board_idx
    clear ball_size_range_margin ball_size_range_margin brightness_th ratio_radius_brightness_spot
    clear b_ck radius b_sik C_ck C_sik RadiusFromImage

    save(individual_mat_file.caliblightpos,'-v7.3');
    
    clear balls balls_im
    clear ref_board_balls balls_center_plane_eqn balls_center_plane_eqn_3_variable
    saveMatFile;

end


function [b_ck, radius, b_sik, ball_im] = get_b_ck_and_b_sik(ball,ball_im,N_lights,...
    ball_size_range_margin, ratio_radius_brightness_spot, brightness_th, b_display)
    %% b_ck        
    temp_im = sum(ball_im.bright_spot_im_whole,3)/N_lights;
    ball_area = ball.b_ball_area;

    sensitivity = 0.9;
    radius = [];
    while (numel(radius)==0)
        [center_trimmed, radius] =...
            imfindcircles(temp_im, ball.b_radius*(1+ball_size_range_margin*[-1,1]), ...
            'ObjectPolarity','dark', 'Sensitivity',sensitivity);
        sensitivity = sensitivity + (1-sensitivity)*0.1;
        if (sensitivity>0.999)
            radius = ball.b_radius;
            center_trimmed = ball.b_ck - ball_area(1:2) + 1;
        end
    end

    b_ck = center_trimmed(1,:)+ball_area(1:2)-1;
    radius = radius(1,:);
    center_trimmed = center_trimmed(1,:);

    if(b_display)
        figure();imshow(temp_im);
        viscircles(center_trimmed(1,:), radius(1,:),'Color','r');
        figure();
    end


    %% b_sik
    circle_mask = zeros(ball_area(3)+1,ball_area(4)+1);
    [ny,nx] = size(circle_mask);
    [X,Y] = meshgrid(1:nx,1:ny);
    circle_ind = ( (X-center_trimmed(1)).^2+(Y-center_trimmed(2)).^2 ) <...
        (ratio_radius_brightness_spot*radius)^2;
    circle_mask(circle_ind) = 1;
    b_sik = ball.b_sik;

    for light_i = 1:N_lights
        im_trimmed_gray = ball_im.bright_spot_im_whole(:,:,light_i).*circle_mask;
        spot_brightness_th = mean( im_trimmed_gray(im_trimmed_gray>brightness_th) );
        [bright_spot_x_trimmed, bright_spot_y_trimmed] =...
            find( im_trimmed_gray>spot_brightness_th );

        if (numel(bright_spot_x_trimmed>0))
            center_spot_brightness =...
                [mean(bright_spot_y_trimmed), mean(bright_spot_x_trimmed)];
            sb_minx = min(bright_spot_x_trimmed);
            sb_miny = min(bright_spot_y_trimmed);
            sb_maxx = max(bright_spot_x_trimmed);
            sb_maxy = max(bright_spot_y_trimmed);

            b_sik(light_i,:) = center_spot_brightness + ball_area(1:2) - 1;

            ball_im.bright_spot_size(1,light_i) = (sb_maxx + sb_maxy - sb_minx - sb_miny)/2;
        

            if(b_display)

                subplot(4,4,light_i); imshow(im_trimmed_gray);
                hold on;         
                scatterPoints(b_ck, 3, 'r');
                scatterPoints(center_spot_brightness, 10, 'g');
                rectangle('Position', [sb_miny, sb_minx, sb_maxy-sb_miny, sb_maxx-sb_minx], 'EdgeColor', 'g');
            end
        end
    end
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