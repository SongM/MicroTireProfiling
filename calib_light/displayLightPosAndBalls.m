function figure_count = displayLightPosAndBalls( LightPos, balls, ref_board_balls, CamCalib, figure_count, light_idx_display, b_show_camera, b_show_ref_board)
%DISPLAYLIGHTPOSANDBALLS Summary of this function goes here
%   Detailed explanation goes here
    
    if (~exist('figure_count','var')),  figure_count = 1;  end
    if (~exist('b_show_camera','var')),  b_show_camera = 0;  end
    if (~exist('b_show_ref_board','var')),  b_show_ref_board = 1;  end
    if (~exist('light_idx_display','var')),  light_idx_display = 1;  end
    
    figure_count = figure_count + 1;
    cam_board_idx = figure_count; hold on;
    if (b_show_ref_board)
        displayCamCalib(CamCalib,ref_board_balls,b_show_camera,cam_board_idx);
    else
        displayCamCalib(CamCalib,[],b_show_camera,cam_board_idx);
    end        

    %% display light pos
    hold on;
    scatter3(LightPos(1,:), LightPos(3,:), -LightPos(2,:) , 10, 'r', 'filled');
    N_lights = size(LightPos,2);
    for light_i=1:N_lights
        text(LightPos(1,light_i), LightPos(3,light_i), -LightPos(2,light_i), num2str(light_i));
    end
    %% display balls
    if(b_show_ref_board)
    ball_N = numel(balls);
    for ball_k = 1:ball_N
        ball = balls(ball_k);
        C_ck = ball.C_ck;
        C_sik = ball.C_sik;
        
        
        

        % viewing vector in {C} (reflectance light)
        C_v_ik = 0 - C_sik; % camera center to the brightest spot
        C_v_ik_unit = normalizeColVector(C_v_ik);
        % halfway vector in {C} surface normal at the brightest spot
        C_h_ik = C_sik - repmat(C_ck,1,N_lights); % brightest point to the ball center
        C_h_ik_unit = normalizeColVector(C_h_ik);
        % light direction in {C} (incidence light)
        C_l_ik = 2*C_h_ik_unit - C_v_ik_unit; % light to the brightest spot;
        C_l_ik_unit = normalizeColVector(C_l_ik);
        
        
        
        RadiusFromImage = ball.C_Radius_FromImage;
        
        drawSphere([C_ck(1),C_ck(3),-C_ck(2)],RadiusFromImage);        
        text(C_ck(1)+10,C_ck(3),-C_ck(2),num2str(ball_k));

%         for light_i = 1:N_lights
        for light_i=light_idx_display
            C_sik_ind = C_sik(:,light_i);
            C_v_ik_unit_ind = C_v_ik_unit(:,light_i);
            C_l_ik_unit_ind = C_l_ik_unit(:,light_i);
            Z_l = 300;
            scale_v = (0-C_sik_ind(3)) / C_v_ik_unit_ind(3);
            scale_l = (Z_l-C_sik_ind(3)) / C_l_ik_unit_ind(3);
            
            line_l = [C_sik_ind+scale_l*C_l_ik_unit_ind, C_sik_ind];
            line_v = [C_sik_ind+scale_v*C_v_ik_unit_ind, C_sik_ind];
            
            
            hold on;
            plot3(line_l(1,:),line_l(3,:),-line_l(2,:),'r');
            plot3(line_v(1,:),line_v(3,:),-line_v(2,:),'g');
            drawnow();
        end

    end
end

