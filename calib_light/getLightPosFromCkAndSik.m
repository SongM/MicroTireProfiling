function [balls, LightPos, error] = getLightPosFromCkAndSik(balls, error)
    N_balls = numel(balls);
    N_lights = numel(balls(1).bright_spot_im_filepath);
    M = cell(N_lights,1);
    for ball_k = 1:N_balls
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

        
        balls(ball_k).C_v_ik_unit = C_v_ik_unit;
        balls(ball_k).C_h_ik_unit = C_h_ik_unit;
        balls(ball_k).C_l_ik_unit = C_l_ik_unit;
        
        % prepare for the equation to calculate LightPos
        for light_i = 1:N_lights
            l_bar = C_l_ik_unit(:,light_i); % light direction
            C_sik_ind = C_sik(:,light_i);   % brightest spot

            A = eye(3)-l_bar*l_bar';
            b = A*C_sik_ind;
            M{light_i} = [ M{light_i};[A,b] ];
        end
    end
    
    % Calculate LightPos
%     LightPos = [];
%     for light_i = 1:N_lights
%         A = M{light_i}(:,1:3);
%         b = M{light_i}(:,4);
%         LightPos(:,light_i) =A\b;
%     end
% 
%     % get the error
%     d_ik = [];
%     for ball_k = 1:N_balls
%         ball = balls(ball_k);
%         C_sik = ball.C_sik;
%         C_l_ik_unit = ball.C_l_ik_unit;
% 
%         for light_i = 1:N_lights
%             LightPos_ind = LightPos(:,light_i); % light position 
%             C_sik_ind = C_sik(:,light_i);   % light direction
%             l_bar = C_l_ik_unit(:,light_i); % brightest spot
% 
%             A = eye(3)-l_bar*l_bar';
%             d_ik(ball_k, light_i, :) = norm(A*(LightPos_ind-C_sik_ind));
%         end
%     end


    d_ik = [];
    LightPos = [];
    large_error_tolerance = 1.8;
    for light_i = 1:N_lights
        b_no_large_error = false;
        active_ball_ind = ones(1,N_balls);
        A = M{light_i}(:,1:3);
        b = M{light_i}(:,4);

        while(~b_no_large_error)
            active_param_ind = find(repmat(active_ball_ind,1,3)==1);
            A_active = A(active_param_ind,:);
            b_active = b(active_param_ind,:);
            LightPos_ind = A_active\b_active;

            d_ik_ind = [];
            for ball_k = 1:N_balls;
                if(active_ball_ind(ball_k)==1)
                    ball = balls(ball_k);
                    C_sik = ball.C_sik; % brightest spot - all
                    C_l_ik_unit = ball.C_l_ik_unit; % light direction - all

                    C_sik_ind = C_sik(:,light_i);   % brightest spot
                    l_bar = C_l_ik_unit(:,light_i); % light direction

                    AA = eye(3)-l_bar*l_bar';
                    d_ik_ind(ball_k,:) = norm(AA*(LightPos_ind-C_sik_ind));
                else
                    d_ik_ind(ball_k,:) = -1;
                end
            end
            large_error_ind = find(d_ik_ind>large_error_tolerance*mean(d_ik_ind(d_ik_ind>0)));
            if(isempty(large_error_ind))
                b_no_large_error = true;
            else
                active_ball_ind(large_error_ind) = 0;
            end
        end
%         d_ik_ind
        d_ik(:,light_i) = d_ik_ind;
        LightPos(:,light_i) = LightPos_ind;

    end
    error.caliblightpos.d_ik = d_ik;
    
    
    error_fig = figure;
    set(error_fig,'Name','light position error');
    d_ik(d_ik<0)=-max(d_ik(:));
    subplot(1,2,1);imagesc(d_ik);axis equal;colorbar;
    xlabel('light_i');
    ylabel('ball_k');
    title({'minus value indicate larger error', 'positive value indicates how far the light', 'position is away from the light direction'})
    
    
    temp_ind = repmat(1:N_lights,N_balls,1);
    subplot(1,2,2);
    boxplot(d_ik(d_ik>0),temp_ind(d_ik>0));
    title('boxplot of distances from each light to its own light positions');



%     ball_N = numel(balls);
%     N_lights = numel(balls(1).bright_spot_im_filepath);
%     M = cell(N_lights,1);
%     for ball_k = 1:ball_N
%         ball = balls(ball_k);
%         C_ck = ball.C_ck;
%         C_sik = ball.C_sik;
%         for light_i = 1:N_lights
%             C_sik_ind = C_sik(:,light_i);
%             s_bar = (C_sik_ind - C_ck)/norm(C_sik_ind - C_ck);
%             A = eye(3)-s_bar*s_bar';
%             b = A*C_ck;
%             M{light_i} = [ M{light_i};[A,b] ];
%         end
%     end
%     
%     LightPos = [];
%     for light_i = 1:N_lights
%         A = M{light_i}(:,1:3);
%         b = M{light_i}(:,4);
%         LightPos(:,light_i) =A\b;
%     end
% 
%     d_ik = [];
%     for ball_k = 1:ball_N
%         ball = balls(ball_k);
%         C_ck = ball.C_ck;
%         C_sik = ball.C_sik;
%         for light_i = 1:N_lights
%             C_li = LightPos(:,light_i);
%             C_sik_ind = C_sik(:,light_i);
%             s_bar = (C_sik_ind - C_ck)/norm(C_sik_ind - C_ck);
%             A = eye(3)-s_bar*s_bar';
%             d_ik(ball_k, light_i, :) = norm(A*(C_li-C_ck));
%         end
%     end
%     error.caliblightpos.d_ik = d_ik;
end
