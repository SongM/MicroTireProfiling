function manuallyChooseBallArea(mat_file)
    load(mat_file);
    
    %%
    % 2.C.1.	Group the pictures.    
    % 2.C.1.i.	The first picture is a picture of checkerboard, which gives the parameters of the checkerboard plane.
    % 2.C.1.ii.	The rest pictures are in a group of 9 pictures. The first picture of the group of pictures is used for find the center of the ball. The other eight are used for detect 

    ball_range_scale = param.caliblightpos.ball_range_scale;
    N_lights = param.N_lights;
    
    task = '2.C.1.GroupPic';
    displayTask(task,3);
    
    [ref_board_balls, balls] = groupPic(folder_path,param.file_name.caliblightpos,param.N_lights);
    job_done = [job_done;task];
    clear task;

    
    %%
    % 2.C.2.	Get the board plane eqn.
    task = '2.C.2.CalibBallPlane';
    displayTask(task,3);
    
    [balls_center_plane_eqn, error, ref_board_balls] =...
        getBallCenterPlane(ref_board_balls, CamCalib,...
        param, error, figure_count);
    
    balls_center_plane_eqn_3_variable = -balls_center_plane_eqn(1:3)/balls_center_plane_eqn(4);

    figure_count = figure_count+1;

    job_done = [job_done;task];
    clear task;
    
    
%%    
    % 2.C.3.	Find the center of the ball (C_k), and the most bright spot (s_i,k).
    task = '2.C.3.Find__C_k__and__s_ik';
    displayTask(task,3);
    
    figure_count = figure_count + 1;

    prompt={'b_ck_x','b_ck_y','b_radius'};
    dlg_title = 'Input';            
    radius = 100;
    center = [0,0];

    color_i = 1;
    
    N_balls = numel(balls);
    for ball_k = 1:N_balls
        light_i = 1;
        %% 2.C.3.i.	Find c_k
        t_im = getImFromPath(balls(ball_k).bright_spot_im_filepath{light_i});
        t_im = t_im/quantile(t_im(:), 1-(1e-2) );
        t_im(t_im>1) = 1;

        figure(figure_count);clf;set(figure_count,'OuterPosition',[10,200,800,600]);
        subplot(1,2,1);imshow(t_im);
        viscircles(center, radius, 'Color', 'r' ,'LineWidth', 1);
        axis on;

        answer = 0;

        while(numel(answer)>0)
            defaultans = {num2str(center(1)), num2str(center(2)), num2str(radius)};
            answer = inputdlg(prompt,dlg_title,1,defaultans);
            if (numel(answer)>0)
                center = [str2num(answer{1}),str2num(answer{2})];
                radius = str2num(answer{3});
                ball_area = [floor(center-radius*ball_range_scale) , repmat(ceil(radius*ball_range_scale*2),1,2 )];
                center_trimmed = center - ball_area(1:2) + 1;

                figure(figure_count);clf;set(figure_count,'OuterPosition',[10,200,800,600]);
                subplot(1,2,1);imshow(t_im);            
                viscircles(center, radius, 'Color', 'r' ,'LineWidth', 1);
                axis on;
                subplot(1,2,2);imshow(t_im( ball_area(2):(ball_area(2)+ball_area(4)), ...
                    ball_area(1):(ball_area(1)+ball_area(3)),...
                    color_i ));
                axis on;
                viscircles(center_trimmed, radius, 'Color', 'r' ,'LineWidth', 1);
            end
        end
        
        balls(ball_k).ball_center_plane_eqn = balls_center_plane_eqn_3_variable;
        balls(ball_k).b_ck = center;
        balls(ball_k).b_radius = radius;
        balls(ball_k).b_ball_area = ball_area;
        %  2.C.3.ii.	Find s_ik 
        b_sik = ones(N_lights,1)*center;
        balls(ball_k).b_sik = b_sik;

    end
    
    for ball_k = 1:N_balls
        fprintf(['working on ball_', num2str(ball_k),'/', num2str(N_balls), '...']);

        ball_area = balls(ball_k).b_ball_area;
        radius = balls(ball_k).b_radius;
        
        ball_im.bright_spot_size = ones(1,N_lights)*radius/2;
        
        ball_im.bright_spot_im_whole = zeros(ball_area(3)+1, ball_area(4)+1, N_lights);  
        
        for light_i = 1:N_lights
            fprintf(['.',num2str(light_i)]);
            clear t_im;

            t_im = getImFromPath(balls(ball_k).bright_spot_im_filepath{light_i});
            t_im = t_im/quantile(t_im(:), 1-(1e-2) );
            t_im(t_im>1) = 1;
            ball_im.bright_spot_im_whole(:,:,light_i) = t_im( ball_area(2):(ball_area(2)+ball_area(4)), ...
                    ball_area(1):(ball_area(1)+ball_area(3)),...
                    color_i );
        end
        fprintf('\n');
        balls_im{ball_k} = ball_im;
    end
    
    clear ans center center_trimmed radius b_sik ball_area ball_range_scale
    clear prompt defaultans dlg_title answer
    clear ball_k N_balls light_i N_lights t_im color_i ball_im
    
    job_done = [job_done;task];
    clear task;
    c = clock;
    disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
    timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];

    individual_mat_file.caliblightpos = [data_folder, '\', case_name, '\processed_data\', 'lightPos_procedure',timestamp,'.mat'];
    save(individual_mat_file.caliblightpos,'-v7.3');
    clear c timestamp;

    clear balls balls_im
    clear ref_board_balls balls_center_plane_eqn balls_center_plane_eqn_3_variable
%     clear ref_board_balls balls_center_plane_eqn ref_board_balls
%     save(lightPos_mat_file,'balls_im','balls','-v7.3');
%     clear 

    
%     clear bs_cam_calib_prepare bs_cam_calib desactivated_images_calib

    saveMatFile;
    
end

function [board, balls] = groupPic(folder_path,light_position_file_name,N_lights)
    folder_name = folder_path;
    im_file_list = dir([folder_name,'\',light_position_file_name,'\*JPG']);
%     im_file_list = dir([folder_name,'\',light_position_file_name,'\*DNG']);
    
    board.filepath = [folder_name,'\',light_position_file_name,'\',im_file_list(1).name];
    N_ball = floor(numel(im_file_list)/(N_lights+1));
    balls = [];

    for ball_i = 1:N_ball
        ball = [];
        ball.ind = ball_i;
        ball.center_im_filepath = [folder_name,'\',light_position_file_name,'\',im_file_list( (ball_i-1)*(N_lights+1) + 2 ).name];
        for light_i=1:N_lights
            ball.bright_spot_im_filepath{light_i,1} = [folder_name,'\',light_position_file_name,'\',im_file_list( (ball_i-1)*(N_lights+1) + 2 + light_i ).name];
        end
        balls=[balls,ball];
    end
end

function [balls_center_plane_eqn, error, ref_board_balls] =...
        getBallCenterPlane(ref_board_balls, CamCalib,...
        param, error, figure_idx)
    
    ref_board_balls = getBoardFromIm(ref_board_balls.filepath,...
        param.calibcam.checkerboard_size, param.calibcam.checkerboard_grid_size,...
        param.calibcam.sigma, param.calibcam.min_l_th, param.calibcam.max_k_th,...
        param.calibcam.var_k_th, param.calibcam.var_l_th,...
        1, 1, param.calibcam.over_s, 0,...
        1);
    
    [ref_board_balls, error_record] = getBoardExtParam(ref_board_balls, CamCalib, 5, 20);
        
%     displayCamCalib(CamCalib,ref_board_balls,1,figure_idx);
    
    d_ball_center_ref_plane = param.caliblightpos.d_ball_center_ref_plane; 
    
    ref_board_balls.plane_eqn = linearFit(ref_board_balls.c_X);
    ref_plane_eqn = [-ref_board_balls.plane_eqn;1]/ref_board_balls.plane_eqn(3);
    balls_center_plane_eqn = [ref_plane_eqn(1:3);...
        ref_plane_eqn(4)-d_ball_center_ref_plane/norm(ref_plane_eqn(1:3))];

    error.caliblightpos.error_ball_plane_x_proj = error_record(:,1:2);
end
