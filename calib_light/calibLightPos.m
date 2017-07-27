function calibLightPos(mat_file)
    load(mat_file);
    
    %%
    % 2.C.1.	Group the pictures.    
    % 2.C.1.i.	The first picture is a picture of checkerboard, which gives the parameters of the checkerboard plane.
    % 2.C.1.ii.	The rest pictures are in a group of 9 pictures. The first picture of the group of pictures is used for find the center of the ball. The other eight are used for detect 

%     blackboard_region = [1286, 1697, 4100, 5842]
    blackboard_region = param.caliblightpos.blackboard_region;
    
    task = '2.C.1.GroupPic';
    displayTask(task,3);
    
    [ref_board_balls, balls] = groupPic(folder_path,param.file_name.caliblightpos,param.N_lights);
    job_done = [job_done;task];
    clear task;

    %% display parameters
    t_im = getImFromPath(balls(1).bright_spot_im_filepath{1});
    t_im = t_im/quantile(t_im(:), 1-(1e-2) );
    t_im(t_im>1) = 1;

    figure();imshow(t_im);
    hold on; rectangle('Position',[blackboard_region(1:2),blackboard_region(3:4)-blackboard_region(1:2)],'EdgeColor','r');
    viscircles((blackboard_region(1:2)+blackboard_region(3:4)) / 2, param.caliblightpos.ball_size_range(1));
    viscircles((blackboard_region(1:2)+blackboard_region(3:4)) / 2, param.caliblightpos.ball_size_range(2));
    drawnow();

    
    
    
    %%
    % 2.C.2.	Get the board plane eqn.
    task = '2.C.2.CalibBallPlane';
    displayTask(task,3);
    
    [balls_center_plane_eqn, error, ref_board_balls] =...
        getBallCenterPlane(ref_board_balls, CamCalib,...
        param, error, figure_count);
    figure_count = figure_count+1;

    job_done = [job_done;task];
    clear task;
    
    
%%    
    % 2.C.3.	Find the center of the ball (C_k), and the most bright spot (s_i,k).
    task = '2.C.3.Find__C_k__and__s_ik';
    displayTask(task,3);
    b_display = 1;

        [balls, figure_count, balls_im] = Find__C_k__and__s_ik...
        (param, CamCalib, balls, balls_center_plane_eqn, blackboard_region, ref_board_balls,...
        b_display, figure_count);
    

    
    
    job_done = [job_done;task];
    clear task b_display;
    %%
    % 2.C.4.	Calculate the light position L_i
    task = '2.C.4.getLightPosFromCkAndSik';
    displayTask(task,3);
    [balls, LightPos, error] = getLightPosFromCkAndSik(balls, error);
    fprintf(['LightPos = ']);
    disp(LightPos);
    fprintf(['d_ik = ']);
    disp(error.caliblightpos.d_ik);
    
    
    

    job_done = [job_done;task];
    clear task;
    c = clock;
    disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
    timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];

    individual_mat_file.caliblightpos = [data_folder, '\', case_name, '\processed_data\', 'lightPos_procedure',timestamp,'.mat'];
    save(individual_mat_file.caliblightpos,'-v7.3');
    clear c timestamp;

    clear ref_board_balls balls balls_im
    clear balls_center_plane_eqn ref_board_balls
%     save(lightPos_mat_file,'balls_im','balls','-v7.3');
%     clear 

    
%     clear bs_cam_calib_prepare bs_cam_calib desactivated_images_calib

    saveMatFile;
end

function [board, balls] = groupPic(folder_path,light_position_file_name,N_lights)
    folder_name = folder_path;
%     im_file_list = dir([folder_name,'\',light_position_file_name,'\*JPG']);
    im_file_list = dir([folder_name,'\',light_position_file_name,'\*DNG']);
    
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
