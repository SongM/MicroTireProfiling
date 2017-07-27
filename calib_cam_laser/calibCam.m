function calibCam(mat_file)
    load(mat_file);
    folder_name_calib = param.calibcam.folder_name_calib;
    checkerboard_size = param.calibcam.checkerboard_size;
    checkerboard_grid_size = param.calibcam.checkerboard_grid_size;
    sigma = param.calibcam.sigma;
    min_l_th = param.calibcam.min_l_th;
    max_k_th = param.calibcam.max_k_th;
    var_k_th = param.calibcam.var_k_th;
    var_l_th = param.calibcam.var_l_th;
    over_s = param.calibcam.over_s;
    %%
    if(0)
        % 2.	Calibrate camera intrinsic and laser plane
        % 2.A.	Camera intrinsic calib
        % 2.A.1.	Test the first image of cam&laser calib, choose checkerboard_size, checkerboard_grid_size, sigma, min_l_th, var_k_th, var_l_th.
        im_ind_need_to_be_checked = 1;
        getBoardFromFolder(folder_name_calib,...
            checkerboard_size, checkerboard_grid_size,...
            sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
            1, 1, over_s, 0,...
            1, figure_count, im_ind_need_to_be_checked);
    end
    
    %%    
    % 2.A.2.	Detect checkerboards in the images of cam&laser calib.
    load(mat_file);
    task = '2.A.2.DetectCheckerBoardForCamCalib';
    displayTask(task,3);
    bs_cam_calib_prepare = getBoardFromFolder(folder_name_calib,...
        checkerboard_size, checkerboard_grid_size,...
        sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
        1, 1, over_s, 0,...
        1);

    clear folder_name_calib;
    clear checkerboard_size checkerboard_grid_size;
    clear sigma min_l_th max_k_th var_k_th var_l_th over_s;

    job_done = [job_done;task];
    clear task;
%     saveMatFile;
    
    %%    
    % 2.A.3.	Calibrate camera intrinsic using the checkerboards detected in 2.a.ii
%     load(mat_file);
    task = '2.A.3.getCamCalibIntrinsics';
    displayTask(task,3);

    winSize = param.calibcam.winSize;
    
    [CamCalib,bs_cam_calib,desactivated_images_calib] =...
    getCalibCamParams(bs_cam_calib_prepare,winSize);
    clear winSize;
    
    % 2.A.3.i.	Error CamCalib, the reprojection error.
    error.calibcam.cc_error = CamCalib.cc_error;
    error.calibcam.fc_error = CamCalib.fc_error;
    error.calibcam.kc_error = CamCalib.kc_error;
    error.calibcam.x_proj_err_std = CamCalib.err_std;
    for board_i = 1:numel(bs_cam_calib)
        error.calibcam.x_proj_err{board_i} = bs_cam_calib(board_i).error.p_x_proj_error;
    end
    clear board_i;

    job_done = [job_done;task];
    clear task;
    c = clock;
    disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
    timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];

    individual_mat_file.camcalib = [data_folder, '\', case_name, '\processed_data\', 'camcalib_procedure',timestamp,'.mat'];
    save(individual_mat_file.camcalib,'-v7.3');
    clear c timestamp;
    clear bs_cam_calib_prepare bs_cam_calib desactivated_images_calib
    saveMatFile;
end
