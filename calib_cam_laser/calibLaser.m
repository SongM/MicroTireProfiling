function calibLaser(mat_file)
    load(mat_file);
    load(individual_mat_file.camcalib, 'bs_cam_calib');

    redness_th = param.caliblaser.redness_th;
    line_err_mse_th = param.caliblaser.line_err_mse_th;
    
    
    % 2.B.1.	Test the first checkerboards contain laser points, choose redness_th, line_err_mse_th
    if(0)
        board_ind_need_to_be_checked = 1;
        PreProcessForLaserCalib(bs_cam_calib,...
            1,redness_th,line_err_mse_th, board_ind_need_to_be_checked);
        clear board_ind_need_to_be_checked;
    end
    
    %%
    % 2.B.2.	Detect laser points in checkerboards detected in 2.a.ii.
    task = '2.B.2.DetectLaserPointOnCheckerBoards';
    displayTask(task,3);
    
    [bs_laser, error.caliblaser.line_err] = PreProcessForLaserCalib(bs_cam_calib,...
        1,redness_th,line_err_mse_th);
    clear redness_th line_err_mse_th
    
    job_done = [job_done;task];
    clear task;
    clear bs_cam_calib_prepare bs_cam_calib desactivated_images_calib

%     saveMatFile;
    
    %%
    % 2.B.3.	Calib laser plane
    task = '2.B.3.CalibLaserPlane';
    displayTask(task,3);

    LaserCalib = calibLaserFromLaserBS(CamCalib,bs_laser,1);
    error.caliblaser.plane_err = LaserCalib.plane_err;
    job_done = [job_done;task];
    clear task;

    c = clock;
    disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
    timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];

    individual_mat_file.lasercalib = [data_folder, '\', case_name, '\processed_data\', 'lasercalib_procedure',timestamp,'.mat'];
    save(individual_mat_file.lasercalib);
    clear c timestamp;

    clear bs_laser

    saveMatFile;
    




end