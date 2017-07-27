function LaserCalib = calibLaserFromLaserBS(CamCalib,bs_laser,b_display)

%     addpath('functions');
%     addpath('preprocessed_data');
%     folder_path = ['../',folder_name,'/'];
%     addpath(folder_path);


    if (b_display)
        figure();hold on;
        displayCamera(CamCalib);
    end
    
    color_map = hsv(numel(bs_laser)+1);
    X_c_laser_all = [];
    board_proc = 1:numel(bs_laser);
    for kk = board_proc
        board = bs_laser(kk);
        % figure();displayBoard2D(board);
        %   displayBoard3D(board,color_map(kk,:));
        X_c_laser = get3DCoorFromPlaneEqnAnd2DCoor(board.camera_plane.plane_eqn, board.laserCalib.p_laser_points, CamCalib);
        if (b_display)
            scatter3(X_c_laser(1,:),X_c_laser(2,:),X_c_laser(3,:),2,color_map(kk,:));
        end
        X_c_laser_all = [X_c_laser_all, X_c_laser];
        LaserCalib.c_calib_points{kk} = X_c_laser;
    end
    % fitting laser plane
%     LaserCalib.c_calib_points = X_c_laser_all;
    [LaserCalib.c_plane_eqn, LaserCalib.plane_err] = linearFit(X_c_laser_all);
    fprintf(1,['laser_plane_std = ', num2str(LaserCalib.plane_err.std),'.\n']);
    
    % for display
    laser_plane_y_min = min(X_c_laser_all(2,:));
    laser_plane_y_max = max(X_c_laser_all(2,:));
    laser_plane_z_min = 300;
    laser_plane_z_max = 1000;

    laser_plane_y = [laser_plane_y_min, laser_plane_y_min, laser_plane_y_max, laser_plane_y_max, laser_plane_y_min];
    laser_plane_z = [laser_plane_z_min, laser_plane_z_max, laser_plane_z_max, laser_plane_z_min, laser_plane_z_min];
    laser_plane_x = (1 - LaserCalib.c_plane_eqn(2)*laser_plane_y - LaserCalib.c_plane_eqn(3)*laser_plane_z)/LaserCalib.c_plane_eqn(1);

    LaserCalib.plane_display = [laser_plane_x;laser_plane_y;laser_plane_z];
    if (b_display)
        plot3(laser_plane_x,laser_plane_y,laser_plane_z,'Color',color_map(end,:));
    end
end