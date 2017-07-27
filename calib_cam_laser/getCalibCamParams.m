function [CamCalib,board_calib,desactivated_images_calib] = getCalibCamParams(calib_cam_board_set,wintx_default,winty_default)
    tic;
    center_optim = 0;
    if (0)
    est_dist = [0;0;0;0;0];
    end

    if (nargin<2),  wintx_default = 5;  end
    if (nargin<3),  winty_default = wintx_default;  end
%     folder_path = ['../',folder_name,'/'];
%     addpath('functions');
%     addpath('preprocessed_data');
%     addpath(folder_path);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    board_calib = calib_cam_board_set;
    %%%%%% load boards for calibration
    for board_i = 1:numel(board_calib)
        board = board_calib(board_i);
        ima_name = board.filepath;

        disp(['Loading board ', num2str(board_i), '/', num2str(numel(board_calib)),...
            ', BoardIdx = ', num2str(board.BoardIdx),...
            ', filepath = ', ima_name, '.']);
        
        
        Ii = double(getImFromPath(ima_name));
        Ii = 0.299 * Ii(:,:,1) + 0.5870 * Ii(:,:,2) + 0.114 * Ii(:,:,3);
        eval(['I_' num2str(board_i) ' = Ii;']);
    end

    fprintf(1,'\ndone\n');

    %%%%%% detect board corner
    [Hcal,Wcal] = size(I_1); 	% size of the calibration image
    [ny,nx] = size(I_1);
    clickname = [];
    map = gray(256);

    dX = board.checkerboard.dX;
    dY = board.checkerboard.dY;
    n_sq_x = board.checkerboard.nx-2;
    n_sq_y = board.checkerboard.ny-2;
    wintx = wintx_default;
    winty = winty_default;

    fprintf(1,'Using (wintx,winty)=(%d,%d) - Window size = %dx%d\n',wintx,winty,2*wintx+1,2*winty+1);
    fprintf(1,['Size of each square: [dX, dY] = [' , num2str(dX), ', ', num2str(dY), ']mm\n']);
    fprintf(1,'Processing image ');

    %%%%%% detect boards' corners
    board_proc = 1:numel(board_calib);
    for kk = board_proc
        fprintf(1,'%d...',kk);

        eval(['I = I_' num2str(kk) ';']);
        board = board_calib(kk);

        [x,X] = getBoardCorners(board, I, 1,...
            wintx, winty, n_sq_x, n_sq_y, dX, dY, map);
        eval(['dX_' num2str(kk) ' = dX;']);
        eval(['dY_' num2str(kk) ' = dY;']);  

        eval(['wintx_' num2str(kk) ' = wintx;']);
        eval(['winty_' num2str(kk) ' = winty;']);

        eval(['x_' num2str(kk) ' = x;']);
        eval(['X_' num2str(kk) ' = X;']);

        eval(['n_sq_x_' num2str(kk) ' = n_sq_x;']);
        eval(['n_sq_y_' num2str(kk) ' = n_sq_y;']);

    end
    fprintf(1,'\ndone\n');

    %%%%%% calibrate boards
    n_ima = numel(board_calib);
    active_images=ones(1,n_ima);
    % est_dist = zeros(5,1);
    go_calib_optim;


    %% check desactive
    if  ~isempty(desactivated_images)
        fprintf(2,'some boards in board_calib is not good\n');
    end    


    %%%%% save calibrated parameter to boards
    for kk = board_proc

        eval(['X_kk = X_', num2str(kk), ';']);
        eval(['Rc_kk = Rc_', num2str(kk), ';']);
        eval(['Tc_kk = Tc_', num2str(kk), ';']);
        Xc_kk = Rc_kk*X_kk + repmat(Tc_kk,1,size(X_kk,2));
        Xc_x = Xc_kk(1,:);
        Xc_y = Xc_kk(2,:);
        Xc_z = Xc_kk(3,:);

        Xp = Xc_x./Xc_z;
        Yp = Xc_y./Xc_z;
        Xp_r = sqrt(Xp.^2+Yp.^2);

        xx = Xp .* (1+kc(1)*Xp_r.^2+kc(2)*Xp_r.^4 + kc(5)*Xp_r.^6) + ...
            2*kc(3)*Xp.*Yp + kc(4)*(Xp_r.^2 + 2*Xp.^2);
        yy = Yp .* (1+kc(1)*Xp_r.^2+kc(2)*Xp_r.^4 + kc(5)*Xp_r.^6) + ...
            kc(3)*(Xp_r.^2 + 2*Yp.^2) + 2*kc(4)*Xp.*Yp;

        xxp = fc(1)*(xx+alpha_c*yy) + cc(1);
        yyp = fc(2)*yy + cc(2);
        y_kk = [xxp;yyp];

        eval(['omc_kk = omc_', num2str(kk), ';']);
        eval(['omc_error_kk = omc_error_', num2str(kk), ';']);
        eval(['Tc_error_kk = Tc_error_', num2str(kk), ';']);

        eval(['x_kk = x_', num2str(kk), ';']);

        ex_kk = y_kk - x_kk;
        
        % other info
        board_calib(kk).otherInfo.energy = board_calib(kk).energy;
        board_calib(kk).otherInfo.idx = board_calib(kk).idx;


        % coor in coordinate systems.
        board_calib(kk).p_x = x_kk;
        board_calib(kk).p_x_proj = y_kk;
        board_calib(kk).b_X = X_kk;
        board_calib(kk).c_X = Xc_kk;
        
        % board plane to camera plane
        board_calib(kk).b2c_omc = omc_kk;
        board_calib(kk).b2c_Rc = Rc_kk;
        board_calib(kk).b2c_Tc = Tc_kk;

        
        % pixel plane
        board_calib(kk).pixel_plane.v1 = board_calib(kk).v1;
        board_calib(kk).pixel_plane.v2 = board_calib(kk).v2;
        board_calib(kk).pixel_plane.coor_x = board_calib(kk).coor_x;
        board_calib(kk).pixel_plane.coor_y = board_calib(kk).coor_y;
        board_calib(kk).pixel_plane.corner = ...
            [board_calib(kk).cornerx;board_calib(kk).cornery];
        
        
        % camera plane
        n_grid_x = board_calib(kk).checkerboard.nx -1;
        n_grid_y = board_calib(kk).checkerboard.ny -1;
            % grid
        board_calib(kk).camera_plane.coor_X =...
            fliplr(reshape(Xc_kk(1,:),n_grid_x,n_grid_y));
        board_calib(kk).camera_plane.coor_Y =...
            fliplr(reshape(Xc_kk(2,:),n_grid_x,n_grid_y));
        board_calib(kk).camera_plane.coor_Z =...
            fliplr(reshape(Xc_kk(3,:),n_grid_x,n_grid_y));
        
            % Vector
        c_grid(:,:,1) = board_calib(kk).camera_plane.coor_X;
        c_grid(:,:,2) = board_calib(kk).camera_plane.coor_Y;
        c_grid(:,:,3) = board_calib(kk).camera_plane.coor_Z;
        
        grid_V1_1 = c_grid(:,1:end-1,:);
        grid_V1_2 = c_grid(:,2:end,:);
        V1_all = grid_V1_2 - grid_V1_1;
        c_V1 = reshape(mean(mean(V1_all)),1,3);
        board_calib(kk).camera_plane.V1 = c_V1;
        
        grid_V2_1 = c_grid(1:end-1,:,:);
        grid_V2_2 = c_grid(2:end,:,:);
        V2_all = grid_V2_2 - grid_V2_1;
        c_V2 = reshape(mean(mean(V2_all)),1,3);
        board_calib(kk).camera_plane.V2 = c_V2;
        
            % Corner
            if (board_calib(kk).BoardIdx == 0)
                n1 = 1;
            else
                n1 = 3;
            end
        c_corner(1,:) = reshape(c_grid(1,1,:),1,3) - c_V1 - n1*c_V2;
        c_corner(2,:) = reshape(c_grid(1,end,:),1,3) + c_V1 - n1*c_V2;
        c_corner(3,:) = reshape(c_grid(end,end,:),1,3) + c_V1 + c_V2;
        c_corner(4,:) = reshape(c_grid(end,1,:),1,3) - c_V1 + c_V2;
        
        board_calib(kk).camera_plane.Corner = c_corner';
       
        [board_calib(kk).camera_plane.plane_eqn,plane_err] = linearFit(Xc_kk);
%         [plane_err]
        if (plane_err.std>0.1)
            disp(['board ', num2str(kk),' plane_err.std > 0.1.']);
        end
        
        % error
        board_calib(kk).error.b2c_Tc_error = Tc_error_kk;
        board_calib(kk).error.b2c_omc_error = omc_error_kk;
        board_calib(kk).error.p_x_proj_error = ex_kk;
    end
    
    
    
    board_calib = rmfield(board_calib,'v1');
    board_calib = rmfield(board_calib,'v2');
    board_calib = rmfield(board_calib,'coor_x');
    board_calib = rmfield(board_calib,'coor_y');
    board_calib = rmfield(board_calib,'cornerx');
    board_calib = rmfield(board_calib,'cornery');
    board_calib = rmfield(board_calib,'energy');
    board_calib = rmfield(board_calib,'idx');
    
    % err_std_calib = err_std;
    desactivated_images_calib = desactivated_images;
    CamCalib.cc = cc;
    CamCalib.fc = fc;
    CamCalib.kc = kc;
    CamCalib.alpha_c = alpha_c;
    CamCalib.KK = KK;
    CamCalib.cc_error = cc_error;
    CamCalib.fc_error = fc_error;
    CamCalib.kc_error = kc_error;
    
    IP = 5*dX*[1 -alpha_c 0;0 1 0;0 0 1]*[1/fc(1) 0 0;0 1/fc(2) 0;0 0 1]*[1 0 -cc(1);0 1 -cc(2);0 0 1]*[0 nx-1 nx-1 0 0 ; 0 0 ny-1 ny-1 0;1 1 1 1 1];
    BASE = 5*dX*([0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]);
    IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);

    CamCalib.display.dX = dX;
    CamCalib.display.IP = IP;
    CamCalib.display.BASE = BASE;
    
    

    
    CamCalib.err_std = err_std;
    
    
    
end