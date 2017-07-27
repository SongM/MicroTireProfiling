% clear all;close all;tic;

% folder_name = '20170208';
function [bs_withExtParam,error_record] = getBoardExtParam(bs, CamCalib, error_scale, wintx_default, winty_default)



    if (~exist('wintx_default','var')),  
        fprintf(2,'wintx_default is assigned as wintx_default = 0.\n')
        wintx_default = 0;  
    end
    if (~exist('winty_default','var')),  
        fprintf(1,'winty_default is assigned as winty_default = wintx_default.\n')
        winty_default = wintx_default;  
    end
    if (~exist('error_scale','var')),  
        fprintf(1,'error_scale is assigned as error_scale = 2.\n')
        error_scale = 2;
    end
    
    wintx = wintx_default;
    winty = winty_default;
    n_sq_x = bs(1).checkerboard.nx-2;
    n_sq_y = bs(1).checkerboard.ny-2;
    dX = bs(1).checkerboard.dX;
    dY = bs(1).checkerboard.dY;
    map = gray(256);
    
    thresh_cond = 1000000;
    
    kc = CamCalib.kc;
    fc = CamCalib.fc;
    cc = CamCalib.cc;
    alpha_c = CamCalib.alpha_c;
    err_std_calib = CamCalib.err_std;
%     
%     addpath('functions');
%     addpath('preprocessed_data');
%     folder_path = ['../',folder_name,'/'];
%     addpath(folder_path);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    board_proc = 1:numel(bs);
%     board_proc = 1:10;
    for kk = board_proc
        board = bs(kk);
        ima_name = board.filepath;
        fprintf(1, ['Loading board ', num2str(kk),...
            '/',num2str(numel(board_proc)),...
            ', BoardIdx = ', num2str(board.BoardIdx),...
            ', filepath = ']);
        disp(ima_name);
        I = double(getImFromPath(ima_name));
        I = 0.299 * I(:,:,1) + 0.5870 * I(:,:,2) + 0.114 * I(:,:,3);
        [x,X] = getBoardCorners(board, I, 0,...
            wintx, winty, n_sq_x, n_sq_y, dX, dY, map);

        bs(kk).x = x;
        bs(kk).X = X;
        eval(['x_', num2str(kk), ' = x;']);
        eval(['X_', num2str(kk), ' = X;']);
    end
%% check desactive


n_ima = numel(board_proc);
active_images=ones(1,n_ima);
check_cond = 1;

desactivated_images = [];
comp_ext_calib;
% if  ~isempty(desactivated_images)
%     fprintf(2,'some boards in board_calib is not good\n');
% end    


for kk = board_proc
    
    eval(['X_kk = X_', num2str(kk), ';']);
    eval(['omc_kk = omc_', num2str(kk), ';']);
    Rc_kk = rodrigues(omc_kk);
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

    Z_A=[Xc_kk(1:2,:)',ones(size(Xc_kk,2),1)];
    Z_b=Xc_z';
    plane_eqn_c = Z_A\Z_b;
    Z_resid = Z_A*plane_eqn_c - Z_b;
    Z_SSresid = sum(Z_resid.^2);
    Z_SStotal = (length(Xc_z')-1) * var(Xc_z);
    Z_rsq = 1 - Z_SSresid/Z_SStotal;
    if (Z_rsq<0.9)
        fprintf(2,'Z_rsp<0.9\n');
        desactivated_images = [desactivated_images , kk];

    end

    

    eval(['x_kk = x_', num2str(kk), ';']);

    ex_kk = y_kk - x_kk;

    bs(kk).Xc = Xc_kk;
    
    bs(kk).xp = y_kk;
    bs(kk).proj_error = ex_kk;

    bs(kk).omc = omc_kk;
    bs(kk).Rc = Rc_kk;
    bs(kk).Tc = Tc_kk;
    bs(kk).Tc = Tc_kk;
    bs(kk).plane_eqn_c = plane_eqn_c;
    
end
 
error_record=[];
% s=2;
bs_withExtParam = [];
for kk=board_proc
    board = bs(kk);
    error_record(kk,:)=[std(board.proj_error'),sum( std(board.proj_error')'>error_scale*err_std_calib) ];
    if (  sum( std(board.proj_error')' > error_scale*err_std_calib )  )
        desactivated_images = [desactivated_images , kk];
    end
    b = [];
    b.BoardIdx = board.BoardIdx;
    b.filepath = board.filepath;
    
    % other info
    b.otherInfo.energy = board.energy;
    b.otherInfo.idx = board.idx;
    b.checkerboard = board.checkerboard;
    
    % coor in coordinate systems.
    b.p_x = board.x;
    b.p_x_proj = board.xp;
    b.b_X = board.X;
    b.c_X = board.Xc;

    % board plane to camera plane
    b.b2c_omc = board.omc;
    b.b2c_Rc = board.Rc;
    b.b2c_Tc = board.Tc;
    
    % pixel plane
    b.pixel_plane.v1 = board.v1;
    b.pixel_plane.v2 = board.v2;
    b.pixel_plane.coor_x = board.coor_x;
    b.pixel_plane.coor_y = board.coor_y;
    b.pixel_plane.corner = [board.cornerx;board.cornery];
    
    % camera plane
    n_grid_x = b.checkerboard.nx -1;
    n_grid_y = b.checkerboard.ny -1;
        % grid
    b.camera_plane.coor_X =...
        fliplr(reshape(b.c_X(1,:),n_grid_x,n_grid_y));
    b.camera_plane.coor_Y =...
        fliplr(reshape(b.c_X(2,:),n_grid_x,n_grid_y));
    b.camera_plane.coor_Z =...
        fliplr(reshape(b.c_X(3,:),n_grid_x,n_grid_y));

        % Vector
    c_grid(:,:,1) = b.camera_plane.coor_X;
    c_grid(:,:,2) = b.camera_plane.coor_Y;
    c_grid(:,:,3) = b.camera_plane.coor_Z;

    grid_V1_1 = c_grid(:,1:end-1,:);
    grid_V1_2 = c_grid(:,2:end,:);
    V1_all = grid_V1_2 - grid_V1_1;
    c_V1 = reshape(mean(mean(V1_all)),1,3);
    b.camera_plane.V1 = c_V1;

    grid_V2_1 = c_grid(1:end-1,:,:);
    grid_V2_2 = c_grid(2:end,:,:);
    V2_all = grid_V2_2 - grid_V2_1;
    c_V2 = reshape(mean(mean(V2_all)),1,3);
    b.camera_plane.V2 = c_V2;

    
        % Corner
    if (b.BoardIdx == 0)
        n1 = 1;
    else
        n1 = 3;
    end
    c_corner(1,:) = reshape(c_grid(1,1,:),1,3) - c_V1 - n1*c_V2;
    c_corner(2,:) = reshape(c_grid(1,end,:),1,3) + c_V1 - n1*c_V2;
    c_corner(3,:) = reshape(c_grid(end,end,:),1,3) + c_V1 + c_V2;
    c_corner(4,:) = reshape(c_grid(end,1,:),1,3) - c_V1 + c_V2;

    b.camera_plane.Corner = c_corner';

    [b.camera_plane.plane_eqn,plane_err] = linearFit(b.c_X);
    if (plane_err.std>0.1)
        disp(['board ', num2str(kk),' plane fitting plane_err.std > 0.1.']);
    end
    
    % error
    b.error.p_x_proj_error = board.proj_error;
    bs_withExtParam = [bs_withExtParam,b];
    
%     figure(kk);
%     displayBoard(rot_board_set(kk));
end
desactivated_images_laser = unique(desactivated_images);

% board_laser_3D_rec = rot_board_set;
bs_withExtParam(desactivated_images_laser) = [];

end