function X_c = get3DCoorFromPlaneEqnAnd2DCoor(c_plane_eqn, p_x, CamCalib)
    % plane eqn: AX + BY + CZ = 1;
    % i_x_undistort (x,y)
    fc = CamCalib.fc;
    cc = CamCalib.cc;
    kc = CamCalib.kc;
    alpha_c = CamCalib.alpha_c;
%     kc = zeros(5,1);
    
    i_x_undistort = normalize_pixel(p_x,fc,cc,kc,alpha_c);
    temp_x = [i_x_undistort;ones(1,size(i_x_undistort,2))];
    Xc_Z = 1./(c_plane_eqn'*temp_x);
    Xc_X = i_x_undistort(1,:) .* Xc_Z;
    Yc_Y = i_x_undistort(2,:) .* Xc_Z;
    X_c = [Xc_X; Yc_Y; Xc_Z];
end

