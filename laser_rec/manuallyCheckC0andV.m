function manuallyCheckC0andV(mat_file, current_laser_i, C0_guess, V_guess)
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');
    
    
    current_lr = laser_results(current_laser_i);
    current_theta = current_lr.theta;
    im = getImFromPath(current_lr.im_path);

    t_im = im/quantile(im(:),0.99);
    t_im(t_im>1) = 1;
    current_xn = [];
    current_xn(1,:) = current_lr.Laser_X_c(1,:)./current_lr.Laser_X_c(3,:);
    current_xn(2,:) = current_lr.Laser_X_c(2,:)./current_lr.Laser_X_c(3,:);
    current_laser_x_pixel = inv_normalize_pixel(current_xn, CamCalib);   
    
    
    
    updateByC0andV(t_im, current_laser_x_pixel, current_theta, current_laser_i, laser_results, Laser_C_ref, CamCalib, C0_guess, V_guess)
    
    
    
    
    
end

function updateByC0andV(t_im, current_laser_x_pixel, current_theta, current_laser_i, laser_results, Laser_C_ref, CamCalib, C0_guess, V_guess)
    C0 = C0_guess;
    V = V_guess/sqrt(sum(V_guess.^2));

    Laser_Cylinder_all = [];
    for laser_i = 1:numel(laser_results)
        lr = laser_results(laser_i);
        theta = lr.theta;
        Laser_X_c = lr.Laser_X_c;    
        Laser_Cylinder = convertEuclideanToCylinder(Laser_X_c,C0,V,Laser_C_ref);
        Laser_Cylinder(2,:) = mod(Laser_Cylinder(2,:) + theta - current_theta, 360);
        Laser_Cylinder_all = [Laser_Cylinder_all, Laser_Cylinder];
    end
    
    Laser_X_c_all_afterRotation = convertCylinderToEuclidean(Laser_Cylinder_all,C0,V,Laser_C_ref);
    xn = [];
    xn(1,:) = Laser_X_c_all_afterRotation(1,:)./Laser_X_c_all_afterRotation(3,:);
    xn(2,:) = Laser_X_c_all_afterRotation(2,:)./Laser_X_c_all_afterRotation(3,:);
    laser_x_p_all_afterRotation = inv_normalize_pixel(xn, CamCalib);   
    
    ind = (laser_x_p_all_afterRotation(1,:)>1) & (laser_x_p_all_afterRotation(1,:)<size(t_im,2)) &...
        (laser_x_p_all_afterRotation(2,:)>1) & (laser_x_p_all_afterRotation(2,:)<size(t_im,1)) &...
        (Laser_X_c_all_afterRotation(3,:)<C0(3));
    
    figure();imshow(t_im);axis on;
    hold on; 
    scatterPoints(laser_x_p_all_afterRotation(:,ind), 1, 'r');
    scatterPoints(current_laser_x_pixel, 2, 'g');
    title(['laser_',num2str(current_laser_i),', theta=', num2str(current_theta),...
        ': C0=[',num2str(C0'),']',', V=[',num2str(V'),']']);
    
end


