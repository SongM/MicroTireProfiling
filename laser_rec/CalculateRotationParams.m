function [laser_results, C0, V] = CalculateRotationParams(laser_results,...
    CamCalib, C0_est, V_est, N_loop, figure_count)
    
%     load(mat_file);
%     load(individual_mat_file.laserrec,'laser_results');
%     
    
%     C0_est = param.laserrec.C0_est;
%     V_est = param.laserrec.V_est;
%     N_loop = param.laserrec.N_loop;

    N_laser = numel(laser_results);
    
    laser_i = 1;
    laser_results(laser_i).d_theta = 0;
    laser_results(laser_i).V = V_est;
    laser_results(laser_i).C0 = C0_est;
    b_debug = 0;
    tic;
    for laser_i = 2:N_laser
        lr = laser_results(laser_i);
        Laser_X_c = lr.Laser_X_c;
        

        matchedPoints1 = lr.matchedPoints1;
        matchedPoints2 = lr.matchedPoints2;
        l1 = normalize_pixel(matchedPoints1.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
        l2 = normalize_pixel(matchedPoints2.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
        l1(:,3) = 1;
        l2(:,3) = 1;
        l1 = l1';
        l2 = l2';

        Z1_est = mean(Laser_X_c(3,:)).*ones(1,size(l1,2));
        Z2_est = Z1_est;
        
        [d_theta, V, C0, figure_count, Z1, Z2, l1, l2] =...
            getRotationAxisAndAngle(C0_est, V_est, Z1_est, Z2_est,...
            Laser_X_c, l1, l2, N_loop, b_debug, figure_count);
        
        t = toc;
        fprintf(['laser_',num2str(laser_i),'/',num2str(N_laser),', d_theta=',num2str(d_theta),...
            ', V=[', num2str(V'), '], C0=[',num2str(C0'),...
            '], N of matched points = ', num2str(matchedPoints2.Count),...
            ', time_remain=',num2str(t/(laser_i-1)*(N_laser-laser_i)),...
            '.\n']);
        laser_results(laser_i).d_theta = d_theta;
        laser_results(laser_i).V = V;
        laser_results(laser_i).C0 = C0;
    end
%     d_thetas = [laser_results.d_theta];
    C0 = mean([laser_results(2:end).C0], 2);
    V = mean([laser_results(2:end).V], 2);
    V = V/sqrt(sum(V.^2));

end

