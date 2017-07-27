function new_laser_results = findLaserGroove...
    (laser_results, C0_est, V_est, CamCalib, tire_surface_region_Z, N_groove, Laser_C_ref)
    new_laser_results = laser_results;
    N_laser = numel(laser_results);
    theta_all = cumsum([laser_results.d_theta])/pi*180;
    for laser_i=1:N_laser
        lr = laser_results(laser_i);
        Laser_X_c = lr.Laser_X_c;
        Laser_Cylinder = double(convertEuclideanToCylinder(Laser_X_c, C0_est, V_est, Laser_C_ref));
        groove_C = findGroovePostion(Laser_Cylinder, Laser_X_c, N_groove, tire_surface_region_Z);
        new_laser_results(laser_i).groove_C = groove_C;
        new_laser_results(laser_i).theta = theta_all(laser_i);
        new_laser_results(laser_i).Laser_Cylinder = Laser_Cylinder;
    end
end


function groove_C = findGroovePostion(Laser_Cylinder, Laser_X_c, N_groove, tire_surface_region_Z)

    groove_Cylinder=[];
    groove_C = [];
    
    Z = Laser_Cylinder(3,:);
    [sorted_Z, sorted_ind] = sort(Z);    
    sorted_Laser_Cylinder = Laser_Cylinder(:,sorted_ind);
    sorted_Laser_Cylinder(:,sorted_Z<tire_surface_region_Z(1)) =[];
    sorted_Z(sorted_Z<tire_surface_region_Z(1)) =[];
    sorted_Laser_Cylinder(:,sorted_Z>tire_surface_region_Z(2)) =[];
    sorted_Z(sorted_Z>tire_surface_region_Z(2)) =[];
    
    
    if (numel(sorted_Z)~=numel(unique(sorted_Z)))
        return;
    end
    [groove_Rho,groove_Z] = findpeaks(-sorted_Laser_Cylinder(1,:),sorted_Z,...
            'NPeaks',N_groove,'MinPeakProminence',4,'SortStr','descend');
    if (numel(groove_Z)<N_groove)
        return;
    end
    [groove_Z, groove_ind] = sort(groove_Z);
    groove_Rho = -groove_Rho(groove_ind);
    for groove_i=1:N_groove
        groove_ind = find(Z==groove_Z(groove_i));
        groove_Cylinder = [groove_Cylinder, Laser_Cylinder(:,groove_ind)];
        groove_C = [groove_C, Laser_X_c(:,groove_ind)];
    end
    
    if(sum(groove_Cylinder(1,:)~=groove_Rho))
        groove_C = [];
        return;
    end
end


