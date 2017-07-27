function stitchPsRes( mat_file, ps_set_res )
    load(mat_file);
    load(individual_mat_file.stitching);
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');
%     map=jet(20);
%     figure();hold on;
%     P_E_after_theta_shift_all = [];
%     P_Tire_all = [];
    
    b_display = 1;

    stretched_factor = max([tire.laser.Tire_slices.R])/180*pi;
    
    tire.tire_region = [29,41,544,633];
    tire.stretched_factor = stretched_factor;
    
    tire_region = tire.tire_region;
    for ps_i = 1:20
        ps = ps_set_res(ps_i);
        theta_shift = ps.theta_shift;
        rerenderIm = ps.rerenderIm;
        
        [Z_E, Laser_E_all] = getZFromSurfaceNormal(ps, C0, V, ps.region_due_to_size_drop, CamCalib, laser_results, Laser_C_ref, 10, 1, 1, 0);
        x = ps.region_due_to_size_drop(1):ps.DROP_RATIO:ps.region_due_to_size_drop(3);
        y = ps.region_due_to_size_drop(2):ps.DROP_RATIO:ps.region_due_to_size_drop(4);
        [NUM_Y, NUM_X] = size(Z_E);
        [xx,yy]=meshgrid(x,y);
        xn = normalize_pixel([xx(:),yy(:)],CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
        X_E = reshape(xn(:,1),NUM_Y,NUM_X).*Z_E;
        Y_E = reshape(xn(:,2),NUM_Y,NUM_X).*Z_E;
        
        X_E = X_E(tire_region(2):tire_region(4),tire_region(1):tire_region(3));
        Y_E = Y_E(tire_region(2):tire_region(4),tire_region(1):tire_region(3));
        Z_E = Z_E(tire_region(2):tire_region(4),tire_region(1):tire_region(3));
        rerenderIm = rerenderIm(tire_region(2):tire_region(4),tire_region(1):tire_region(3),:);
        
        
        P_E = [X_E(:),Y_E(:),Z_E(:)];
        
        P_Cylinder = convertEuclideanToCylinder(P_E,C0,V,Laser_C_ref);
        P_Cylinder(2,:) = mod(P_Cylinder(2,:) + theta_shift + 180,360) - 180;

        P_E_after_theta_shift = convertCylinderToEuclidean(P_Cylinder,C0,V,Laser_C_ref);
        P_Tire = convertEuclideanToCylinder(P_E_after_theta_shift,tire_axis.C0,tire_axis.V,Laser_C_ref);

        P_Tire_stretched = P_Tire([2,3,1],:);
        P_Tire_stretched(1,:) = P_Tire_stretched(1,:)*stretched_factor;
        
        
%         P_E_after_theta_shift_all = [P_E_after_theta_shift_all, P_E_after_theta_shift];
%         P_Tire_all = [P_Tire_all,P_Tire];
        
        if (b_display)
            figure();hold on;pcshow(P_E,reshape(rerenderIm,numel(rerenderIm)/3,3));
            hold on;scatterPoints(Laser_E_all);
            title('P_E');view(90,180);

            figure();pcshow(P_Tire_stretched'); colorbar; colormap jet;
            caxis([tire.laser.groove_ref_R-1, max([tire.laser.Tire_slices.R]+2)]);
            title('P_Tire_stretched');view(0,-90);

            figure();imshow(rerenderIm/quantile(rerenderIm(:),0.99))        

    %         figure(1000);hold on;pcshow(P_E_after_theta_shift',reshape(rerenderIm,numel(rerenderIm)/3,3));
    %         figure();pcshow(P_Tire([2,3,1],:)'); colorbar; colormap colorcube;
    %         caxis([tire.laser.groove_ref_R-1, max([tire.laser.Tire_slices.R]+2)]);
    %         
    %         figure();pcshow(P_Tire',reshape(rerenderIm,numel(rerenderIm)/3,3));

        end
        
        
        
        tire.ps_set(ps_i).im_filepath = ps.im_filepath;
        tire.ps_set(ps_i).region = ps.region;
        tire.ps_set(ps_i).whiteboard_region = ps.whiteboard_region;
        tire.ps_set(ps_i).original_region = ps.original_region;
        tire.ps_set(ps_i).region_due_to_size_drop = ps.region_due_to_size_drop;        
        tire.ps_set(ps_i).DROP_RATIO = ps.DROP_RATIO;
        tire.ps_set(ps_i).stiching_tire_region = tire_region;
        tire.ps_set(ps_i).color_vec_map = ps.color_vec_map(tire_region(2):tire_region(4),tire_region(1):tire_region(3),:);
        tire.ps_set(ps_i).rerenderIm = ps.rerenderIm(tire_region(2):tire_region(4),tire_region(1):tire_region(3),:);
        tire.ps_set(ps_i).theta_shift = ps.theta_shift;
        tire.ps_set(ps_i).n = ps.n(tire_region(2):tire_region(4),tire_region(1):tire_region(3),:);
        tire.ps_set(ps_i).rho = ps.rho(tire_region(2):tire_region(4),tire_region(1):tire_region(3));
        
        tire.ps_set(ps_i).P_E = P_E;
        tire.ps_set(ps_i).Laser_E_all = Laser_E_all;
            
        tire.ps_set(ps_i).P_E_after_theta_shift = P_E_after_theta_shift;
        tire.ps_set(ps_i).P_Tire = P_Tire;
        tire.ps_set(ps_i).stretched_factor = stretched_factor;
        tire.ps_set(ps_i).P_Tire_stretched = P_Tire_stretched;
        
        clear ps theta_shift rerenderIm Z_E Laser_E_all x y NUM_Y NUM_X xx yy
        clear xn X_E Y_E P_E P_Cylinder P_E_after_theta_shift P_Tire P_Tire_stretched


    end           

    clear tire_region stretched_factor ps_set_res laser_results ps_i b_display
    save(individual_mat_file.stitching);
    clear tire;
    saveMatFile();

end

