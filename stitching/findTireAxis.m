function findTireAxis(mat_file)
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');

    N_groove = param.stitching.N_groove;
    
    
    for groove_i = 1:N_groove
        groove_pixel_Cylinder_all{groove_i} = [];
    end
    for surface_i = 1:(N_groove+1)
        surface_pixel_Cylinder_all{surface_i} = [];
    end
    
    
%     Laser_Cylinder_all = [];
    Laser_E_all = [];

    for laser_i = 1:numel(laser_results)
        lr = laser_results(laser_i);
%         Laser_Cylinder_all = [Laser_Cylinder_all,lr.Laser_Cylinder];
        Laser_E_all = [Laser_E_all,lr.Laser_E];
        if isempty(lr.groove)
            continue;
        end
        for groove_i = 1:N_groove
            groove_pixel_Cylinder_all{groove_i} = ...
                [groove_pixel_Cylinder_all{groove_i}, lr.groove(groove_i).groove_pixels];
        end
        for surface_i = 1:(N_groove+1)
            surface_pixel_Cylinder_all{surface_i} = ...
                [surface_pixel_Cylinder_all{surface_i}, lr.surface(surface_i).surface_pixels];
        end
    end

    %% convert grooves and surfaces into euclidean coordinate system
    for groove_i = 1:N_groove
        groove_pixel_E_all{groove_i} = convertCylinderToEuclidean...
            (groove_pixel_Cylinder_all{groove_i}, C0, V, Laser_C_ref);
    end
    for surface_i = 1:(N_groove+1)
        surface_pixel_E_all{surface_i} = convertCylinderToEuclidean...
            (surface_pixel_Cylinder_all{surface_i}, C0, V, Laser_C_ref);
    end
    

        
    %% calculate rotation axis from groove 3
%     figure();hold on;
%     for groove_i = 1:N_groove
        groove_i = 3;
        [groove_ref_cylinder_Model,groove_ref_R,tire_axis.err,tire_axis.C0,tire_axis.V] = cylinderFit(groove_pixel_E_all{groove_i}, C0, V, 5);
%         plot(temp_cylinder_Model);
%         scatterPoints(groove_pixel_E_all{groove_i});
%         fprintf(['groove_',num2str(groove_i),': R=',num2str(R), ', C0=[,',num2str(C0'),'], ']);
%         displayError(err.resid);
%         fprintf('.\n');
%     end
    


    %% convert everything into tire coordinate system        
    
    Laser_Tire_all = convertEuclideanToCylinder(Laser_E_all,tire_axis.C0,tire_axis.V,Laser_C_ref);
    for groove_i = 1:N_groove
        groove_pixel_Tire_all{groove_i} = convertEuclideanToCylinder...
            (groove_pixel_E_all{groove_i},tire_axis.C0,tire_axis.V,Laser_C_ref);
    end
    for surface_i = 1:(N_groove+1)
        surface_pixel_Tire_all{surface_i} = convertEuclideanToCylinder...
            (surface_pixel_E_all{surface_i},tire_axis.C0,tire_axis.V,Laser_C_ref);
    end
    
    Tire_Z_min_all = min(Laser_Tire_all(3,:));
    Tire_Z_max_all = max(Laser_Tire_all(3,:));
    N_slices = 3000;
    N_th = max(20,size(Laser_Tire_all,2)/N_slices/5);
    Tire_Z_slice = linspace(Tire_Z_min_all,Tire_Z_max_all,N_slices+1);
    Tire_slices = [];
   
    Laser_Tire_groove_and_surface = [surface_pixel_Tire_all{:}];
    for slice_i = 1:N_slices
        Tire_Z_min_slice = Tire_Z_slice(slice_i);
        Tire_Z_max_slice = Tire_Z_slice(slice_i+1);

        temp_Tire = Laser_Tire_groove_and_surface;
        temp_ind = (temp_Tire(3,:)>Tire_Z_min_slice) &... 
            (temp_Tire(3,:)<Tire_Z_max_slice);
        temp_Tire_slice = temp_Tire(:,temp_ind);

        if (size(temp_Tire_slice,2)<N_th)
            temp_Tire = Laser_Tire_all;
            temp_ind = (temp_Tire(3,:)>Tire_Z_min_slice) &... 
                (temp_Tire(3,:)<Tire_Z_max_slice);
            temp_Tire_slice = temp_Tire(:,temp_ind);
        end

        if(size(temp_Tire_slice,2)>N_th)
            Tire_slices(slice_i).Z = (Tire_Z_min_slice+Tire_Z_max_slice)/2;
            Tire_slices(slice_i).R = mean(temp_Tire_slice(1,:));
        end
    end
    figure();
    scatter([Tire_slices.Z],[Tire_slices.R],3,'r','filled');axis equal;
    
    
    %%
    
%     Laser_E_all_pc = pointCloud(Laser_E_all');
% %     Laser_Cylinder_all_pc = pointCloud(Laser_Cylinder_all');
%     Laser_Tire_all_pc = pointCloud(Laser_Tire_all');
%     
%     for groove_i = 1:N_groove
% %         groove_pixel_Cylinder_all_pc{groove_i} = pointCloud(groove_pixel_Cylinder_all{groove_i}');
%         groove_pixel_E_all_pc{groove_i} = pointCloud(groove_pixel_E_all{groove_i}');
% %         groove_pixel_Tire_all_pc{groove_i} = pointCloud(groove_pixel_Tire_all{groove_i}');
%     end
%     for surface_i = 1:(N_groove+1)
% %         surface_pixel_Cylinder_all_pc{surface_i} = pointCloud(surface_pixel_Cylinder_all{surface_i}');
%         surface_pixel_E_all_pc{surface_i} = pointCloud(surface_pixel_E_all{surface_i}');
%         surface_pixel_Tire_all_pc{surface_i} = pointCloud(surface_pixel_Tire_all{surface_i}');
%     end
       

    figure();
    pcshow(Laser_E_all','r');
    hold on;
    pcshow(groove_pixel_E_all{3}','b');
    for surface_i = 1:(N_groove+1)
        pcshow(surface_pixel_E_all{surface_i}','g','MarkerSize',20);
    end
    plot(groove_ref_cylinder_Model)
    title('euclidean coordinate system');



    
    order = [2,3,1];
    figure();
    pcshow(Laser_Tire_all(order,:)','r');
    hold on;
    pcshow(groove_pixel_Tire_all{3}(order,:)','b','MarkerSize',20);
    for surface_i = 1:(N_groove+1)
        pcshow(surface_pixel_Tire_all{surface_i}(order,:)','g','MarkerSize',20);
    end
    title('Tire Cylinder coordinate system');
    fill3([max(groove_pixel_Tire_all{3}(2,:)),max(groove_pixel_Tire_all{3}(2,:)),min(groove_pixel_Tire_all{3}(2,:)),min(groove_pixel_Tire_all{3}(2,:))],...
        [min(groove_pixel_Tire_all{3}(3,:)),max(groove_pixel_Tire_all{3}(3,:)),max(groove_pixel_Tire_all{3}(3,:)),min(groove_pixel_Tire_all{3}(3,:))],...
        [groove_ref_R,groove_ref_R,groove_ref_R,groove_ref_R],'b');
    view(-30,40)    
    
    
    
    
    tire.laser.N_slices = N_slices;
    tire.laser.Tire_slices = Tire_slices;
    tire.laser.groove_ref_cylinder_Model = groove_ref_cylinder_Model;
    tire.laser.groove_ref_R = groove_ref_R;
    tire.laser.tire_axis = tire_axis;
    tire.laser.Z_range = [Tire_Z_min_all,Tire_Z_max_all];
    tire.laser.Laser_Tire_all = Laser_Tire_all;
    tire.laser.groove_pixel_Tire_all = groove_pixel_Tire_all;
    tire.laser.surface_pixel_Tire_all = surface_pixel_Tire_all;
    
    
    clear surface_i groove_i N_groove laser_i lr Laser_E_all
    clear groove_pixel_Cylinder_all surface_pixel_Cylinder_all
    clear groove_pixel_E_all surface_pixel_E_all
    clear groove_pixel_Tire_all surface_pixel_Tire_all
    clear groove_ref_cylinder_Model groove_ref_R Laser_Tire_all
    clear groove_pixel_Tire_all_pc surface_pixel_Tire_all_pc
    clear Laser_Tire_groove_and_surface
    clear temp_Tire temp_ind temp_Tire_slice slice_i N_slices N_th
    clear Tire_Z_min_all Tire_Z_max_all Tire_Z_slice Tire_Z_min_slice Tire_Z_max_slice
    clear Laser_E_all_pc groove_pixel_E_all_pc surface_pixel_E_all_pc
    clear laser_results Laser_Tire_all Tire_slices
    
    save(individual_mat_file.stitching,'-v7.3');
    clear tire;
    saveMatFile();
    

end


