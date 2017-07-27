function seperateGrooveAndSurfacePixel(mat_file)
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');
    
    b_debug = 0;
%%
    reference_groove_idx = 3;
%     reference_groove_region_Z = [-15,0];
    reference_groove_region_Z = param.stitching.ref_groove_region_Z;
    tire_surface_region_Z = param.stitching.tire_surface_region_Z;
    groove_region_half_size_Z = param.stitching.groove_region_half_size_Z;
    N_groove = param.stitching.N_groove;
    fitting_err_max_th = param.stitching.fitting_err_max_th;
    if(0)
        laser_i = 1;
        lr = laser_results(laser_i);
        im = getImFromPath(lr.im_path);
        
        Laser_Cylinder = double(lr.Laser_Cylinder);

        Rho = Laser_Cylinder(1,:);
        Z = Laser_Cylinder(3,:);
        Rho_min = min(Rho);
        Rho_max = max(Rho);
        
        figure();subplot(1,2,1);scatter(Z,Rho,2,'r','filled');axis equal;
        rectangle('Position',[tire_surface_region_Z(1),Rho_min,tire_surface_region_Z(2)-tire_surface_region_Z(1),Rho_max-Rho_min],'EdgeColor','g');
        rectangle('Position',[reference_groove_region_Z(1),Rho_min,reference_groove_region_Z(2)-reference_groove_region_Z(1),Rho_max-Rho_min],'EdgeColor','b');
        
        x_pixel = convert3DXCto2DPixel(lr.Laser_X_c,CamCalib);
        subplot(1,2,2);imshow(im/quantile(im(:),0.99));hold on;scatter(x_pixel(1,:),x_pixel(2,:),2,'g','filled');
    end
%%
    for laser_i = 1:numel(laser_results)
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf(['     laser_',num2str(laser_i),':\n']);

        lr = laser_results(laser_i);

        Laser_Cylinder = double(lr.Laser_Cylinder);

        Rho = Laser_Cylinder(1,:);
        Z = Laser_Cylinder(3,:);
        
        [groove_Rho,groove_Z] = findGroovePostion(Rho,Z,N_groove,tire_surface_region_Z);
        if (numel(groove_Rho)<N_groove)
            fprintf(2,'xxxxxxxxxxxxxxxxxxxxxx.......error.......xxxxxxxxxxxxxxxxxxxxxx\n'); continue;
        end
        [groove_Z, groove_ind] = sort(groove_Z);
        groove_Rho = -groove_Rho(groove_ind);
        clear groove_ind;
%%
        groove = getGrooveDetails(Laser_Cylinder, Rho, Z, N_groove, groove_Rho, groove_Z, groove_region_half_size_Z);
        mean_ref_groove_Z = mean(groove(reference_groove_idx).groove_pixels(3,:));
        if ( mean_ref_groove_Z < reference_groove_region_Z(1) |...
            mean_ref_groove_Z > reference_groove_region_Z(2) )
            fprintf(2,'xxxxxxxxx.......error, groove 3 is not in the region.......xxxxxxxxx\n'); continue;
        end
        clear mean_ref_groove_Z

%%
        surface = getSurfaceDetails(Laser_Cylinder, Rho, Z, N_groove, ...
            tire_surface_region_Z, groove_Z,...
            fitting_err_max_th, groove_region_half_size_Z, b_debug);
        
        laser_results(laser_i).groove = groove;
        laser_results(laser_i).surface = surface;


        if (b_debug)
            figure();
            displaySeperationResults(Z, Rho, N_groove, groove_Z, groove_Rho, groove, surface);
        end
        
        clear lr Laser_Cylinder groove surface Z Rho groove_Z groove_Rho
    end
    clear laser_i b_debug tire_surface_region_Z groove_region_half_size_Z;
    clear reference_groove_idx reference_groove_region_Z N_groove fitting_err_max_th
    save(individual_mat_file.laserrec);
    clear laser_results;
    saveMatFile
end

function [groove_Rho,groove_Z] = findGroovePostion(Rho, Z, N_groove, tire_surface_region_Z)
    [sorted_Z, sorted_ind] = sort(Z);    
    sorted_Rho = Rho(sorted_ind);

    sorted_Rho(sorted_Z<tire_surface_region_Z(1)) =[];
    sorted_Z(sorted_Z<tire_surface_region_Z(1)) =[];
    sorted_Rho(sorted_Z>tire_surface_region_Z(2)) =[];
    sorted_Z(sorted_Z>tire_surface_region_Z(2)) =[];
    if (numel(sorted_Z)~=numel(unique(sorted_Z)))
        groove_Rho = [];groove_Z = [];
        return;
    end
    [groove_Rho,groove_Z] = findpeaks(-sorted_Rho,sorted_Z,...
            'NPeaks',N_groove,'MinPeakProminence',4,'SortStr','descend');
end

function groove = getGrooveDetails(Laser_Cylinder, Rho, Z, N_groove,...
    groove_Rho, groove_Z, groove_region_half_size_Z)
    for groove_i = 1:N_groove
        groove_ind = (Z>(groove_Z(groove_i)-groove_region_half_size_Z))...
            & (Z<(groove_Z(groove_i)+groove_region_half_size_Z));
        mean_groove_Rho = mean(Rho(groove_ind));
        groove_ind = groove_ind & (Rho<mean_groove_Rho);
        groove(groove_i).ind = groove_i;
        groove(groove_i).peak = [groove_Rho(groove_i),groove_Z(groove_i)];
        groove(groove_i).groove_pixels = Laser_Cylinder(:,groove_ind);
        fprintf(['groove_',num2str(groove_i),', peak=[',num2str(groove(groove_i).peak),...
            '], groove_Rho_th=', num2str(mean_groove_Rho),...
            ', mean_groove_Rho=,' num2str(mean(groove(groove_i).groove_pixels(1,:))),'.\n']);
    end
end

function surface = getSurfaceDetails(Laser_Cylinder, Rho, Z, N_groove, ...
    tire_surface_region_Z, groove_Z,...
    fitting_err_max_th, groove_region_half_size_Z, b_debug)

    surface_edge_Z = [tire_surface_region_Z(1),groove_Z,tire_surface_region_Z(2)];
    for surface_i = 1:(N_groove+1)
        if (b_debug)
            fprintf(['surface_',num2str(surface_i),':\n']);
        end
        surface_region_Z = surface_edge_Z([surface_i,surface_i+1]) - groove_region_half_size_Z*[-1,1];
        temp_surface_pixel_ind = Z>surface_region_Z(1)...
            & (Z<surface_region_Z(2));
%             Rho_region_min = min(Rho(temp_surface_pixel_ind));
%             Rho_region_max = max(Rho(temp_surface_pixel_ind));

        temp_surface_pixel_ind = temp_surface_pixel_ind & Rho>mean(Rho(temp_surface_pixel_ind));

        tsp_Rho = Rho(temp_surface_pixel_ind);
        tsp_Z = Z(temp_surface_pixel_ind);
        while(numel(tsp_Rho)>1)
            [~, fitting_err] = linearFit([tsp_Rho;tsp_Z]);
            ind = abs(fitting_err.resid)>fitting_err.std;
            tsp_Rho(ind)=[];
            tsp_Z(ind)=[];
            if(b_debug)
                fprintf(['N_pixel=', num2str(numel(tsp_Rho)), ...
                    ', fitting_err: std=', num2str(fitting_err.std), ', max=', num2str(fitting_err.max),'.\n']);
            end
            if ((fitting_err.std<0.2)||(fitting_err.max<fitting_err_max_th))
                break;
            end
        end

        fitting_line_param = linearFit([tsp_Rho;tsp_Z]);
        dist2FittingLine = abs(fitting_line_param'*[Rho;Z]-1)/sqrt(sum(fitting_line_param.^2));
        surface_pixel_ind = temp_surface_pixel_ind & (dist2FittingLine<fitting_err_max_th);


        surface(surface_i).ind = surface_i;
        surface(surface_i).surface_region_Z = surface_region_Z;
        surface(surface_i).fitting_line_param = fitting_line_param;
        surface(surface_i).fitting_line_end_points = [surface_region_Z;...
            (1-fitting_line_param(2)*surface_region_Z)/fitting_line_param(1)];
        surface(surface_i).surface_pixels = Laser_Cylinder(:,surface_pixel_ind);

    end
end

function displaySeperationResults(Z, Rho, N_groove, groove_Z, groove_Rho, groove, surface)
    scatter(groove_Z,groove_Rho,'b+');
    hold on;scatter(Z,Rho,2,'r','filled');
    axis equal;

    for groove_i = 1:N_groove
        groove_pixels = groove(groove_i).groove_pixels;
        scatter(groove_pixels(3,:),groove_pixels(1,:),1,'g','filled');
    end

    for surface_i = 1:(N_groove+1)
        surface_pixels = surface(surface_i).surface_pixels;
        fitting_line_end_points = surface(surface_i).fitting_line_end_points;
        plot(fitting_line_end_points(1,:),fitting_line_end_points(2,:),'b');
        scatter(surface_pixels(3,:),surface_pixels(1,:),1,'g','filled');
    end
end