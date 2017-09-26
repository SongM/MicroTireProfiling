function [tire_axis, T_groove_range_Z_ind] = findTireAxis(groove_i, laserrec, param, C2ObjRotAxis, Laser_C_ref)

%     N_groove = param.stitching.N_groove;
    dist_between_groove = 20;
   
%% get Tire C0 and V, and find the outlier for groove classification.
%     groove_i = 4;
    t_M_groove_bottom = [];
    t_M_groove = [];
    for laser_i = 1:numel(laserrec)
        lr = laserrec(laser_i);
        if isempty(lr.groove)
            continue;
        end
        groove_pixels = lr.groove(groove_i).groove_pixels;
        groove_bottom_R_th = min(groove_pixels(1,:))+0.1*(max(groove_pixels(1,:))-min(groove_pixels(1,:)));
        groove_bottom_ind = groove_pixels(1,:)<groove_bottom_R_th;
        t_M_groove_bottom = [t_M_groove_bottom, groove_pixels(:,groove_bottom_ind)];
        t_M_groove = [t_M_groove, groove_pixels];
    end
    
    ind = zeros(1,size(t_M_groove_bottom,2));
%     t_Z = t_Obj_groove'*C2ObjRotAxis.V;
    % find outlier of groove classification.
    while(1)
        t_Z = t_M_groove_bottom(3,~ind);
        mean_tt_Z = mean(t_Z);
        if (max(abs(t_Z-mean(t_Z)))<dist_between_groove)
            break;
        end
        ind = ind | (abs(t_M_groove_bottom(3,:)-mean_tt_Z)>std(t_Z));
    end
    % get tire_axis.V
    t_Obj_groove_bottom = convertCylinderToEuclidean(t_M_groove_bottom,C2ObjRotAxis,Laser_C_ref);
    [VV,VV_err] = linearFit(t_Obj_groove_bottom(:,~ind));
    VV_unit = VV/sqrt(sum(VV.^2))*sign(VV(2));
    t_Obj_groove_bottom_min = min(t_Obj_groove_bottom');
    t_Obj_groove_bottom_max = max(t_Obj_groove_bottom');
    t_Obj_groove_bottom_plane = [...
        t_Obj_groove_bottom_min(1),0,t_Obj_groove_bottom_min(3);...
        t_Obj_groove_bottom_min(1),0,t_Obj_groove_bottom_max(3);...
        t_Obj_groove_bottom_max(1),0,t_Obj_groove_bottom_max(3);...
        t_Obj_groove_bottom_max(1),0,t_Obj_groove_bottom_min(3);...
        t_Obj_groove_bottom_min(1),0,t_Obj_groove_bottom_min(3)];
    t_Obj_groove_bottom_plane(:,2) = (1-t_Obj_groove_bottom_plane(:,[1,3])*VV([1,3]))/VV(2);
    figure();scatterPoints(t_Obj_groove_bottom(:,~ind),5);axis equal;rotate3d on;
    hold on;fill3(t_Obj_groove_bottom_plane(:,1),t_Obj_groove_bottom_plane(:,2),t_Obj_groove_bottom_plane(:,3),...
        'b', 'FaceAlpha', 0.2);
    

    % get tire_axis.C0
    [groove_ref_Tire_Model,groove_ref_R,tire_axis_err,tire_axis] =...
        cylinderFit(t_Obj_groove_bottom(:,~ind), C2ObjRotAxis.t, VV_unit, 5, Laser_C_ref);
    
%     
    groove_ind = zeros(1,size(t_M_groove,2));
%     t_Z = t_Obj_groove'*C2ObjRotAxis.V;
    % find outlier of groove classification.
    while(1)
        t_Z = t_M_groove(3,~groove_ind);
        mean_tt_Z = mean(t_Z);
        if (max(abs(t_Z-mean(t_Z)))<dist_between_groove)
            break;
        end
        groove_ind = groove_ind | (abs(t_M_groove(3,:)-mean_tt_Z)>std(t_Z));
    end

    
    t_Obj_groove = convertCylinderToEuclidean(t_M_groove,C2ObjRotAxis,Laser_C_ref); 
    t_T_groove = convertEuclideanToCylinder(t_Obj_groove(:,~groove_ind), tire_axis, Laser_C_ref);

    min_t_groove_Tire_Z = min(t_T_groove(3,:));
    max_t_groove_Tire_Z = max(t_T_groove(3,:));
    
    fprintf(['groove_',num2str(groove_i), ', C0=[', num2str(tire_axis.t'),...
        '], VV=[',num2str(VV_unit'),'], groove_Z_range=[',...
        num2str([min_t_groove_Tire_Z,max_t_groove_Tire_Z]), '].\n']);
    T_groove_range_Z_ind = [min_t_groove_Tire_Z,max_t_groove_Tire_Z];
end


