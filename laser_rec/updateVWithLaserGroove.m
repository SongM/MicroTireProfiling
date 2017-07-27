function V = updateVWithLaserGroove(t_laser_results, N_groove, C0, V, Laser_C_ref, b_display)
    if ~exist('b_display','var')
        b_display = 1;
    end

    groove_C_all = cell(N_groove, 1);
    for laser_i = 1:numel(t_laser_results)
        lr = t_laser_results(laser_i);
        if (numel(lr.groove_C)==0)
            continue;
        end
        for groove_i = 1:N_groove
            groove_C_all{groove_i}=[groove_C_all{groove_i},[lr.groove_C(:,groove_i);lr.theta]];
        end
    end
    

    for groove_i = 2:(N_groove-1)
        fit_V(:,groove_i-1) = fitV(groove_C_all{groove_i}, C0, V, Laser_C_ref, b_display);
    end
    V = mean(fit_V,2);
    V = V/sqrt(sum(V.^2));



end

function fitted_V = fitV(t_C, C0, V, Laser_C_ref, b_display)
        
    t_C_theta = t_C(4,:);
    P = t_C(1:3,:)-repmat(C0,1,220);
    Cylinder_Z = P'*V;
    Cylinder_R = sqrt(sum((P-V*(P'*V)').^2));

    ind = zeros(numel(t_C_theta),1);
    while(max(abs(Cylinder_Z(~ind) - mean(Cylinder_Z(~ind))))>4*mean(abs(Cylinder_Z(~ind) - mean(Cylinder_Z(~ind)))))
        ind = ind | (abs(Cylinder_Z - mean(Cylinder_Z(~ind)))>4*mean(abs(Cylinder_Z(~ind) - mean(Cylinder_Z(~ind)))));
    end

    
    
        if(b_display)
            figure();
            subplot(1,5,1);scatter(t_C(4,~ind),t_C(1,~ind),5,'r','filled');
            subplot(1,5,2);scatter(t_C(4,~ind),t_C(2,~ind),5,'r','filled');
            subplot(1,5,3);scatter(t_C(4,~ind),t_C(3,~ind),5,'r','filled');
            subplot(1,5,4);scatter(t_C(4,~ind),Cylinder_Z(~ind),5,'r','filled');title('Cylinder_Z');
            subplot(1,5,5);scatter(t_C(4,~ind),Cylinder_R(~ind),5,'r','filled');title('Cylinder_Z');
        end
        
        V=[0;1;0.4];
        V=V/sqrt(sum(V.^2));
        
        P_Cylinder1 = convertEuclideanToCylinder(t_C(1:3,:), C0, V, Laser_C_ref);
        P_Cylinder1(2,:) = P_Cylinder1(2,:) + t_C_theta;
        P_Cylinder1(:,ind) = [];
        
        figure(); scatterPoints(P_Cylinder1([2,3],:),5);title(num2str(V'));
       
       
       
        P_Cartesian1 = convertCylinderToEuclidean(P_Cylinder1,C0,V,Laser_C_ref);
        
        
        delta_V = inv(P_Cartesian_Cylinder1*P_Cartesian_Cylinder1')*P_Cartesian_Cylinder1*ones(size(P_Cartisian_Cylinder1,2),1);
        
        V1 = V;
        delta_Z1 = P'*V1 - mean(P'*V1);
        std(delta_Z1)
        delta_V = inv(P*P')*P*delta_Z;
        

        t_E = convertCylinderToEuclidean(t_C(:,~ind),C0,V,Laser_C_ref);
        [plane_param,plane_err] = linearFit(t_E);
        fitted_V = plane_param/sqrt(sum(plane_param.^2))*sign(plane_param(2));
end



