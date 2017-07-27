function points_euclidean = convertCylinderToEuclidean...
    (P_cylinder, center, rotation_vector, point_that_indicate_0_phi_E)

    %% check parameters;
    
    C = center;
    V = rotation_vector;
    if (size(P_cylinder,1)~=3)
        if (size(P_cylinder,2)==3)
            P_cylinder = P_cylinder';
        else
            fprintf(2,'points dementions wrong');
            return;
        end
    end

    if ~( min(size(V)==[3,1]) )
        if ( min(size(V)==[1,3]) )
            V = V';
        else
            fprintf(2,'rotation_vector dementions wrong');
            return;
        end
    end

    if ~( min(size(C)==[3,1]) )
        if ( min(size(C)==[1,3]) )
            C = C';
        else
            fprintf(2,'center_coordinate dementions wrong');
            return;
        end
    end

    if (~exist('point_that_indicate_0_phi_E','var')),  
        fprintf('point_that_indicate_0_phi_E required');
    end


    %% calculate phi=0 direction

    new_0_phi = point_that_indicate_0_phi_E - C;
    Height_0_phi = V'*new_0_phi;
    new_0_phi_remain = new_0_phi - V*Height_0_phi;
    new_0_phi_direction = normalizeColVector(new_0_phi_remain);

    %%
    Rho = P_cylinder(1,:);
    Phi_degree = P_cylinder(2,:);
    Height = P_cylinder(3,:);
    
    new_P_direction = rotateAroundAxis(new_0_phi_direction,[0;0;0],V,-Phi_degree/180*pi);
    new_P_remain = new_P_direction .* repmat(Rho,3,1);
    new_P = new_P_remain + V*Height;
    points_euclidean = new_P + repmat(C,1,size(P_cylinder,2));

end

