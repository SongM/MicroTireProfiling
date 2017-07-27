function P_cylinder = convertEuclideanToCylinder...
    (points_euclidean, center, rotation_vector, point_that_indicate_0_phi_E)

    %% check parameters;
    
    P_E = points_euclidean;
    C = center;
    V = rotation_vector;
    if (size(P_E,1)~=3)
        if (size(P_E,2)==3)
            P_E = P_E';
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
        point_that_indicate_0_phi_E = P_E(:,1);  
    end


    %% calculate phi=0 direction
    new_0_phi = point_that_indicate_0_phi_E - C;
    Height_0_phi = V'*new_0_phi;
    new_0_phi_remain = new_0_phi - V*Height_0_phi;
    new_0_phi_direction = normalizeColVector(new_0_phi_remain);

    %%
    X = P_E(1,:)-C(1);
    Y = P_E(2,:)-C(2);
    Z = P_E(3,:)-C(3);
    new_P = [X;Y;Z];
    Height = V'*new_P;
    new_P_remain = new_P - V*Height;
    [new_P_direction,Rho]=normalizeColVector(new_P_remain);

    cosPhi = new_0_phi_direction'*new_P_direction;
    cosPhi(cosPhi>1) = 1;
    cosPhi(cosPhi<-1) = -1;
    sinPhi = cross( new_P_direction, repmat(new_0_phi_direction,1,size(new_P_direction,2)) )...
        ./ repmat(V,1,size(new_P_direction,2));
    sinPhi = sinPhi(2,:);
%     sign_Phi = sign((new_P_direction(1,:)-new_0_phi_direction(1)));

    sign_Phi = sign(sinPhi);
    Phi_degree = acos(cosPhi).*sign_Phi/pi*180;

    P_cylinder = [Rho;Phi_degree;Height];

end

