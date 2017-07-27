function [cylinder_Model,R0,err,C0,V] = cylinderFit(P_E, est_C0, est_V, N_loop)
    if ~exist('N_loop','var')
        N_loop = 1;
    end
    P_E_sqr = sum(P_E.^2)';
    
    C0 = est_C0;
    V = est_V;
    P_E__V = P_E'*V;

    % assume V remain the same as rotation axis
    % only calculate C0
    for i=1:N_loop
        A(:,[1,2]) = 2* P_E([1,3],:)';
        A(:,3) = 1;
        b = P_E_sqr - P_E__V.^2 + 2*P_E__V*(C0'*V);
        res = inv(A'*A)*A'*b;

        R0 = sqrt(res(3)+ C0'*C0 - (C0'*V)^2);

        C0(1) = res(1);
        C0(3) = res(2);
        resid = sqrt(P_E_sqr - 2*P_E'*C0 + C0'*C0 - P_E__V.^2 + 2*P_E__V*(C0'*V) - (C0'*V)^2) - R0;
        err = displayError(resid,0);
    end
    
    P_Tire = convertEuclideanToCylinder(P_E, C0, V);
    Tire_Z_min = min(P_Tire(3,:));
    Tire_Z_max = max(P_Tire(3,:));
    cylinder_param = [C0+Tire_Z_max*V; C0+Tire_Z_min*V; R0]';
    cylinder_Model=cylinderModel(cylinder_param);
end

    