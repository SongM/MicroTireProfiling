function xyzC_tire_original_est = initTireSurface(param,CamCalib)
% xyzC = getWhiteBoardxyzC(plane_eqn,CamCalib,region)
    region = param.tire_region;

    R = param.tire.original_est.R;
    C0 = param.tire.original_est.C;
    V = param.tire.original_est.V;

    x = region(1):region(3);
    y = region(2):region(4);
    NUM_Y = numel(y);
    NUM_X = numel(x);
    [xx,yy] = meshgrid(x,y);
    [xn] = normalize_pixel([xx(:),yy(:)]',CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);


    tmp = xn(:,1);
    x = xn(1,:);
    y = xn(2,:);
    u = V(1);
    v = V(2);
    w = V(3);
    C0_1 = C0(1);
    C0_2 = C0(2);
    C0_3 = C0(3);

    xu_yv_w = x*u + y*v + w;
    c1u_c2v_c3w = C0'*V;

    a1 = x - xu_yv_w*u;
    a2 = y - xu_yv_w*v;
    a3 = 1 - xu_yv_w*w;
    b1 = c1u_c2v_c3w*u - C0_1;
    b2 = c1u_c2v_c3w*v - C0_2;
    b3 = c1u_c2v_c3w*w - C0_3;
    %% ax^2 + 2bx = c;
    %% x = -(b+sqrt(ac+b^2))/a
    a = a1.^2 + a2.^2 + a3.^2;
    b = a1*b1 + a2*b2 + a3*b3;
    c = R^2 - b1^2 - b2^2 - b3^2;
    Z = -(b+sqrt(c*a+b.^2))./a;
%     Z = -plane_eqn(4)./(plane_eqn(1)*xn(1,:)+plane_eqn(2)*xn(2,:)-1);
    xyzC_tire_original_est(:,:,1) = reshape(Z.*xn(1,:),NUM_Y,NUM_X);
    xyzC_tire_original_est(:,:,2) = reshape(Z.*xn(2,:),NUM_Y,NUM_X);
    xyzC_tire_original_est(:,:,3) = reshape(Z,NUM_Y,NUM_X);
end

