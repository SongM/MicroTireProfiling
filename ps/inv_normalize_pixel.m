function [x_kk] = inv_normalize_pixel(xn, CamCalib)

    fc = CamCalib.fc;
    cc = CamCalib.cc;
    kc = CamCalib.kc;
    alpha_c = CamCalib.alpha_c;
    k1 = kc(1);
    k2 = kc(2);
    k3 = kc(5);
    p1 = kc(3);
    p2 = kc(4);


    flag_transpose = false;
    if (size(xn,1)~=2)
        xn = xn';
        flag_transpose = true;
    end
    
    n = size(xn,2);
    xn_x = xn(1,:);
    xn_y = xn(2,:);
    r_2 = xn_x.^2 + xn_y.^2;

    k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
    delta_x = [2*p1*xn_x.*xn_y + p2*(r_2 + 2*xn_x.^2);...
        p1 * (r_2 + 2*xn_y.^2)+2*p2*xn_x.*xn_y];

    x_kk = (xn .* (ones(2,1)*k_radial) + delta_x).* (fc * ones(1,n))+ cc*ones(1,n);

    if(flag_transpose)
        x_kk = x_kk';
    end
end
