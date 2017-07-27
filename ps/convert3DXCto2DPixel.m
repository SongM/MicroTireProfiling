function x_pixel = convert3DXCto2DPixel(X_C,CamCalib)
    xn(1,:) = X_C(1,:) ./ X_C(3,:);
    xn(2,:) = X_C(2,:) ./ X_C(3,:);
    x_pixel = inv_normalize_pixel(xn,CamCalib);

end