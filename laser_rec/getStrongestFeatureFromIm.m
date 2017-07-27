function [filtered_im, p_laser_x ,Laser_X_c, features, validPoints, sat_val]...
    = getStrongestFeatureFromIm(im, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, b_display, figure_count, SAT_PERCENT)

    sat_val = quantile(im(:),1-SAT_PERCENT);
    temp_im = im2double(rgb2gray(im)) / sat_val;
    temp_im(temp_im>1) = 1;
    filtered_im = imfilter(temp_im,filter);
    [p_laser_x,Laser_X_c] = getOneLaserScanResult(im,LaserCalib,CamCalib,...
    laser_region,quantile(im(:),0.999),0);

    
    points = detectSURFFeatures(filtered_im);
    dist = pdist2(p_laser_x',points.Location);
    dist_min = min(dist);
    ind = dist_min<dist_th;
    ind = ind' | (points.Location(:,1)<tire_region(1)) | (points.Location(:,1)>tire_region(3))...
        | (points.Location(:,2)<tire_region(2)) | (points.Location(:,2)>tire_region(4));
    points(ind) = [];
    
    
    strongest = selectStrongest(points,N_strongest);
    [features,validPoints] = extractFeatures(filtered_im,strongest); 
    if (b_display);
        figure(figure_count);
        imshow(filtered_im);
        hold on;
        scatterPoints(validPoints.Location,5,'r');
    end
end