function getStrongestFeatures(mat_file)
    load(mat_file);
    
    sigma = param.laserrec.sigma;
    dist_th = param.laserrec.dist_th;
    tire_region = param.tire_region;
    laser_region = param.laserrec.laser_region;
    
    N_strongest = param.laserrec.N_strongest;
    folder_name = folder_path;
    

    
    filter=fspecial('Gaussian',4*sigma+1,sigma);
    N_laser = numel(laser_file_list);

    b_display = 0;

    laser_results = [];
    tic;
    for laser_i = 1:N_laser
        SAT_PERCENT=0.01;

        im_path = [folder_name,'\',param.file_name.laserandps,'\',laser_file_list(laser_i).name];
        im = getImFromPath(im_path);
        
        [~, p_laser_x, Laser_X_c, features, validPoints, sat_val]...
            = getStrongestFeatureFromIm(im, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, b_display, figure_count, SAT_PERCENT);
        figure_count = figure_count + b_display;
        
        laser_results(laser_i).ind = laser_i;
        laser_results(laser_i).im_path = im_path;
        laser_results(laser_i).SAT_PERCENT = SAT_PERCENT;
        laser_results(laser_i).sat_val = sat_val;
        

        laser_results(laser_i).p_laser_x = p_laser_x;
        laser_results(laser_i).Laser_X_c = Laser_X_c;
        laser_results(laser_i).validPoints = validPoints;
        laser_results(laser_i).features = features;

        time = toc;
        fprintf(['working on ' num2str(laser_i),'/',num2str(N_laser), ...
            ', sat_val = ', num2str(sat_val),...
            ', time remain = ',num2str((N_laser-laser_i)*time/laser_i),'.\n']);

        clear SAT_PERCENT im_path im Laser_X_c features validPoints sat_val time
%         save(individual_mat_file.laserrec,'laser_results','-v7.3');
    end
    
    clear sigma dist_th tire_region laser_region N_strongest
    clear folder_name filter laser_i N_laser b_display
    save(individual_mat_file.laserrec,'-v7.3');
    clear laser_results;
    saveMatFile;
end

% function [filtered_im, Laser_X_c, features, validPoints, sat_val]...
%     = getStrongestFeatureFromIm(im, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, b_display, figure_count, SAT_PERCENT)
% 
%     sat_val = quantile(im(:),1-SAT_PERCENT);
%     temp_im = im2double(rgb2gray(im)) / sat_val;
%     temp_im(temp_im>1) = 1;
%     filtered_im = imfilter(temp_im,filter);
%     [p_laser_x,Laser_X_c] = getOneLaserScanResult(im,LaserCalib,CamCalib,...
%     laser_region,quantile(im(:),0.999),0);
% 
%     
%     points = detectSURFFeatures(filtered_im);
%     dist = pdist2(p_laser_x',points.Location);
%     dist_min = min(dist);
%     ind = dist_min<dist_th;
%     ind = ind' | (points.Location(:,1)<tire_region(1)) | (points.Location(:,1)>tire_region(3))...
%         | (points.Location(:,2)<tire_region(2)) | (points.Location(:,2)>tire_region(4));
%     points(ind) = [];
%     
%     
%     strongest = selectStrongest(points,N_strongest);
%     [features,validPoints] = extractFeatures(filtered_im,strongest); 
%     if (b_display);
%         figure(figure_count);
%         imshow(filtered_im);
%         hold on;
%         scatterPoints(validPoints.Location,5,'r');
%     end
% end


