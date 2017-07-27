function laserReconstruct(mat_file)
    load(mat_file);
    task = '3.A.1.GetRotationAxisAndAngle';
    displayTask(task,3);
    
    
    tire_region = param.tire_region;

    
    sigma = param.laserrec.sigma;
    dist_th = param.laserrec.dist_th;
    laser_region = param.laserrec.laser_region;
    N_strongest = param.laserrec.N_strongest;
    C0_est = param.laserrec.C0_est;
    V_est = param.laserrec.V_est;
    N_loop = param.laserrec.N_loop;

    
    b_debug = 0;
    SAT_PERCENT=0.01;
    
    c = clock;
    disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
    timestamp = ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5))];
    individual_mat_file.laserrec = [data_folder, '\', case_name, '\processed_data\', 'laserRec_result',timestamp,'.mat'];
    clear c timestamp;

    
    
    folder_name = folder_path;
    filter=fspecial('Gaussian',4*sigma+1,sigma);
    N_laser = numel(laser_file_list);
    laser_results = [];
    laser_i = 106;
        im1 = getImFromPath([folder_name,'\',param.file_name.laserandps,'\',laser_file_list(laser_i).name]);
        [filtered_im1, Laser_X_c1, features1, validPoints1]...
        = getStrongestFeatures(im1, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, 0, figure_count, SAT_PERCENT);
    
        laser_results(laser_i).ind = laser_i;
        laser_results(laser_i).im_path = [folder_name,'\',param.file_name.laserandps,'\',laser_file_list(laser_i).name];
        laser_results(laser_i).Laser_X_c = Laser_X_c1;
        laser_results(laser_i).d_theta = 0;
        laser_results(laser_i).V = V_est;
        laser_results(laser_i).C0 = C0_est;
        laser_results(laser_i).NofMatchedPoints_withPreviousFrame = N_strongest;

        fprintf(['laser_',num2str(laser_i),'/',num2str(numel(laser_file_list)),', d_theta=0',...
            ', V_est=[', num2str(V_est'), '], C0_est=[',num2str(C0_est'),...
            '], N_loop=',num2str(N_loop),'.\n']);
    
        save(individual_mat_file.laserrec,'laser_results','-v7.3');


    tic;
    while (laser_i < numel(laser_file_list))
        laser_i = laser_i + 1;
        im2 = getImFromPath([folder_name,'\',param.file_name.laserandps,'\',laser_file_list(laser_i).name]);
        [filtered_im2, Laser_X_c2, features2, validPoints2]...
            = getStrongestFeatures(im2, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, 0, figure_count, SAT_PERCENT);

        indexPairs = matchFeatures(features1,features2);
        matchedPoints1 = validPoints1(indexPairs(:,1),:);
        matchedPoints2 = validPoints2(indexPairs(:,2),:);
%         figure(); showMatchedFeatures(filtered_im1,filtered_im2,matchedPoints1,matchedPoints2);

        [l1, l2, matchedPoints1, matchedPoints2] = getRidOfWrongMatches(matchedPoints1, matchedPoints2, CamCalib);

        if(b_debug)
            figure_count = figure_count + 1;    figure(figure_count); 
            showMatchedFeatures(filtered_im1,filtered_im2,matchedPoints1,matchedPoints2);
            title(['laser_',num2str(laser_i)]);
        end
        
        features1 = features2;
        validPoints1 = validPoints2;
        filtered_im1 = filtered_im2; 

        
        Z1_est = mean(Laser_X_c2(3,:)).*ones(1,size(l1,2));
        Z2_est = Z1_est;
        [d_theta, V, C0, figure_count, Z1, Z2, l1, l2] =...
            getRotationAxisAndAngle(C0_est, V_est, Z1_est, Z2_est,...
            Laser_X_c2, l1, l2, N_loop, b_debug, figure_count);
        
        t = toc;
        fprintf(['laser_',num2str(laser_i),'/',num2str(numel(laser_file_list)),', d_theta=',num2str(d_theta),...
            ', V=[', num2str(V'), '], C0=[',num2str(C0'),...
            '], N of matched points = ', num2str(matchedPoints2.Count),...
            ', time_remain=',num2str(t/(laser_i-1)*(N_laser-laser_i)),...
            '.\n']);
        
        laser_results(laser_i).ind = laser_i;
        laser_results(laser_i).im_path = [folder_name,'\',param.file_name.laserandps,'\',laser_file_list(laser_i).name];
        laser_results(laser_i).Laser_X_c = Laser_X_c2;
        laser_results(laser_i).d_theta = d_theta;
        laser_results(laser_i).V = V;
        laser_results(laser_i).C0 = C0;
        laser_results(laser_i).NofMatchedPoints_withPreviousFrame = matchedPoints2.Count;

        save(individual_mat_file.laserrec,'laser_results','-v7.3');

    end
    clear SAT_PERCENT
    clear d_theta C0 V l1 l2 Z1 Z2 V_est C0_est Z1_est Z2_est
    clear b_debug dist_th sigma filter N_loop N_strongest t
    clear laser_file_list laser_i N_laser Laser_X_c1 Laser_X_c2
    clear features1 features2 im1 im2 filtered_im1 filtered_im2 
    clear indexPairs validPoints1 validPoints2 matchedPoints1 matchedPoints2
    clear laser_region tire_region
    
    job_done = [job_done;task];
    clear task;

    
    save(individual_mat_file.laserrec,'-v7.3');

    clear laser_results
    saveMatFile();
    
end

function [filtered_im, Laser_X_c, features, validPoints]...
    = getStrongestFeatures(im, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, b_display, figure_count, SAT_PERCENT)

    temp_im = im2double(rgb2gray(im)) / quantile(im(:),1-SAT_PERCENT);
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


function [l1, l2, matchedPoints1, matchedPoints2] = getRidOfWrongMatches...
    (matchedPoints1, matchedPoints2, CamCalib)

    l1 = normalize_pixel(matchedPoints1.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
    l2 = normalize_pixel(matchedPoints2.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
    
 
    temp = l1-l2;
    
    temp_l = sqrt(sum(temp.^2,2));
    temp_angle = atan2(temp(:,1),temp(:,2));
    ind = (temp_l > (median(temp_l)*2)) | (temp_l < (median(temp_l)/2)) |...
        (temp_angle > (median(temp_angle)+0.05*pi)) | (temp_angle < (median(temp_angle)-0.05*pi));
    
    l1(ind,:) = [];
    l2(ind,:) = [];
    matchedPoints1(ind) = [];
    matchedPoints2(ind) = [];
    
    l1(:,3) = 1;
    l2(:,3) = 1;
    l1 = l1';
    l2 = l2';

end