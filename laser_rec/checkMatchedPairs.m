function checkMatchedPairs(mat_file, laser_i)
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');

    sigma = param.laserrec.sigma;
    dist_th = param.laserrec.dist_th;
    tire_region = param.tire_region;
    laser_region = param.laserrec.laser_region;
    N_strongest = param.laserrec.N_strongest;

    filter = fspecial('Gaussian',4*sigma+1,sigma);

    
    lr1 = laser_results(laser_i-1);
    lr2 = laser_results(laser_i);

    im1 = getImFromPath(lr1.im_path);    
    im2 = getImFromPath(lr2.im_path);

    indexPairs = laser_results(laser_i).indexPairs;
    matchedPoints1 = lr2.matchedPoints1;
    matchedPoints2 = lr2.matchedPoints2;
    
    SAT_PERCENT1 = lr1.SAT_PERCENT;
    SAT_PERCENT2 = lr2.SAT_PERCENT;
    
    [filtered_im1, sat_val1] = getFilteredIm(im1, filter, SAT_PERCENT1);
    [filtered_im2, sat_val2] = getFilteredIm(im2, filter, SAT_PERCENT2);
    
    
    figure_count = figure_count + 1;    figure(figure_count); 
    showMatchedFeatures(filtered_im1,filtered_im2,matchedPoints1,matchedPoints2);
    title(['laser_',num2str(laser_i)]);
    figure_count = figure_count + 1;    figure(figure_count); 
    subplot(1,2,1); imshow(filtered_im1);
    subplot(1,2,2); imshow(filtered_im2);
    
    
    prompt = {'SAT_PERCENT1','SAT_PERCENT2'};
    dlg_title = 'Input';

    answer = 0;
    b_display = 1;
    while(numel(answer)>0)
        defaultans = {num2str(SAT_PERCENT1), num2str(SAT_PERCENT2)};
        answer = inputdlg(prompt,dlg_title,1,defaultans);
        if (numel(answer)>0)
            SAT_PERCENT1 = str2num(answer{1});
            SAT_PERCENT2 = str2num(answer{2});

            figure_count = figure_count + 1;
            [filtered_im1, ~, features1, validPoints1, sat_val1]...
                = getStrongestFeatureFromIm(im1, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, b_display, figure_count, SAT_PERCENT1);
            title(['sat_val=',num2str(sat_val1)]); drawnow();
            figure_count = figure_count + 1;
            [filtered_im2, ~, features2, validPoints2, sat_val2]...
                = getStrongestFeatureFromIm(im2, filter, N_strongest, dist_th, LaserCalib, CamCalib, laser_region, tire_region, b_display, figure_count, SAT_PERCENT2);
            title(['sat_val=',num2str(sat_val2)]); drawnow();

            indexPairs = matchFeatures(features1,features2);
            
            matchedPoints1 = validPoints1(indexPairs(:,1),:);
            matchedPoints2 = validPoints2(indexPairs(:,2),:);
            l1 = normalize_pixel(matchedPoints1.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
            l2 = normalize_pixel(matchedPoints2.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);

            temp = l1-l2;

            temp_l = sqrt(sum(temp.^2,2));
            temp_angle = atan2(temp(:,1),temp(:,2));
            ind = (temp_l > (median(temp_l)*2)) | (temp_l < (quantile(temp_l,0.25)/2)) |...
                (temp_angle > (median(temp_angle)+0.05*pi)) | (temp_angle < (median(temp_angle)-0.05*pi));
            
            while (std(temp_l(~ind))/mean(temp_l(~ind))>0.1)
                ind = ind | (temp_l>mean(temp_l(~ind)));
            end
            
            indexPairs(ind,:) = [];
            
            matchedPoints1 = validPoints1(indexPairs(:,1),:);
            matchedPoints2 = validPoints2(indexPairs(:,2),:);
            
            figure_count = figure_count + 1;    figure(figure_count); 
            showMatchedFeatures(filtered_im1,filtered_im2,matchedPoints1,matchedPoints2);
            title(['laser_',num2str(laser_i),', NofMatchedPairs = ', num2str(size(indexPairs,1))]);
       
        
        end
    end
%     laser_results(laser_i).indexPairs = indexPairs;

    laser_results(laser_i).matchedPoints1 = matchedPoints1;
    laser_results(laser_i).matchedPoints2 = matchedPoints2;
    fprintf(['laser_',num2str(laser_i),': N_pairs = ',...
        num2str(size(laser_results(laser_i).indexPairs,1)),...
        ', sat_val = ', num2str(laser_results(laser_i).sat_val), '.\n']);

    clear laser_i sigma dist_th tire_region laser_region filter N_strongest
    clear lr1 lr2 im1 im2 matchedPoints1 matchedPoints2 
    clear SAT_PERCENT1 SAT_PERCENT2 filtered_im1 filtered_im2 sat_val1 sat_val2
    clear prompt dlg_title answer b_display defaultans
    clear features1 validPoints1 features2 validPoints2 
    clear l1 l2 temp temp_l temp_angle ind indexPairs

    
    save(individual_mat_file.laserrec,'-v7.3');
    clear laser_results;
    saveMatFile;
end

function [filtered_im, sat_val] = getFilteredIm(im, filter, SAT_PERCENT)
    sat_val = quantile(im(:),1-SAT_PERCENT);
    temp_im = im2double(rgb2gray(im)) / sat_val;
    temp_im(temp_im>1) = 1;
    filtered_im = imfilter(temp_im,filter);
end

 