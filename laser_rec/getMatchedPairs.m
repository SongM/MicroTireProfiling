function getMatchedPairs(mat_file)
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');
    laser_i = 1;
    laser_results(laser_i).matched=0;
    laser_results(laser_i).indexPairs=[];
    laser_results(laser_i).matchedPoints1=[];
    laser_results(laser_i).matchedPoints2=[];
    fprintf(['laser_',num2str(laser_i),': N_pairs = ',...
        num2str(size(laser_results(laser_i).indexPairs,1)),...
        ', sat_val = ', num2str(laser_results(laser_i).sat_val), '.\n']);
    N_laser = numel(laser_results);
    for laser_i = 2:N_laser
        lr1 = laser_results(laser_i-1);
        lr2 = laser_results(laser_i);
        indexPairs = getPairsInd(lr1, lr2, CamCalib);
%         laser_results(laser_i).indexPairs = indexPairs;
        
        laser_results(laser_i).matchedPoints1 = lr1.validPoints(indexPairs(:,1),:);
        laser_results(laser_i).matchedPoints2 = lr2.validPoints(indexPairs(:,2),:);

        fprintf(['laser_',num2str(laser_i),': N_pairs = ',...
        num2str(size(laser_results(laser_i).matchedPoints1,1)),...
        ', sat_val = ', num2str(laser_results(laser_i).sat_val), '.\n']);
        clear lr1 lr2 indexPairs
    end
    clear laser_i N_laser
    save(individual_mat_file.laserrec,'-v7.3');
    clear laser_results;
    saveMatFile;
end

% 
% function indexPairs = getPairsInd(lr1, lr2, CamCalib)
% 
%     features1 = lr1.features;
%     features2 = lr2.features;
%     validPoints1 = lr1.validPoints;
%     validPoints2 = lr2.validPoints;
% 
%     indexPairs = matchFeatures(features1,features2);
% 
%     matchedPoints1 = validPoints1(indexPairs(:,1),:);
%     matchedPoints2 = validPoints2(indexPairs(:,2),:);
% 
%     % 
%     l1 = normalize_pixel(matchedPoints1.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
%     l2 = normalize_pixel(matchedPoints2.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
%      
%     temp = l1-l2;
%     
%     temp_l = sqrt(sum(temp.^2,2));
%     temp_angle = atan2(temp(:,1),temp(:,2));
%     ind = (temp_l > (median(temp_l)*2)) | (temp_l < (median(temp_l)/2)) |...
%         (temp_angle > (median(temp_angle)+0.05*pi)) | (temp_angle < (median(temp_angle)-0.05*pi));
%     
%     indexPairs(ind,:) = [];
% 
% end
