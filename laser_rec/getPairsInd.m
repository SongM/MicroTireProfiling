function indexPairs = getPairsInd(lr1, lr2, CamCalib)

    features1 = lr1.features;
    features2 = lr2.features;
    validPoints1 = lr1.validPoints;
    validPoints2 = lr2.validPoints;

    indexPairs = matchFeatures(features1,features2);

    matchedPoints1 = validPoints1(indexPairs(:,1),:);
    matchedPoints2 = validPoints2(indexPairs(:,2),:);

    % 
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

end
