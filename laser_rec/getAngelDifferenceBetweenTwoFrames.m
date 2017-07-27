function angleDifference = getAngelDifferenceBetweenTwoFrames(mat_file, ind1, ind2, b_display)
    
    load(mat_file);
    load(individual_mat_file.laserrec,'laser_results');

    lr1 = laser_results(ind1);
    lr2 = laser_results(ind2);
    
    indexPairs = getPairsInd(lr1, lr2, CamCalib);
    
    matchedPoints1 = lr1.validPoints(indexPairs(:,1),:);
    matchedPoints2 = lr2.validPoints(indexPairs(:,2),:);

    Laser_X_c = lr1.Laser_X_c;

    l1 = normalize_pixel(matchedPoints1.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
    l2 = normalize_pixel(matchedPoints2.Location,CamCalib.fc,CamCalib.cc,CamCalib.kc,CamCalib.alpha_c);
    l1(:,3) = 1;
    l2(:,3) = 1;
    l1 = l1';
    l2 = l2';
    
    Z1_est = mean(Laser_X_c(3,:)).*ones(1,size(l1,2));
    Z2_est = Z1_est;
        
    [angleDifference, ~, ~, ~, ~, ~, ~, ~] =...
        getRotationAxisAndAngle(C0, V, Z1_est, Z2_est,...
        Laser_X_c, l1, l2, 1, 0, figure_count);
    
    fprintf(['laser_',num2str(ind1),'->laser_', num2str(ind2),...
        ', N of matched pairs = ', num2str(size(indexPairs,1)),...
        ', angle_difference = ', num2str(angleDifference), '.\n']);
    if(b_display)
        im1 = getImFromPath(lr1.im_path);
        im1 = im1/quantile(im1(:), 0.99);
        im2 = getImFromPath(lr2.im_path);
        im2 = im2/quantile(im2(:), 0.99);

        figure();
        showMatchedFeatures(im1,im2,matchedPoints1,matchedPoints2);
        title([num2str(ind1), '->', num2str(ind2)]);

        figure();
        subplot(1,2,1);imshow(im1);title(num2str(ind1));
        subplot(1,2,2);imshow(im2);title(num2str(ind2));
    end


end