function figure_count = displaySurfaceAndLightPos(xC, yC, zC, LightPos, CamCalib, figure_count)


    temp_xyzC_tire_original_est(:,1) = reshape(xC, 1, numel(xC));
    temp_xyzC_tire_original_est(:,2) = reshape(yC, 1, numel(yC));
    temp_xyzC_tire_original_est(:,3) = reshape(zC, 1, numel(zC));


    pc_temp_xyzC_tire_original_est = pointCloud( [temp_xyzC_tire_original_est(:,1),temp_xyzC_tire_original_est(:,3),-temp_xyzC_tire_original_est(:,2)] );
    pc_temp_xyzC_tire_original_est_down_size = pcdownsample(pc_temp_xyzC_tire_original_est,'gridAverage',1);

    figure_count = figure_count + 1;
    figure_count = displayLightPosAndBalls( LightPos, [], [], CamCalib , figure_count , 1,0);
    figure(figure_count);
    hold on; pcshow(pc_temp_xyzC_tire_original_est_down_size);

end

