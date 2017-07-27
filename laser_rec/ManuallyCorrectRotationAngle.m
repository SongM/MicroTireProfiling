function [laser_results, manualRotationCorrection] =...
    ManuallyCorrectRotationAngle(mat_file, laser_results_beforeCorrection, manualRotationCorrection,...
    b_display)

    d_thetas_beforeCorrection = [laser_results_beforeCorrection.d_theta];
    thetas_beforeCorrection = cumsum(d_thetas_beforeCorrection);

    laser_results = laser_results_beforeCorrection;

    for i = 1:size(manualRotationCorrection,1)
        manualRotationCorrection(i,4) = 0;
        ind1 = manualRotationCorrection(i,1);
        ind2 = manualRotationCorrection(i,2);

        manualRotationCorrection(i,5) = ...
            thetas_beforeCorrection(ind2) - thetas_beforeCorrection(ind1);

        manualRotationCorrection(i,6) = getAngelDifferenceBetweenTwoFrames(mat_file, ind1, ind2, 0) + manualRotationCorrection(i,3);
        correction_scale = manualRotationCorrection(i,6)/manualRotationCorrection(i,5);
        manualRotationCorrection(i,7) = correction_scale;

        for laser_i = (ind1+1):(ind2)
            laser_results(laser_i).d_theta = laser_results_beforeCorrection(laser_i).d_theta * correction_scale;
        end
    end

        ind1 = 1;
        ind2 = manualRotationCorrection(1,1);
        correction_scale = manualRotationCorrection(1,7);
        for laser_i = (ind1+1):(ind2)
            laser_results(laser_i).d_theta = laser_results_beforeCorrection(laser_i).d_theta * correction_scale;
        end
        ind1 = manualRotationCorrection(end,1);
        ind2 = numel(laser_results);
        correction_scale = manualRotationCorrection(end,7);
        for laser_i = (ind1+1):(ind2)
            laser_results(laser_i).d_theta = laser_results_beforeCorrection(laser_i).d_theta * correction_scale;
        end
        
        
    if (b_display)
        figure();
        % plot(cumsum([laser_results_forC0andV.d_theta]),'r');
        hold on;
        plot(thetas_beforeCorrection,'b');
        plot(cumsum([laser_results.d_theta]),'g');
        legend('beforeCorrection', 'afterCorrection');
        scatter(manualRotationCorrection(:,1),zeros(size(manualRotationCorrection,1),1),'r*');
        scatter(manualRotationCorrection(:,2),zeros(size(manualRotationCorrection,1),1),'r*');
        disp(manualRotationCorrection);
    end


end


