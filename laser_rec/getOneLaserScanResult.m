function [p_laser_x,Laser_X_c] = getOneLaserScanResult(im,LaserCalib,CamCalib,laser_region,redness_th,b_display)

    corner_minx = laser_region(1);
    corner_maxx = laser_region(3);
    corner_miny = laser_region(2);
    corner_maxy = laser_region(4);
    
    cornerx = [corner_minx,corner_minx,corner_maxx,corner_maxx];
    cornery = [corner_miny,corner_maxy,corner_maxy,corner_miny];
    trim_roi = roipoly(im,cornerx, cornery);
    
    trim_im = repmat(trim_roi,1,1,3).*im;
    trim_im_r_capped = trim_im(:,:,1)- max(trim_im(:,:,2:3),[],3);
    trim_im_r_capped = trim_im_r_capped/max(trim_im_r_capped(:));
    trim_im_r_capped(trim_im_r_capped<0) = 0;
    trim_im_r_capped(trim_im(:,:,1)>0.9999) = 1;

%     figure();
% imagesc(trim_im_r_capped);
% colorbar
% colormap colorcube
% colormap jet
    
    
    
    r_x = [];   r_y = [];
    for i=1:size(trim_im_r_capped,1)
        trim_im_r_capped_line = trim_im_r_capped(i,:);
        brightness_th = mean(trim_im_r_capped_line(trim_im_r_capped_line>redness_th));
        bright_pixel = find(trim_im_r_capped_line>brightness_th);
        if (numel(bright_pixel)>5)
            diff_bright_pixel = bright_pixel(2:end) - bright_pixel(1:end-1) - 1;
            if (sum(diff_bright_pixel) < 5)
                r_y = [r_y;i];
%                 r_x = [r_x;mean(bright_pixel)];
                r_x = [r_x;mean(bright_pixel.*trim_im_r_capped_line(bright_pixel))/mean(trim_im_r_capped_line(bright_pixel))];
%                 fprintf(['r_y=',num2str(i),', diff=',num2str(mean(bright_pixel)-mean(bright_pixel.*trim_im_r_capped_line(bright_pixel))/mean(trim_im_r_capped_line(bright_pixel))),'.\n']);
            end
        end
    end
    p_laser_x = [r_x';r_y'];
    Laser_X_c = get3DCoorFromPlaneEqnAnd2DCoor(LaserCalib.c_plane_eqn, p_laser_x, CamCalib);

    
    if(b_display)
        figure();
        subplot(1,2,1); hold on; axis equal; view([0,-90]); axis([-200,100,-200,200,500,1500])
        scatter3(Laser_X_c(1,:),Laser_X_c(2,:),Laser_X_c(3,:),1);
        subplot(1,2,2);imshow(im);
        hold on;plot([cornerx,cornerx(1)],[cornery,cornery(1)],'b');

        drawnow;
    end
end