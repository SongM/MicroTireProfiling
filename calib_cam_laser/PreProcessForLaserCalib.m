function [bs_laser, error_laser_line] = PreProcessForLaserCalib(bs_cam_calib,b_display,redness_th,line_err_std_th,board_ind)
    tic;
    
    if (~exist('board_ind','var')),  board_ind = 0;  end
    B_DISPLAY = 0;
    
    fprintf([num2str(numel(bs_cam_calib)), ' boards are used here to detect laser points on checkerboards.\n']);
    fprintf(['  redness_th=', num2str(redness_th), '.\n']);
    fprintf(['  line_err_std_th=', num2str(line_err_std_th), '.\n']);
    
    board_proc = 1:numel(bs_cam_calib);
    bs_laser = [];
    error_laser_line = [];
    
    
    
    if (sum(board_ind)~=0),     board_proc = board_ind;    b_display = 1;    B_DISPLAY = 1; end

    for kk = board_proc
        fprintf(1,['[laser calib] board ', num2str(kk), '/', num2str(numel(board_proc)), ', ']);
        board = bs_cam_calib(kk);
        im = im2double(getImFromPath(board.filepath));
        p_cornerx = board.pixel_plane.corner(1,:);
        p_cornery = board.pixel_plane.corner(2,:);

        board_roi = roipoly(im , p_cornerx, p_cornery);
        board_im = repmat(board_roi,1,1,3).*im;
        board_im_r_capped = board_im(:,:,1) - max(board_im(:,:,2:3),[],3);
        board_im_r_capped = board_im_r_capped/max(board_im_r_capped(:));
        board_im_r_capped(board_im_r_capped<0) = 0;
        brightness_th = mean(board_im_r_capped(board_im_r_capped>redness_th));
        board_im_r_capped(board_im(:,:,1)>0.9999) = 1;
%         figure(); imshow(board_im_r_capped);
        r_x = [];
        r_y = [];
        if(B_DISPLAY);
            figure();displayBoard2D(board);hold on;
        end
            
        for i=1:size(board_im_r_capped,1)
            board_im_r_capped_line = board_im_r_capped(i,:);
%             brightness_th = mean(board_im_r_capped_line(board_im_r_capped_line>redness_th));
            bright_pixel = find(board_im_r_capped_line>brightness_th);
            if (numel(bright_pixel)>5)
                diff_bright_pixel = bright_pixel(2:end) - bright_pixel(1:end-1) - 1;
                if (sum(diff_bright_pixel) == 0)
                    r_y = [r_y;i];
                    r_x = [r_x;mean(bright_pixel)];
                    
                    if(B_DISPLAY);
                        scatter(bright_pixel,i*ones(1,numel(bright_pixel)),1,'g');
                    end
                    
                    
                end
            end
        end
    % line_eqn: x=ay+b;
        if (numel(r_x)<10)
            disp('no line on the board');
            if (b_display)
                figure();displayBoard2D(board);hold on;
            end
            continue;
        end
        [line_eqn, line_err] = linearFit([r_x';r_y']);    
        error_laser_line = [error_laser_line,line_err];
        
        fit_line_y = [min(r_y); max(r_y)];
        fit_line_x = (1-line_eqn(2)*fit_line_y)/line_eqn(1);
%         fit_line_x = [fit_line_y, ones(numel(fit_line_y),1)]*line_eqn;
        fit_line_length = sqrt( (fit_line_x(2)-fit_line_x(1))^2 +...
            (fit_line_y(2)-fit_line_y(1))^2);

        [board_x,board_y]=find(board_roi==1);
        board_x_range = max(board_x)-min(board_x);
        board_y_range = max(board_y)-min(board_y);

        if (b_display)
            figure();displayBoard2D(board);hold on;
            plot(fit_line_x,fit_line_y,'g','LineWidth',2);
            scatter(r_x,r_y,1,'b');
        end

        if (fit_line_length<0.6*min(board_x_range,board_y_range))
            disp('the laser line on the board is too short');
            continue;
        end
%         [line_err]
        if (line_err.std>line_err_std_th)
            fprintf(1, ['not use for laser calib, line_err_std=', ...
                num2str(line_err.std), '.\n']);
            continue;
        end    
        fprintf(1, ['line_err_std=', num2str(line_err.std),'.\n']);
        
        
        board.laserCalib.p_laser_points = [r_x';r_y'];
        board.laserCalib.p_laser_line = line_eqn;
        board.laserCalib.p_laser_line_end_points = [fit_line_x';fit_line_y'];
        bs_laser = [bs_laser,board];
    end
    
    fprintf([num2str(numel(bs_laser)), ' boards are used for laser calibration.\n']);
    toc;

    
    
end