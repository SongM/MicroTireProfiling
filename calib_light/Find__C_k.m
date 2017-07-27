function [b_ck, ball_area, center_trimmed, radius, figure_count] =...
    Find__C_k(temp_im, resize_scale, ball_size_range, ball_range_scale,...
    b_display, figure_count)
    %% 
    % to reduce the calculation, first resize the image to find the area
    % where the ball exists, then find the ball in the orginal-size image,
    % but got rid of the rest of the image besides the ball area.
    b_debug = 0;
    temp_im_resize_r = imresize(temp_im(:,:,1), resize_scale);
    
    sensitivity = 0.9;
    [center_scaled, radius_scaled] =...
        imfindcircles(temp_im_resize_r, ball_size_range*resize_scale, ...
        'ObjectPolarity','dark', 'Sensitivity',sensitivity);
    while (numel(radius_scaled)==0)
        sensitivity = sensitivity + (1-sensitivity)*0.1;
        [center_scaled, radius_scaled] =...
            imfindcircles(temp_im_resize_r, ball_size_range*resize_scale, ...
            'ObjectPolarity','dark', 'Sensitivity',sensitivity);
    end

    
    
    
    
%     [center_scaled, radius_scaled] = imfindcircles(temp_im_resize_r,...
%         ball_size_range*resize_scale, 'ObjectPolarity','dark', 'Sensitivity',0.9);
%     if (numel(radius_scaled)==0)
%         fprintf(2, ['no circle is detected, N_circle = \n']);
%     end
    
    if (b_debug)
        if (numel(radius_scaled)>1)
            fprintf(1, ['more than 1 circles are detected, N_circle = ',...
                num2str(numel(radius_scaled)), '\n']);
        end
    end
    if(b_display)
        figure_count = figure_count + 1;
        figure(figure_count); 
        subplot(1,2,1);
        imshow(temp_im_resize_r);
        if (b_debug)
            viscircles(center_scaled, radius_scaled,'Color','b');
        end
        viscircles(center_scaled(1,:), radius_scaled(1,:),'Color','r');
        title(['(',num2str((center_scaled(1,:)-1)/resize_scale+1),...
            '), r=',num2str(radius_scaled(1,:)/resize_scale)]);
    end
    center_scaled = center_scaled(1,:);
    radius_scaled = radius_scaled(1,:);
    
    ball_area = [ floor( center_scaled-radius_scaled*ball_range_scale ), ...
        repmat( ceil(radius_scaled*ball_range_scale*2),1,2 ) ]/resize_scale;
%%
    
    temp_im_trimmed = temp_im( ball_area(2):(ball_area(2)+ball_area(4)),...
        ball_area(1):(ball_area(1)+ball_area(3)));
    
%     filter_N = ceil(mean(ball_area(3:4))/10);
%     temp_im_trimmed_filtered = imfilter(temp_im_trimmed, fspecial('gaussian',filter_N,filter_N/2));
%     temp_im_trimmed_filtered = temp_im_trimmed_filtered/max(temp_im_trimmed_filtered(:));
%     
%     [center_trimmed, radius] =...
%         imfindcircles(t_temp_im_trimmed, ball_size_range, ...
%         'ObjectPolarity','dark', 'Sensitivity',0.9,'EdgeThreshold',0);

    t_temp_im_trimmed = temp_im_trimmed/max(temp_im_trimmed(:))*3;
    t_temp_im_trimmed(t_temp_im_trimmed>1) = 1;
    
    sensitivity = 0.9;
    [center_trimmed, radius] =...
        imfindcircles(t_temp_im_trimmed, ball_size_range, ...
        'ObjectPolarity','dark', 'Sensitivity',sensitivity);
    while (numel(radius)==0)
        sensitivity = sensitivity + (1-sensitivity)*0.1;
        [center_trimmed, radius] =...
            imfindcircles(t_temp_im_trimmed, ball_size_range, ...
            'ObjectPolarity','dark', 'Sensitivity',sensitivity);
    end
    
%     if (numel(radius)==0)
%         fprintf(2, ['no circle is detected, N_circle = \n']);
%     end
    
    if (b_debug)
        if (numel(radius)>1)
            fprintf(1, ['more than 1 circles are detected, N_circle = ',...
                num2str(numel(radius)), '\n']);
        end
    end
    
    if(b_display)
        subplot(1,2,2);imshow(t_temp_im_trimmed);
        if (b_debug)
            viscircles(center_trimmed, radius,'Color','b');
        end
        viscircles(center_trimmed(1,:), radius(1,:),'Color','r');
        title({['sensitivity = ', num2str(sensitivity)], ...
            ['(',num2str(center_trimmed(1,:)+ball_area(1:2)-1),...
            '), r=',num2str(radius(1,:))]});
    end    
    
    b_ck = center_trimmed(1,:)+ball_area(1:2)-1;
    radius = radius(1,:);
    center_trimmed = center_trimmed(1,:);
end
