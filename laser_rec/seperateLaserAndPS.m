function [laser_file_list, ps_set_file_list] = seperateLaserAndPS(folder_path, laserandps_file_name)
    folder_name = folder_path;
    im_file_list = dir([folder_name, '\', laserandps_file_name, '\*JPG']);
%     im_file_list = dir([folder_name, '\', laserandps_file_name, '\*DNG']);
    N_im_laser_before_rotate = 5;
    N_im_laser_after_rotate = 5;
    N_im_extra = 1;
    N_im_ps_image = 16;
    N_im_per_set = N_im_laser_before_rotate + N_im_laser_after_rotate + N_im_extra + N_im_ps_image;
    
    if (mod(numel(im_file_list),N_im_per_set)~=0)
        fprintf(2,['some set is imcomplete, N_im_per_set=', num2str(N_im_per_set), ', N_im_total=', num2str(numel(im_file_list))]);
        return;
    end
    N_set = numel(im_file_list)/N_im_per_set;
    fprintf(['N_set=',num2str(N_set),...
        ', N_im_laser_before_rotate=',num2str(N_im_laser_before_rotate), ...
        ', N_im_laser_after_rotate=',num2str(N_im_laser_after_rotate), ...
        ', N_im_extra=',num2str(N_im_extra), ...
        ', N_im_ps_image=',num2str(N_im_ps_image), '.\n']);
    laser_file_list = [];
    for set_i = 1:N_set
        laser_file_list = [laser_file_list;...
            im_file_list((1:N_im_laser_before_rotate) + (set_i-1)*N_im_per_set);...
            im_file_list((1:N_im_laser_after_rotate) - N_im_laser_after_rotate + set_i*N_im_per_set)] ;
        ps_ind_file_list = [];
        ps_ind_file_list.laser_im_file = im_file_list( (set_i-1)*N_im_per_set + N_im_laser_before_rotate);
        ps_ind_file_list.extra_im_file = im_file_list( (set_i-1)*N_im_per_set + N_im_laser_before_rotate + 1);
        ps_ind_file_list.im_file_list = im_file_list( (1:N_im_ps_image) + (set_i-1)*N_im_per_set + N_im_laser_before_rotate + N_im_extra);
        ps_set_file_list(set_i) = ps_ind_file_list;
    end
    
    if(0)
        for i=1:16
            temp_im = getImFromPath([folder_path,'/',param.file_name.laserandps,'/',ps_set_file_list(1).im_file_list(i).name]);
            figure();imshow(temp_im/quantile(temp_im(:),1-1e-3));
            title(['quantile 1e-3 = ', num2str(quantile(temp_im(:),1-1e-3))]);
        end
        
        for i=1:numel(laser_file_list)
            temp_im = getImFromPath([folder_path,'/',param.file_name.laserandps,'/',laser_file_list(i).name]);
            figure();imshow(temp_im/quantile(temp_im(:),1-1e-2));
            title(['quantile 1e-2 = ', num2str(quantile(temp_im(:),1-1e-2))]);
        end
    end
    
end