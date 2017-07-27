function compareFinalRes(final_res_after, final_res_before, tire_after, tire_before)
    fr_a_P_Tire = final_res_after.grid_P_Tire;
    fr_b_P_Tire = final_res_before.grid_P_Tire;
    
    fr_a_Im = final_res_after.grid_rerenderIm;
    fr_b_Im = final_res_before.grid_rerenderIm;
    
    fr_a_color_vec = final_res_after.grid_color_vec_map;
    fr_b_color_vec = final_res_before.grid_color_vec_map;
    figure();
    subplot(2,1,1),imshow(fr_b_Im/quantile(fr_b_Im(:),0.99));
    subplot(2,1,2),imshow(fr_a_Im/quantile(fr_b_Im(:),0.99));
    
    
%     temp_fr_b_Im = fr_b_Im(10:78,1:500,1);
%     temp_fr_a_Im = fr_a_Im(10:78,1:500,1);
%     
%     
%     temp_fr_b_Im_BW = edge(temp_fr_b_Im,'Sobel',0.02,'vertical');
%     temp_fr_a_Im_BW = edge(temp_fr_a_Im,'Sobel',0.02,'vertical');
%     figure();
%     subplot(2,2,1);imshow(temp_fr_b_Im_BW)
%     subplot(2,2,2);imshow(temp_fr_a_Im_BW)
%     subplot(2,2,3);imshow(temp_fr_b_Im/quantile(temp_fr_b_Im(:),0.99));
%     subplot(2,2,4);imshow(temp_fr_a_Im/quantile(temp_fr_a_Im(:),0.99));
    
    shift = 7000-100;
    fr_a_color_vec = fr_a_color_vec(:,[(shift+1):end,1:shift],:);
    fr_a_P_Tire = fr_a_P_Tire(:,[(shift+1):end,1:shift],:);
    fr_a_Im = fr_a_Im(:,[(shift+1):end,1:shift],:);
    figure();
    subplot(2,1,1),imshow(fr_b_Im/quantile(fr_b_Im(:),0.99));
    subplot(2,1,2),imshow(fr_a_Im/quantile(fr_b_Im(:),0.99));
    
    figure();imshow(fr_b_Im(:,3486:4386,:)/quantile(fr_b_Im(:),0.99))
    figure();imshow(fr_a_Im(:,573:1473,:)/quantile(fr_b_Im(:),0.99))
    
%     t_fr_a_P_Tire = fr_a_P_Tire;
%     t_fr_b_P_Tire = fr_b_P_Tire;
  %%
  
    
    b_groove_ref_R = tire_before.laser.groove_ref_R;
    b_max_tire_R = max([tire_before.laser.Tire_slices.R])+2;
    b_groove_Z_range = [min(tire_before.laser.groove_pixel_Tire_all{3}(3,:)),...
    max(tire_before.laser.groove_pixel_Tire_all{3}(3,:))];

    fr_b_depth = zeros(size(fr_b_P_Tire,1), size(fr_b_P_Tire,2));
    
    for i=1:size(fr_b_P_Tire,2)
        t_fr_b_P_Tire_i = squeeze(fr_b_P_Tire(:,i,:));
        ind = find( (t_fr_b_P_Tire_i(:,3)>b_groove_Z_range(1))&...
            (t_fr_b_P_Tire_i(:,3)<b_groove_Z_range(2))&...
            (t_fr_b_P_Tire_i(:,1)>0) );
        
        if (numel(ind)>0)
            t_fr_b_P_Tire_i_groove = t_fr_b_P_Tire_i(ind,:);
            [min_groove_R,min_groove_ind] = min(t_fr_b_P_Tire_i_groove(:,1));

            groove_ind = ind(min_groove_ind);
            scaler = min_groove_R/b_groove_ref_R;

            t_fr_b_depth_i = t_fr_b_P_Tire_i(:,1)/scaler;
            t_fr_b_depth_i(t_fr_b_depth_i>b_max_tire_R) = 0;
            
            t_fr_b_depth_i = t_fr_b_depth_i - b_groove_ref_R;
            t_fr_b_depth_i(t_fr_b_depth_i<0) = 0;
            fr_b_depth(:,i) = t_fr_b_depth_i;
            fprintf(['i=',num2str(i),', groove_ind=', num2str(groove_ind),...
                ', min_groove_R=',num2str(min_groove_R),...
                ', scaler=', num2str(scaler), '.\n']);

        end
    end
    figure();imagesc(fr_b_depth);axis equal;colormap jet;
    caxis([0,b_max_tire_R-b_groove_ref_R]);colorbar;
        title('depth before');

    

%%
    
    a_groove_ref_R = tire_after.laser.groove_ref_R;
    a_max_tire_R = max([tire_after.laser.Tire_slices.R])+2;
    a_groove_Z_range = [min(tire_after.laser.groove_pixel_Tire_all{3}(3,:)),...
    max(tire_after.laser.groove_pixel_Tire_all{3}(3,:))];

    fr_a_depth = zeros(size(fr_a_P_Tire,1), size(fr_a_P_Tire,2));
    
    for i=1:size(fr_a_P_Tire,2)
        t_fr_a_P_Tire_i = squeeze(fr_a_P_Tire(:,i,:));
        ind = find( (t_fr_a_P_Tire_i(:,3)>a_groove_Z_range(1))&...
            (t_fr_a_P_Tire_i(:,3)<a_groove_Z_range(2))&...
            (t_fr_a_P_Tire_i(:,1)>0) );
        
        if (numel(ind)>0)
            t_fr_a_P_Tire_i_groove = t_fr_a_P_Tire_i(ind,:);
            [min_groove_R,min_groove_ind] = min(t_fr_a_P_Tire_i_groove(:,1));

            groove_ind = ind(min_groove_ind);
            scaler = min_groove_R/a_groove_ref_R;

            t_fr_a_depth_i = t_fr_a_P_Tire_i(:,1)/scaler;
            t_fr_a_depth_i(t_fr_a_depth_i>a_max_tire_R) = 0;
            
            t_fr_a_depth_i = t_fr_a_depth_i - a_groove_ref_R;
            t_fr_a_depth_i(t_fr_a_depth_i<0) = 0;
            fr_a_depth(:,i) = t_fr_a_depth_i;
            fprintf(['i=',num2str(i),', groove_ind=', num2str(groove_ind),...
                ', min_groove_R=',num2str(min_groove_R),...
                ', scaler=', num2str(scaler), '.\n']);

        end
    end
    figure();imagesc(fr_a_depth);axis equal;colormap jet;
    caxis([0,a_max_tire_R-a_groove_ref_R]);colorbar;
    title('depth after');
    
    
    figure();imagesc(fr_b_depth(:,3486:4386));axis equal;colormap jet;
    caxis([0,b_max_tire_R-b_groove_ref_R]);colorbar;
    title('depth before');

    figure();imagesc(fr_a_depth(:,573:1473));axis equal;colormap jet;
    caxis([0,a_max_tire_R-a_groove_ref_R]);colorbar;
    title('depth after');
    
    depth_diff = fr_b_depth - fr_a_depth;
    ind = (fr_b_depth>1) & (fr_a_depth>1);
    depth_diff(~ind) = 0;
    figure();imagesc(depth_diff);axis equal;colormap jet;
    
    t_std = [];
    t_N_minus = [];
    for shift_i = 2:size(fr_b_depth,2)
        t_fr_a_depth = fr_a_depth(:,[shift_i:end,1:(shift_i-1)]);
        depth_diff = fr_b_depth - t_fr_a_depth;
        ind = (fr_b_depth>3) & (t_fr_a_depth>3);
        depth_diff(~ind) = 0;
        t_std = [t_std,std(depth_diff(ind))];
        t_N_minus = [t_N_minus, sum(depth_diff(ind)<0)];
    end
    
    
    shift_i = 70;
        t_fr_a_depth = fr_a_depth(:,[shift_i:end,1:(shift_i-1)]);
    depth_diff = fr_b_depth - t_fr_a_depth;
    ind = (fr_b_depth>3) & (t_fr_a_depth>3);
    depth_diff(~ind) = 0;

    figure();imagesc(depth_diff);axis equal;colormap jet;
    
    t_diff=[];
    for i=1:3000
        if isempty(tire_before.laser.Tire_slices(i).R)
            continue;
        end
        if isempty(tire_after.laser.Tire_slices(i).R)
            continue;
        end
        t_diff(i) = tire_before.laser.Tire_slices(i).R - tire_after.laser.Tire_slices(i).R - b_groove_ref_R + a_groove_ref_R;
    end
    
    
%         t_std = [t_std,std(depth_diff(ind))];
%         t_N_minus = [t_N_minus, sum(depth_diff(ind)<0)];

end

