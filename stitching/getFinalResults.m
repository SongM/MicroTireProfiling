function final_res = getFinalResults(tire)
    
    
    P_Tire_all = [];
    P_E_after_theta_shift_all = [];
    color_vec_map_all = [];
    rerenderIm_all = [];

    P_Tire_stretched_all = [];
    
    for ps_i = 1:20
    
        tire_ps = tire.ps_set(ps_i);

        P_Tire = tire_ps.P_Tire;

        P_E_after_theta_shift = tire_ps.P_E_after_theta_shift;
        color_vec_map = reshape(tire_ps.color_vec_map,numel(tire_ps.color_vec_map)/3,3)';
        rerenderIm = reshape(tire_ps.rerenderIm,numel(tire_ps.rerenderIm)/3,3)';
        P_Tire_stretched = tire_ps.P_Tire_stretched;

        P_Tire_all = [P_Tire_all, P_Tire];

        P_E_after_theta_shift_all = [P_E_after_theta_shift_all, P_E_after_theta_shift];
        color_vec_map_all = [color_vec_map_all, color_vec_map];
        rerenderIm_all = [rerenderIm_all, rerenderIm];
        P_Tire_stretched_all = [P_Tire_stretched_all, P_Tire_stretched];
    end
    
    grid_x_min = min(P_Tire_stretched_all(1,:));
    grid_x_max = max(P_Tire_stretched_all(1,:));

    grid_y_min = min(P_Tire_stretched_all(2,:));
    grid_y_max = max(P_Tire_stretched_all(2,:));

    grid_y_NUM = size(tire_ps.color_vec_map,1);
    
    grid_y = linspace(grid_y_min, grid_y_max, grid_y_NUM);
    grid_size = grid_y(2)-grid_y(1);

    grid_x_NUM = round((grid_x_max-grid_x_min)/grid_size)+1;
    grid_x_NUM = 7237;
    
    
    grid_x = linspace(grid_x_min, grid_x_max, grid_x_NUM);
    
%     [grid_xx, grid_yy] = meshgrid(grid_x, grid_y);

    
    
    
    temp_p = [];
    temp_p(1,:) = P_Tire_stretched_all(1,:) - grid_x_min;
    temp_p(2,:) = P_Tire_stretched_all(2,:) - grid_y_min;
    round_temp_p = round(temp_p/grid_size) + 1;
    round_error = sum( (round_temp_p - temp_p/grid_size-1).^2);
    


    
    grid_map_ind = zeros(grid_y_NUM,grid_x_NUM);
    error_map = 2*ones(grid_y_NUM,grid_x_NUM);

    for p_i = 1:size(P_Tire_stretched_all,2)
        jj = round_temp_p(1,p_i);
        ii = round_temp_p(2,p_i);
        
        if (error_map(ii,jj)>round_error(p_i))
            grid_map_ind(ii,jj) = p_i;
            error_map(ii,jj) = round_error(p_i);
        end
    end
    

% 
%     grid_X = zeros(grid_y_NUM,grid_x_NUM);
%     grid_X(grid_map_ind>0) = P_Tire_stretched_all(1,grid_map_ind(grid_map_ind>0));
%    
%     grid_Y = zeros(grid_y_NUM,grid_x_NUM);
%     grid_Y(grid_map_ind>0) = P_Tire_stretched_all(2,grid_map_ind(grid_map_ind>0));
    
%     P_E_after_theta_shift_all = [];


    grid_P_Tire = zeros(grid_y_NUM,grid_x_NUM,3);
    for i=1:3
        temp_grid = zeros(grid_y_NUM,grid_x_NUM);
        temp_grid(grid_map_ind>0) = P_Tire_all(i,grid_map_ind(grid_map_ind>0));
        grid_P_Tire(:,:,i) = temp_grid;
    end
    
    grid_P_E_after_theta_shift = zeros(grid_y_NUM,grid_x_NUM,3);
    for i=1:3
        temp_grid = zeros(grid_y_NUM,grid_x_NUM);
        temp_grid(grid_map_ind>0) = P_E_after_theta_shift_all(i,grid_map_ind(grid_map_ind>0));
        grid_P_E_after_theta_shift(:,:,i) = temp_grid;
    end
    
    grid_color_vec_map = zeros(grid_y_NUM,grid_x_NUM,3);
    for i=1:3
        temp_grid = zeros(grid_y_NUM,grid_x_NUM);
        temp_grid(grid_map_ind>0) = color_vec_map_all(i,grid_map_ind(grid_map_ind>0));
        grid_color_vec_map(:,:,i) = temp_grid;
    end
    

    grid_rerenderIm = zeros(grid_y_NUM,grid_x_NUM,3);
    for i=1:3
        temp_grid = zeros(grid_y_NUM,grid_x_NUM);
        temp_grid(grid_map_ind>0) = rerenderIm_all(i,grid_map_ind(grid_map_ind>0));
        grid_rerenderIm(:,:,i) = temp_grid;
    end
    
    grid_P_Tire_stretched = zeros(grid_y_NUM,grid_x_NUM,3);
    for i=1:3
        temp_grid = zeros(grid_y_NUM,grid_x_NUM);
        temp_grid(grid_map_ind>0) = P_Tire_stretched_all(i,grid_map_ind(grid_map_ind>0));
        grid_P_Tire_stretched(:,:,i) = temp_grid;
    end
    
    grid_Z = zeros(grid_y_NUM,grid_x_NUM);
    grid_Z(grid_map_ind>0) = P_Tire_stretched_all(3,grid_map_ind(grid_map_ind>0));
    figure();imagesc(grid_Z);axis equal;colorbar;colormap colorcube;
    caxis([tire.laser.groove_ref_R-1, max([tire.laser.Tire_slices.R]+2)]);
    figure();imshow(grid_color_vec_map);
    figure();imshow(grid_rerenderIm/quantile(grid_rerenderIm(:),0.99))

    
    
    final_res.grid_P_Tire = grid_P_Tire;
    final_res.grid_P_E_after_theta_shift = grid_P_E_after_theta_shift;
    final_res.grid_color_vec_map = grid_color_vec_map;
    final_res.grid_rerenderIm = grid_rerenderIm;
    final_res.grid_P_Tire_stretched = grid_P_Tire_stretched;
    
    
end