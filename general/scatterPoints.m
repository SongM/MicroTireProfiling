function scatterPoints(points,scatter_size,scatter_color,b_filled)
    

    if (~exist('scatter_size','var')),  scatter_size = 1;  end
    if (~exist('scatter_color','var')),  scatter_color = [1,0,0];  end
    if (~exist('b_filled','var')),  b_filled = 1;  end

    
    [nx,ny] = size(points);
    
    if b_filled
        if (nx == 2)
            scatter(points(1,:),points(2,:),scatter_size,scatter_color,'filled');
            return
        end
        if (nx == 3)
            scatter3(points(1,:),points(2,:),points(3,:),scatter_size,scatter_color,'filled');
            return
        end
        if (ny == 2)
            scatter(points(:,1),points(:,2),scatter_size,scatter_color,'filled');
            return
        end
        if (ny == 3)
            scatter3(points(:,1),points(:,2),points(:,3),scatter_size,scatter_color,'filled');
            return
        end
    else
        if (nx == 2)
            scatter(points(1,:),points(2,:),scatter_size,scatter_color);
            return
        end
        if (nx == 3)
            scatter3(points(1,:),points(2,:),points(3,:),scatter_size,scatter_color);
            return
        end
        if (ny == 2)
            scatter(points(:,1),points(:,2),scatter_size,scatter_color);
            return
        end
        if (ny == 3)
            scatter3(points(:,1),points(:,2),points(:,3),scatter_size,scatter_color);
            return
        end
    end
end

