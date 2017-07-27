function [x,X] = getBoardCorners(board, I, b_display,...
    wintx, winty, n_sq_x, n_sq_y, dX, dY, map)
    % 4 extreme corner
    S_xxx = board.coor_x([1,end-n_sq_x,end,n_sq_x+1]);
    S_yyy = board.coor_y([1,end-n_sq_x,end,n_sq_x+1]);
    [Xc,good,bad,type] = cornerfinder([S_xxx;S_yyy],I,winty,wintx); % the four corners
    x = Xc(1,:)';
    y = Xc(2,:)';

    ind = [2 3 4 1];
    x = x(ind);
    y = y(ind);

    %%%%%%%% use homography to find the all the corners useing the 4 extreme
    %%%%%%%% corners.
    % Compute the inside points through computation of the planar homography (collineation)
    a00 = [x(1);y(1);1];
    a10 = [x(2);y(2);1];
    a11 = [x(3);y(3);1];
    a01 = [x(4);y(4);1];
    % Compute the planar collineation: (return the normalization matrix as well)
    [Homo,Hnorm,inv_Hnorm] = compute_homography([a00 a10 a11 a01],[0 1 1 0;0 0 1 1;1 1 1 1]);
    % Build the grid using the planar collineation:
    x_l = ((0:n_sq_x)'*ones(1,n_sq_y+1))/n_sq_x;
    y_l = (ones(n_sq_x+1,1)*(0:n_sq_y))/n_sq_y;
    pts = [x_l(:) y_l(:) ones((n_sq_x+1)*(n_sq_y+1),1)]';
    XX = Homo*pts;
    XX = XX(1:2,:) ./ (ones(2,1)*XX(3,:));
    % Complete size of the rectangle
    W = n_sq_x*dX;
    L = n_sq_y*dY;

    %%%%%%%%%%%%%%%%%%%%%%%% ADDITIONAL STUFF IN THE CASE OF HIGHLY DISTORTED IMAGES %%%%%%%%%%%%%
    quest_distort = 0;

    %%%%%%%%%%%%%%%%%%%%% END ADDITIONAL STUFF IN THE CASE OF HIGHLY DISTORTED IMAGES %%%%%%%%%%%%%
    Np = (n_sq_x+1)*(n_sq_y+1);

    grid_pts = cornerfinder(XX,I,winty,wintx); %%% Finds the exact corners at every points!
    %save all_corners x y grid_pts
    grid_pts = grid_pts - 1; % subtract 1 to bring the origin to (0,0) instead of (1,1) in matlab (not necessary in C)



    Xi = reshape(([0:n_sq_x]*dX)'*ones(1,n_sq_y+1),Np,1)';
    Yi = reshape(ones(n_sq_x+1,1)*[n_sq_y:-1:0]*dY,Np,1)';
    Zi = zeros(1,Np);

    Xgrid = [Xi;Yi;Zi];


% All the point coordinates (on the image, and in 3D) - for global optimization:

    x = grid_pts;
    X = Xgrid;


    if (b_display)
            ind_corners = [1 n_sq_x+1 (n_sq_x+1)*n_sq_y+1 (n_sq_x+1)*(n_sq_y+1)]; % index of the 4 corners
        ind_orig = (n_sq_x+1)*n_sq_y + 1;
        xorig = grid_pts(1,ind_orig);
        yorig = grid_pts(2,ind_orig);
        dxpos = mean([grid_pts(:,ind_orig) grid_pts(:,ind_orig+1)]');
        dypos = mean([grid_pts(:,ind_orig) grid_pts(:,ind_orig-n_sq_x-1)]');

        x_box_kk = [grid_pts(1,:)-(wintx+.5);grid_pts(1,:)+(wintx+.5);grid_pts(1,:)+(wintx+.5);grid_pts(1,:)-(wintx+.5);grid_pts(1,:)-(wintx+.5)];
        y_box_kk = [grid_pts(2,:)-(winty+.5);grid_pts(2,:)-(winty+.5);grid_pts(2,:)+(winty+.5);grid_pts(2,:)+(winty+.5);grid_pts(2,:)-(winty+.5)];

        %%%%%%%%%%%%%% display the corners
        x1= x(1); x2 = x(2); x3 = x(3); x4 = x(4);
        y1= y(1); y2 = y(2); y3 = y(3); y4 = y(4);

        % Find center:
        p_center = cross(cross([x1;y1;1],[x3;y3;1]),cross([x2;y2;1],[x4;y4;1]));
        x5 = p_center(1)/p_center(3);
        y5 = p_center(2)/p_center(3); 
        % center on the X axis:
        x6 = (x3 + x4)/2;
        y6 = (y3 + y4)/2;
        % center on the Y axis:
        x7 = (x1 + x4)/2;
        y7 = (y1 + y4)/2;
        % Direction of displacement for the X axis:
        vX = [x6-x5;y6-y5];
        vX = vX / norm(vX);
        % Direction of displacement for the X axis:
        vY = [x7-x5;y7-y5];
        vY = vY / norm(vY);
        % Direction of diagonal:
        vO = [x4 - x5; y4 - y5];
        vO = vO / norm(vO); 
        delta = 30;

        figure(3);
        imshow(I); colormap(map); hold on;
        plot(grid_pts(1,:)+1,grid_pts(2,:)+1,'r+');
        plot(x_box_kk+1,y_box_kk+1,'-b');
        plot(grid_pts(1,ind_corners)+1,grid_pts(2,ind_corners)+1,'mo');
        plot(xorig+1,yorig+1,'*m');
        h = text(xorig+3*delta*vO(1),yorig+3*delta*vO(2),num2str(board.BoardIdx));
        set(h,'Color','y','FontSize',20);
        h2 = text(dxpos(1)+delta*vX(1),dxpos(2)+delta*vX(2),'dX');
        set(h2,'Color','g','FontSize',14);
        h3 = text(dypos(1)+delta*vY(1),dypos(2)+delta*vY(2),'dY');
        set(h3,'Color','g','FontSize',14);
        xlabel('Xc (in camera frame)');
        ylabel('Yc (in camera frame)');
        title('Extracted corners');
        zoom on;
        drawnow;
        hold off;
    end
end


