function displayBoard2D(board)

    n_sq_x = board.checkerboard.nx - 2;
    
    grid_pts = board.p_x;
    
    x_orig = grid_pts(1,end-n_sq_x)+1;
    y_orig = grid_pts(2,end-n_sq_x)+1;
    
    cornerx = [board.pixel_plane.corner(1,:),board.pixel_plane.corner(1,1)];
    cornery = [board.pixel_plane.corner(2,:),board.pixel_plane.corner(2,1)];

        
    
    delta = 50;
    im = getImFromPath(board.filepath);
    imshow(im);
    hold on;
   
    plot(grid_pts(1,:)+1,grid_pts(2,:)+1,'r+');
    plot(board.p_x_proj(1,:)+1,board.p_x_proj(2,:)+1,'yo');
    s = 30;
    quiver(grid_pts(1,:)+1,grid_pts(2,:)+1,...
        s*board.error.p_x_proj_error(1,:),s*board.error.p_x_proj_error(2,:),0,'Color','g');

    plot(x_orig,y_orig,'*m');
    h = text(x_orig-delta,y_orig-delta,num2str(board.BoardIdx));
    set(h,'Color','y','FontSize',20);
    h2 = text(x_orig-delta,y_orig+delta,'dX');
    set(h2,'Color','g','FontSize',14);
    h3 = text(x_orig+delta,y_orig-delta,'dY');
    set(h3,'Color','g','FontSize',14);

    plot(cornerx,cornery,'*-','Color',[1,0.5,0],'LineWidth',1);

    title(['std error: ', num2str(std(board.error.p_x_proj_error'))]);
    zoom on;
    drawnow;
    hold off;
end

