function board_set = getBoardFromIm(im_path,...
    checkerboard_size, checkerboard_grid_size,...
    sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
    b_first_grid_color, possible_checkerboard_num, over_s, b_boardidx,...
    b_display, figure_index)

    N_figure = 1;
    if (~exist('b_display','var')),  b_display = 0;  end
    if (~exist('figure_index','var')),  figure_index = 1:N_figure;    end
    if ( (numel(figure_index)<N_figure) && (figure_index~=0) ),  
        disp(['N_figure = ', num2str(N_figure)]);
        figure_index = 1:N_figure;    
    end

    disp('working im = ');
    disp(im_path);

    if (b_boardidx)
        fprintf('the task is to find board for rotation calib.\n');
    else
        fprintf('the task is to find main checkerboard.\n');
    end


    fprintf(['  checkerboard_size = ', num2str(checkerboard_size),'.\n',...
        '  checkerboard_grid_size = ', num2str(checkerboard_grid_size),'.\n',...
        '  sigma = ', num2str(sigma),'.\n',...
        '  min_l_th = ', num2str(min_l_th),'.\n',...
        '  max_k_th = ', num2str(max_k_th),'.\n',...
        '  var_k_th = ', num2str(var_k_th),'.\n',...
        '  var_l_th = ', num2str(var_l_th),'.\n',...
        '  over_s = ', num2str(over_s),'.\n']);

    if b_first_grid_color
        fprintf('  first grid is black.\n');
    else
        fprintf('  first grid is white.\n');
    end
    
    %%
    dX = checkerboard_grid_size;
    dY = checkerboard_grid_size;
    checkerboard.dX = dX;
    checkerboard.dY = dY;
    nx = checkerboard_size(1);
    ny = checkerboard_size(2);
    checkerboard.nx = nx;
    checkerboard.ny = ny;
    temp_x = 0:(ny-1-1);
    temp_y = 0:(nx-1-1);
    [X_x,X_y] = meshgrid(temp_x,temp_y);
    X(:,:,1) = X_x;
    X(:,:,2) = X_y;
    X(:,:,3) = 0;
    b_X = reshape(X,(nx-1)*(ny-1),3)' * checkerboard_grid_size;
    checkerboard.b_X = b_X;

    %%
    board_set = [];

    im = getImFromPath(im_path);
    boards = checkerboard_detector(im, checkerboard_size, b_display,...
        sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
        b_first_grid_color, possible_checkerboard_num, over_s, b_boardidx);

    fprintf([num2str(numel(boards)), ' board(s) are detected.\n']);
    for board_i = 1:numel(boards)
        board = boards(board_i);
        board.filepath = im_path;
        board.checkerboard = checkerboard;
        board_set = [board_set,board];
    end    
end