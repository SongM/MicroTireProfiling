function board_set = getBoardFromFolder(folder_name,...
    checkerboard_size, checkerboard_grid_size,...
    sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
    b_first_grid_color, possible_checkerboard_num, over_s, b_boardidx,...
    b_display, figure_index, im_ind)

    tic;
    N_figure = 1;
    if (~exist('b_display','var')),  b_display = 0;  end
    if (~exist('figure_index','var')),  figure_index = 1:N_figure;    end
    if ( (numel(figure_index)<N_figure) && (figure_index~=0) ),  
        disp(['N_figure = ', num2str(N_figure)]);
        figure_index = 1:N_figure;    
    end
    
    if (~exist('im_ind','var')),  im_ind = 0;  end


    disp('working folder = ');
    disp(folder_name);
    
    

    if (b_boardidx)
        task = 'laser';
    else
        task = 'cam_laser_calib';
%         fprintf([num2str(N_image),' images are used for camera calibration.\n']);
    end
    
%     im_file_list = dir([folder_name,'\',task,'*mat']);
    im_file_list = dir([folder_name,'\',task,'\*JPG']);
%     im_file_list = dir([folder_name,'\',task,'\*DNG']);
    N_image = numel(im_file_list);
    im_proc = 1:N_image;
    if (sum(im_ind)~=0),    N_image = numel(im_ind);   im_proc = im_ind;    b_display = 1;    end

    fprintf( [num2str(N_image), ' images are used for ', task, ' calibration.\n'] );
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
    
    board_set = [];
    for im_i=im_proc
        fprintf(1, [task,', working on image_', num2str(im_i), '/', num2str(N_image),'. ']);
        im_name = im_file_list(im_i).name;
        im_path = [folder_name,'\',task,'\',im_name];

        im = getImFromPath(im_path);

        boards = checkerboard_detector(im, checkerboard_size, b_display,...
            sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
            b_first_grid_color, possible_checkerboard_num, over_s, b_boardidx);
        
        
        time = toc;
        fprintf([num2str(numel(boards)), ' board(s) are detected. ',...
            'time remain = ',num2str((N_image-im_i)*time/im_i),'s.\n']);
        for board_i = 1:numel(boards)
            board = boards(board_i);
            board.filepath = im_path;
            board.checkerboard = checkerboard;
            board_set = [board_set,board];
        end
    end
    
    fprintf([num2str(numel(board_set)), 'boards are detected in ', num2str(N_image), ' images.' ]);
    toc;
end