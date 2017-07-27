function boards = checkerboard_detector(im, checkerboard_size, ...
    b_display, sigma, min_l_th, max_k_th, var_k_th, var_l_th,...
    b_first_grid_color, possible_checkerboard_num, over_s, b_boardidx)
    
    B_DISPLAY = 0;
    B_DEBUG = 0;
% 
%     if (nargin < 7),    var_l_th = 100; end,
%     if (nargin < 6),    var_k_th = 0.1; end,
%     if (nargin < 5),    min_l_th = 2;   end,
%     if (nargin < 4),    sigma = 2;  end,
%     if (nargin < 3),    b_display = 0;  end,

    if (b_boardidx)
        up_percent = 0.4;
        I = im(1:(size(im,1)*up_percent),:,2);
    else
        I = im(:,:,2);
    end
    
    % get possible corner candidates
    [points, scores] = getPointsFromImage(I, sigma, ...
        checkerboard_size, possible_checkerboard_num, over_s, B_DISPLAY);
    [lines_idx, lines_energy] = getLinesFromPoints(points, scores,...
        checkerboard_size, min_l_th, max_k_th, var_k_th, var_l_th, B_DISPLAY, I);
    boards = getBoardsFromLines(points, scores, lines_idx, lines_energy,...
        checkerboard_size, var_k_th, var_l_th, B_DISPLAY, I, B_DEBUG); 
    
    if (numel(boards) > 0)
        for board_i = 1:numel(boards)
            board = boards(board_i);
            flag_board_isvalid(board_i) = checkBoardsFirstGrid(board, I , b_first_grid_color, B_DEBUG);
        end
        
        boards = boards(flag_board_isvalid==1);
        for board_i = 1:numel(boards)
            boards(board_i).BoardIdx = getBoardIdx(boards(board_i),I,b_boardidx);
            if (b_boardidx)
                boards(board_i).cornerx([1,2]) = boards(board_i).cornerx([1,2]) - boards(board_i).v2(1)*2;
                boards(board_i).cornery([1,2]) = boards(board_i).cornery([1,2]) - boards(board_i).v2(2)*2;
            end
        end
    end
    
    if (b_display)

        figure();imshow(I);hold on;
        scatter(points(:,1),points(:,2),scores/max(scores(:))*20);
        for i = 1:size(lines_idx,1)
            text(points(lines_idx(i,1),1), points(lines_idx(i,1),2),...
                num2str(i), 'Color', 'y', 'FontSize', 12);
%             scatter(points(lines_idx(i,:),1),points(lines_idx(i,:),2),10,'y','filled');
            plot(points(lines_idx(i,:),1),points(lines_idx(i,:),2),'y');
        end


        for board_i = 1:numel(boards)
            board = boards(board_i);
            cornerx = [board.cornerx,board.cornerx(1)];
            cornery = [board.cornery,board.cornery(1)];

            plot(cornerx,cornery,'*-','Color',[1,0.5,0],'LineWidth',1);
            text(board.cornerx(1),board.cornery(1),...
                num2str(board.BoardIdx),'Color',[1,0.5,0],'FontSize',12);
        end
        drawnow;
    end


end


function  [points, scores] = ...
    getPointsFromImage(I, sigma, ...
    checkerboard_size, possible_checkerboard_num, over_s, b_display)
    
    [cxy, ~, ~, ~] = ...
    vision.internal.calibration.checkerboard.secondDerivCornerMetric(I, sigma);

    peakThreshold = 0.2;
    points = [];
    while (size(points,1)<checkerboard_size(1)*checkerboard_size(2)*over_s*possible_checkerboard_num)
        points = vision.internal.calibration.checkerboard.find_peaks(cxy, peakThreshold);
        peakThreshold = peakThreshold/1.2;
    end
%     
%     while (size(points,1)>checkerboard_size(1)*checkerboard_size(2)*(possible_checkerboard_num+1)*10)
%         points = vision.internal.calibration.checkerboard.find_peaks(cxy, peakThreshold);
%         peakThreshold = peakThreshold*1.2;
%     end
    
    scores = cxy(sub2ind(size(cxy), points(:, 2), points(:, 1)));
    
    if (b_display)
         figure();imshow(I);hold on;
         scatter(points(:,1),points(:,2),scores/max(scores(:))*20);
    end
   
    
end

function [lines_idx, lines_energy] =...
    getLinesFromPoints(points, scores,...
    checkerboard_size, min_l_th, max_k_th, var_k_th, var_l_th, b_display, I)

    if (b_display)
        figure();imshow(I);hold on;
        scatter(points(:,1),points(:,2),scores/max(scores(:))*20);        
        for i = 1:size(points,1)
            text(points(i,1),points(i,2),num2str(i),'Color','red','FontSize',10);
        end
    end
    
    % rearrage the points wrt to their scores, so that we can start
    % searching the line from the most possible corner.
    point_seedsIdx = 1:size(points, 1);
    [~, sortedIdx] = sort(scores(point_seedsIdx), 'descend');
    point_seedsIdx = point_seedsIdx(sortedIdx);
    
    lines_energy = [];
    lines_idx = [];
    for point_seed_i = 1:numel(point_seedsIdx)
        point_seed_idx = point_seedsIdx(point_seed_i);
        if (point_seed_idx == 0)
            continue;
        end
        line = getLineFromPointsWRTSeed(point_seed_idx,...
            points, scores, checkerboard_size(2)-1,...
            min_l_th, max_k_th, var_k_th, var_l_th);
        if (line.energy > 0)
            lines_energy = [lines_energy; line.energy];
            lines_idx = [lines_idx; line.idx];
            % remove the line points from seed candidates.
            for line_idx_i = line.idx
                point_seedsIdx( point_seedsIdx == line_idx_i ) = 0;
            end
            if (b_display)
                text(points(line.idx(1),1), points(line.idx(1),2),...
                    num2str(size(lines_idx,1)), 'Color', 'y', 'FontSize', 12);
%                 scatter(points(line.idx,1),points(line.idx,2),10,'y','filled');
                plot(points(line.idx,1),points(line.idx,2),'y');
            end
        end
    end
end

function line =...
    getLineFromPointsWRTSeed(point_seed_idx,...
    points, scores, checkerboard_size_x,...
    min_l_th, max_k_th, var_k_th, var_l_th)

    N_grid_x = checkerboard_size_x-1;
    factor_l_th = 1.1;
    % calculate k and l of all the points to the seed point
    seed_coor = points(point_seed_idx,:);
    seed_x = seed_coor(1);
    seed_y = seed_coor(2);

    points_x = points(:,1);
    points_y = points(:,2);
    
    k = (points_y - seed_y) ./ (points_x - seed_x + realmin);
    l = sqrt( (points_x - seed_x).^2 + (points_y - seed_y).^2 );
    
    % Get rid of outliers wrt the seed
    cand_idx = ones(1,numel(k));
        % distance far away from the seed
    temp_l = sort(l);
    dist_base = mean(temp_l(2:4));
    temp_l_th = factor_l_th * N_grid_x * dist_base;
    
    
    cand_idx(l<min_l_th) = 0;
    cand_idx(l>temp_l_th) = 0;
        % slope larger than max_k_th
    cand_idx(abs(k)>max_k_th) = 0;
%     cand_idx(k<-1) = 0;
        % the seed itself is not counted
    cand_idx(point_seed_idx) = 0;
    
    l_cand = l(cand_idx==1);
    
    % if no enough candidate, then return.
    if (numel(l_cand)<N_grid_x)
        line_energy = -1;
        line_idx = [];
    else
        if (numel(l_cand)==N_grid_x)
            s_cand = scores(cand_idx==1);
            line_energy = min(s_cand);
            unsorted_line_idx = [point_seed_idx , find(cand_idx==1)];
        else
            % find the most possible points for the line which contains the
            % seed, basically use the information from the similarity of k.
            k_cand = k(cand_idx==1);
            s_cand = scores(cand_idx==1);
            dist_k = dist(k_cand');
            temp_sorted_dist_k = sort(dist_k);
                % It is assumed that the ks should be similar to each
                % other, so that larger dist_k can filtered.
            dist_k_th = mean(mean( temp_sorted_dist_k(1:N_grid_x,:) ));
            dist_k_filtered_result = ones(size(temp_sorted_dist_k));
            dist_k_filtered_result(dist_k>dist_k_th) = 0;
                % the scores are taken into consideration.
            scores_of_dist_k_filtered_result =...
                dist_k_filtered_result .* repmat(s_cand,1,numel(s_cand));
                % find the most possible line_mate of each k.
            [scores_result_val,scores_result_idx] =...
                sort(scores_of_dist_k_filtered_result,'descend');
                % the scores of the seed itself is taken into
                % consideration. and only the first N_grid_x scores are
                % taken into consideration.
            scores_for_each_k =...
                min(scores_result_val(1:N_grid_x,:)).*s_cand';
                % choose the max scores as the vector.
            [line_energy,winner_of_k] = max(scores_for_each_k);
            idx = find(cand_idx==1);
            unsorted_line_idx =...
                [point_seed_idx , idx(scores_result_idx(1:N_grid_x,winner_of_k))];
        end
        % sort the line points based on their x-coordinate.
        [~,sorted_line_idx_order] = sort(points(unsorted_line_idx,1));
        line_idx = unsorted_line_idx(sorted_line_idx_order);
        % verify the line
        points_line = points(line_idx,:);
        vector = points_line(2:(N_grid_x+1),:)-points_line(1:N_grid_x,:);
        vector_k = atan( vector(:,2) ./ (vector(:,1)+realmin) );
        vector_l = sqrt( vector(:,1).^2 + vector(:,2).^2 );
        var_vector_k = var(vector_k);
        var_vector_l = var(vector_l);
        if (mean(vector_l)<min_l_th)
            line_energy = -1;
        end
        if (var_vector_k>var_k_th)
            line_energy = -1;
        end
        if (var_vector_l>var_l_th)
            line_energy = -1;
        end
    end
    line.energy=line_energy;
    line.idx=line_idx;
end

function [boards, boards_idx_of_line] =...
    getBoardsFromLines(points, scores, lines_idx, lines_energy, ...
    checkerboard_size, var_k_th, var_l_th, b_display, I, B_DEBUG) 

    if (b_display)
        figure();imshow(I);hold on;
        scatter(points(:,1),points(:,2),scores/max(scores(:))*20);        
        for i = 1:size(points,1)
            text(points(i,1),points(i,2),num2str(i),'Color','red','FontSize',10);
        end
        
        for i = 1:size(lines_idx,1)
            text(points(lines_idx(i,1),1), points(lines_idx(i,1),2),...
                num2str(i), 'Color', 'y', 'FontSize', 12);
            scatter(points(lines_idx(i,:),1),points(lines_idx(i,:),2),10,'y','filled');
        end
        
    end

    line_seedsIdx = 1:size(lines_idx,1);
    
    boards = [];
    boards_energy = [];
    boards_idx_of_line = [];
    for line_seed_i = 1:numel(line_seedsIdx)
        line_seed_idx = line_seedsIdx(line_seed_i);
        if (line_seed_idx==0)
            continue;
        end
        
        [board , board_idx_of_line] = getBoardFromLinesWRTSeed(...
        line_seed_idx, lines_idx, lines_energy, points,...
        checkerboard_size-1, var_k_th, var_l_th, B_DEBUG);

        if (board.energy > 0)
            boards_energy = [boards_energy; board.energy];
            boards_idx_of_line = [boards_idx_of_line, board_idx_of_line];
            boards = [boards,board];
            % remove the line points from seed candidates.
            for board_idx_i = board_idx_of_line'
                line_seedsIdx( line_seedsIdx == board_idx_i ) = 0;
            end
            if (b_display)
                cornerx = [board.cornerx,board.cornerx(1)];
                cornery = [board.cornery,board.cornery(1)];
                plot(cornerx,cornery,'*-','Color',[1,0.5,0],'LineWidth',1);
            end
        end   
    end
end
        
function [board , board_idx_of_line] =...
    getBoardFromLinesWRTSeed(line_seed_idx, lines_idx, lines_energy,...
    points,  checkerboard_size, var_k_th, var_l_th, B_DEBUG)

if B_DEBUG
    fprintf(['\n line_seed_idx=',num2str(line_seed_idx),', ']);
end
    %%%
    N_grid_x = checkerboard_size(2)-1;
    N_grid_y = checkerboard_size(1)-1;

    % only the first point of each line will be calculated for the vector
    % of the other direction.
    first_point_array = points(lines_idx(:,1),:);

    % calculate k and l of all the first-points to the seed first-point
    seed_coor = first_point_array(line_seed_idx,:);
    seed_x = seed_coor(1);
    seed_y = seed_coor(2);
    points_x = first_point_array(:,1);
    points_y = first_point_array(:,2);
    
    k = (points_y - seed_y) ./ (points_x - seed_x + realmin);
    l = sqrt( (points_x - seed_x).^2 + (points_y - seed_y).^2);
    % Get rid of outliers wrt the seed
    cand_idx = ones(1,numel(k));
        % slope less than 1 (45 degree)
    cand_idx( (k<1) & (k>-1) ) = 0;
    l_cand = l(cand_idx==1);

    % if no enough candidates, return -1    
    if (numel(l_cand)<N_grid_y)
        board.energy = -1;
        board_idx_of_line = [];
    else
        if ( numel(l_cand) == N_grid_y )
            s_cand = lines_energy(cand_idx==1);
            board.energy = sum(s_cand);
            unsorted_board_idx_of_line = [line_seed_idx ; find(cand_idx==1)'];
        else    
            % find the most possible first_points to form the line which
            % contains the seed. Use theta instead of k because the line is
            % vertiacl.
            theta_cand = atan(k(cand_idx==1));
            theta_cand(theta_cand<0) = theta_cand(theta_cand<0) + pi;
            s_cand = lines_energy(cand_idx==1);            
            dist_theta = dist(theta_cand');
            temp_sorted_dist_theta = sort(dist_theta);
                % It is assumed that the thetas should be similar to each
                % other, so that larger dist_theta can filtered.
            dist_theta_th = mean(mean(...
                temp_sorted_dist_theta(1:N_grid_y,:) ));
            dist_theta_filtered_result = ones(size(temp_sorted_dist_theta));
            dist_theta_filtered_result(dist_theta>dist_theta_th) = 0;
                % the scores are taken into consideration.
            scores_of_dist_theta_filtered_result =...
                dist_theta_filtered_result .* repmat(s_cand,1,numel(s_cand));
                % find the most possible line_mate of each k.
            [scores_result_val,scores_result_idx] =...
                sort(scores_of_dist_theta_filtered_result,'descend');
                % the scores of the seed itself is taken into
                % consideration. and only the first N_grid_y scores are
                % taken into consideration.
            scores_for_each_theta =...
                sum(scores_result_val(1:N_grid_y,:)).*s_cand';
                % choose the max scores as the vector.
            [board.energy, winner_of_theta] = max(scores_for_each_theta);
            idx = find(cand_idx==1);
            unsorted_board_idx_of_line =...
                [line_seed_idx ; idx(scores_result_idx(1:N_grid_y,winner_of_theta))'];
        end
                
        % sort the lines based on their y-coordinate.
        [~,sorted_board_idx_order] =...
            sort(first_point_array(unsorted_board_idx_of_line,2));
        % get the board idx
        board_idx_of_line = unsorted_board_idx_of_line(sorted_board_idx_order);
        board_idx_of_point = lines_idx(board_idx_of_line,:);

        % verify the board
        points_x = points(:,1);
        points_y = points(:,2);


        points_board_x = points_x(board_idx_of_point);
        points_board_y = points_y(board_idx_of_point);
        v1_x = points_board_x(:, 2:(N_grid_x+1))...
            - points_board_x(:, 1:N_grid_x);
        v1_y = points_board_y(:, 2:(N_grid_x+1))...
            - points_board_y(:, 1:N_grid_x);
        v2_x = points_board_x(2:(N_grid_y+1), :)...
            - points_board_x(1:N_grid_y, :);
        v2_y = points_board_y(2:(N_grid_y+1), :)...
            - points_board_y(1:N_grid_y, :);

        v1_theta = atan(v1_y ./ (v1_x+realmin) );
        v2_theta = atan(v2_y ./ (v2_x+realmin) );
        if (mean(abs(v2_theta(:))) > pi/2*0.9)
            v2_theta(v2_theta<0) = v2_theta(v2_theta<0)+pi;
        end
        v1_l = sqrt(v1_x.^2 + v1_y.^2);
        v2_l = sqrt(v2_x.^2 + v2_y.^2);

        var_v1_theta = var( v1_theta(:) );
        var_v2_theta = var( v2_theta(:) );
        var_v1_l = var( v1_l(:) );
        var_v2_l = var( v2_l(:) );
        if ( var_v1_theta > atan(var_k_th) )
            if B_DEBUG
                fprintf(['\n line_idx', num2str(board_idx_of_line'),...
                    ', var_v1_theta=',num2str(var_v1_theta),'\n']);
            end
            board.energy = -1;return;
        end
        if ( var_v1_l > var_l_th )
            if B_DEBUG
                fprintf(['\n line_idx', num2str(board_idx_of_line'),...
                    ', var_v1_l=',num2str(var_v1_l),'\n']);
            end
            board.energy = -1;return;
        end
        if ( var_v2_theta > atan(var_k_th) )
            if B_DEBUG
                fprintf(['\n line_idx', num2str(board_idx_of_line'),...
                    ', var_v2_theta=',num2str(var_v2_theta),'\n']);
            end
            board.energy = -1;return;
        end
        if ( var_v2_l > var_l_th)
            if B_DEBUG           
                fprintf(['\n line_idx', num2str(board_idx_of_line'),...
                    ', var_v2_l=',num2str(var_v2_l),'\n']);
            end
            board.energy = -1;return;
        end
        board.idx = board_idx_of_point;
        board.v1(1) = mean( v1_x(:) );
        board.v1(2) = mean( v1_y(:) );

        board.v2(1) = mean( v2_x(:) );
        board.v2(2) = mean( v2_y(:) );

        board.coor_x = points_board_x;
        board.coor_y = points_board_y;

        % get board corners
        cornerx(1) = points_board_x(1,1) - board.v1(1) - board.v2(1);
        cornery(1) = points_board_y(1,1) - board.v1(2) - board.v2(2);
        cornerx(2) = points_board_x(1,end) + board.v1(1) - board.v2(1);
        cornery(2) = points_board_y(1,end) + board.v1(2) - board.v2(2);
        cornerx(3) = points_board_x(end,end) + board.v1(1) + board.v2(1);
        cornery(3) = points_board_y(end,end) + board.v1(2) + board.v2(2);
        cornerx(4) = points_board_x(end,1) - board.v1(1) + board.v2(1);
        cornery(4) = points_board_y(end,1) - board.v1(2) + board.v2(2);
        board.cornerx = cornerx;
        board.cornery = cornery;
    end
end


function flag_board_isvalid = checkBoardsFirstGrid(board, I , b_first_grid_color, B_DEBUG)
    % check if the first grid is black to determin if the board is valid.
    % p1 p2 p3 p4
    % p5 p6 p7 p8
    % p9 p10
    p1_x = board.coor_x(1,1);
    p2_x = board.coor_x(1,2); 
    p3_x = board.coor_x(1,3); 
    p4_x = board.coor_x(1,4); 
    p5_x = board.coor_x(2,1); 
    p6_x = board.coor_x(2,2); 
    p7_x = board.coor_x(2,3); 
    p8_x = board.coor_x(2,4); 
    p9_x = board.coor_x(3,1); 
    p10_x = board.coor_x(3,2); 
    
    p1_y = board.coor_y(1,1);
    p2_y = board.coor_y(1,2); 
    p3_y = board.coor_y(1,3); 
    p4_y = board.coor_y(1,4); 
    p5_y = board.coor_y(2,1); 
    p6_y = board.coor_y(2,2); 
    p7_y = board.coor_y(2,3); 
    p8_y = board.coor_y(2,4); 
    p9_y = board.coor_y(3,1); 
    p10_y = board.coor_y(3,2); 
    % grid11 grid12
    % grid21
    grid11_x = [p1_x, p2_x, p6_x, p5_x];
    grid11_y = [p1_y, p2_y, p6_y, p5_y];
    grid11_roi = roipoly(I, grid11_x, grid11_y);
    grid11_color = mean(I(grid11_roi));

    grid12_x = [p2_x, p3_x, p7_x, p6_x];
    grid12_y = [p2_y, p3_y, p7_y, p6_y];
    grid12_roi = roipoly(I, grid12_x, grid12_y);
    grid12_color = mean(I(grid12_roi));

    grid21_x = [p5_x, p6_x, p10_x, p9_x];
    grid21_y = [p5_y, p6_y, p10_y, p9_y];
    grid21_roi = roipoly(I, grid21_x, grid21_y);
    grid21_color = mean(I(grid21_roi));
    if b_first_grid_color
        flag_board_isvalid = (grid11_color<grid12_color)&(grid11_color<grid21_color);
    else
        flag_board_isvalid = (grid11_color>grid12_color)&(grid11_color>grid21_color);
    end
    if (~flag_board_isvalid)
        if (B_DEBUG)
            fprintf('did not pass the first grid check.\n');
        end
    end
end
% 
function BoardIdx = getBoardIdx(board,I,b_boardidx)
    if(b_boardidx)
        boardIdx_grid_N = 3;
        BoardIdx = 0;
        for boardIdx_grid_i = 1:boardIdx_grid_N
            %   1  2 
            %    56       
            %    78xXxX
            %   3  4
            %
            %xXxX0XxXx
            
            p0(1) = board.coor_x(1,boardIdx_grid_i*2+3);
            p0(2) = board.coor_y(1,boardIdx_grid_i*2+3);
            t = 0.3;
            v1 = board.v1;
            v2 = board.v2;
            p5 = p0 + 0*v1 - 3*v2;
            p6 = p0 + 1*v1 - 3*v2;
            p7 = p0 + 0*v1 - 2*v2;
            p8 = p0 + 1*v1 - 2*v2;
            p1 = p5 - t*v1 - t*v2;
            p2 = p6 + t*v1 - t*v2;
            p3 = p7 - t*v1 + t*v2;
            p4 = p8 + t*v1 + t*v2;
            in = [p5;p6;p8;p7];
            in_roi = roipoly(I, in(:,1), in(:,2));
            in_color = mean(I(in_roi));

            out = [p1;p2;p4;p3];
            out_roi = roipoly(I, out(:,1), out(:,2));
            out_color = mean(I(out_roi));
            BoardIdx = BoardIdx*2;
            if(in_color< out_color)
                BoardIdx = BoardIdx+1;
            end
        end
    else
        BoardIdx = 0;
    end
end
