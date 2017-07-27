function displayBoard3D(board,color)
    if (nargin<2)
        color = [1,0,0];
    end
    corner3 = board.camera_plane.Corner;
    corner3(:,5) = corner3(:,1);

        c_grid(:,:,1) = board.camera_plane.coor_X;
        c_grid(:,:,2) = board.camera_plane.coor_Y;
        c_grid(:,:,3) = board.camera_plane.coor_Z;
    
    plot3(corner3(1,:),corner3(2,:),corner3(3,:), 'Color', color);hold on;
    plot3(board.c_X(1,:),board.c_X(2,:),board.c_X(3,:),'o', 'Color', color);
    plot3(c_grid(1,1,1),c_grid(1,1,2),c_grid(1,1,3),'+','Color', color,'LineWidth',2);
    plot3(c_grid(:,:,1),c_grid(:,:,2),c_grid(:,:,3), 'Color', color);
    plot3(c_grid(:,:,1)',c_grid(:,:,2)',c_grid(:,:,3)', 'Color', color);

    if (board.BoardIdx ~= 0)
        c_V1 = board.camera_plane.V1;
        c_V2 = board.camera_plane.V2;
%         nx = board.checkerboard.nx;
        for square_i = 1:3
            square_corner = [];
            square_corner(1,:) = reshape(c_grid(1,square_i*2+3,:),1,3) + 1*c_V1 - 2*c_V2;
            square_corner(2,:) = reshape(c_grid(1,square_i*2+3,:),1,3) + 1*c_V1 - 3*c_V2;
            square_corner(3,:) = reshape(c_grid(1,square_i*2+3,:),1,3) + 0*c_V1 - 3*c_V2;
            square_corner(4,:) = reshape(c_grid(1,square_i*2+3,:),1,3) + 0*c_V1 - 2*c_V2;
            square_corner(5,:) = square_corner(1,:);
            square_corner = square_corner';
            B = board.BoardIdx;
            if     (floor(mod(B,2^(4-square_i))/(2^(3-square_i))))
                fill3(square_corner(1,:),square_corner(2,:),square_corner(3,:), color);
            else
                plot3(square_corner(1,:),square_corner(2,:),square_corner(3,:), 'Color', color);
            end
        end
    end
end