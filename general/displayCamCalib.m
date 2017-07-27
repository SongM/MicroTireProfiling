function displayCamCalib(CamCalib,board_calib,show_camera,figure_index)

    if (~exist('board_calib','var')),  
        board_calib = [];  
    end
    
    if (~exist('show_camera','var')),  show_camera = 1;  end
    if (~exist('figure_index','var')),  figure_index = 100;    end

    %%%%%%%%%%%%%%%%%%%% SHOW EXTRINSIC RESULTS %%%%%%%%%%%%%%%%%%%%%%%%
    % Color code for each image:
    colors = 'brgkcm';
    
    fc = CamCalib.fc;
    cc = CamCalib.cc;
    % kc = CamCalib.fc;
    alpha_c = CamCalib.alpha_c;
    % KK = CamCalib.KK;

    
    if (isempty(board_calib))
        dX = CamCalib.display.dX;
        dY = dX;
        IP = CamCalib.display.IP;
        BASE = CamCalib.display.BASE;
        
    else
        dX = board_calib(1).checkerboard.dX;
        dY = board_calib(1).checkerboard.dY;
        n_sq_x = board_calib(1).checkerboard.nx-2;
        n_sq_y = board_calib(1).checkerboard.ny-2;
        [ny,nx,~]=size(getImFromPath(board_calib(1).filepath));

        IP = 5*dX*[1 -alpha_c 0;0 1 0;0 0 1]*[1/fc(1) 0 0;0 1/fc(2) 0;0 0 1]*[1 0 -cc(1);0 1 -cc(2);0 0 1]*[0 nx-1 nx-1 0 0 ; 0 0 ny-1 ny-1 0;1 1 1 1 1];
        BASE = 5*dX*([0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]);
        IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);
    end


    figure(figure_index);
    a = 50;
    b = 20;
    if show_camera,
        figure(figure_index);
        plot3(BASE(1,:),BASE(3,:),-BASE(2,:),'b-','linewidth',2);
        hold on;
        plot3(IP(1,:),IP(3,:),-IP(2,:),'r-','linewidth',2);
        text(6*dX,0,0,'X_c');
        text(-dX,5*dX,0,'Z_c');
        text(0,0,-6*dX,'Y_c');
        text(-dX,-dX,dX,'O_c');
    else
        figure(figure_index);
        clf;
        hold on;
    end;

    for kk = 1:numel(board_calib)
        board = board_calib(kk);
        XX_kk = board.b_X;
        N_kk = size(XX_kk,2);

        omc_kk = board.b2c_omc;
        Tc_kk = board.b2c_Tc;
        R_kk = rodrigues(omc_kk);

        YY_kk = R_kk * XX_kk + Tc_kk * ones(1,length(XX_kk));

        uu = [-dX;-dY;0]/2;
        uu = R_kk * uu + Tc_kk; 

        YYx = zeros(n_sq_x+1,n_sq_y+1);
        YYy = zeros(n_sq_x+1,n_sq_y+1);
        YYz = zeros(n_sq_x+1,n_sq_y+1);

        YYx(:) = YY_kk(1,:);
        YYy(:) = YY_kk(2,:);
        YYz(:) = YY_kk(3,:);
    
        figure(figure_index);
        hhh= mesh(YYx,YYz,-YYy);
        set(hhh,'edgecolor',colors(rem(kk-1,6)+1),'linewidth',1); %,'facecolor','none');
        %plot3(YY_kk(1,:),YY_kk(3,:),-YY_kk(2,:),['o' colors(rem(kk-1,6)+1)]);
        text(uu(1),uu(3),-uu(2),num2str(kk),'fontsize',14,'color',colors(rem(kk-1,6)+1));
    end

    figure(figure_index);rotate3d on;
    axis('equal');
    title('Extrinsic parameters (camera-centered)');
    
    view(a,b);
    grid on;
    hold off;
    axis vis3d;
    axis tight;
    set(figure_index,'color',[1 1 1]);
    if ~show_camera,
        xlabel('X_c');
        ylabel('Z_c');
        zlabel('<-- Y_c');
    end;
    set(figure_index,'Name','3D');
    if show_camera,
        h_switch2 = uicontrol('Parent',figure_index,'Units','normalized', 'Callback','displayCamCalib(CamCalib,board_calib,0);', 'Position',[1-.30 0  .30  .04],'String','Remove camera reference frame','fontsize',8,'fontname','clean','Tag','Pushbutton1');
    else
        h_switch2 = uicontrol('Parent',figure_index,'Units','normalized', 'Callback','displayCamCalib(CamCalib,board_calib,1);', 'Position',[1-.30 0  .30  .04],'String','Add camera reference frame','fontsize',8,'fontname','clean','Tag','Pushbutton1');
    end;

    figure(figure_index);
    rotate3d on;
end