function displayCamera(CamCalib)
    dX = CamCalib.display.dX;
    IP = CamCalib.display.IP;
    BASE = CamCalib.display.BASE;

    hold on;
    plot3(BASE(1,:),BASE(2,:),BASE(3,:),'b-','linewidth',2);
    plot3(IP(1,:),IP(2,:),IP(3,:),'r-','linewidth',2);
    text(6*dX,0,0,'X_c');
    text(-dX,0,5*dX,'Z_c');
    text(0,6*dX,0,'Y_c');
    text(-dX,dX,-dX,'O_c');
    axis equal;rotate3d on;
    grid on;
    view(0,-90);

    axis vis3d;
    axis tight;
    axis equal;
end