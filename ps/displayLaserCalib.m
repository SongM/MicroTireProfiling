function displayLaserCalib(CamCalib,LaserCalib)

color_map = hsv(numel(LaserCalib.c_calib_points)+1);
board_proc = 1:numel(LaserCalib.c_calib_points);

displayCamCalib(CamCalib,[],1,101);
hold on;
for kk = board_proc
    X_c_laser = LaserCalib.c_calib_points{kk};
    scatter3(X_c_laser(1,:),X_c_laser(3,:),-X_c_laser(2,:),2,color_map(kk,:));
end

plot3(LaserCalib.plane_display(1,:),LaserCalib.plane_display(3,:),-LaserCalib.plane_display(2,:),'Color',color_map(end,:));
end