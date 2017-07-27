function display_surface_normal_map()
step = 0.001;
x = -1.2:step:1.2;
y = -1.2:step:1.2;
[xx,yy] = meshgrid(x,y);

R = 1;
temp = R^2 - xx.^2 - yy.^2;
temp(temp<0) = 0;
zz = sqrt(temp);


% figure();
% surf(xx,yy,zz);


p = zz(2:end,:) - zz(1:(end-1),:);
p(end+1,:) = 0;
q = zz(:,2:end) - zz(:,1:(end-1));
q(:,end+1) = 0;
nC = [];
nC(:,:,1) = p./sqrt(p.^2+q.^2+step.^2);
nC(:,:,2) = q./sqrt(p.^2+q.^2+step.^2);
nC(:,:,3) = step./sqrt(p.^2+q.^2+step.^2);
% figure();
% subplot(1,2,1);
% imshow((-nC+1)/2)
ind = (zz==0);
nC2 = zeros(size(zz,1),size(zz,2),3);
xx(ind) = 0;
yy(ind) = 0;
zz(ind) = 1;
nC2(:,:,1) = -xx/R;
nC2(:,:,2) = -yy/R;
nC2(:,:,3) = -zz/R;
% subplot(1,2,2);
imshow((-nC2+1)/2)

end