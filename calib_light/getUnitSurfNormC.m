%%%%%%%%%%%%%%%%%
% Get unit surface normal in {C}
% by Boren Li
% 2016.04.13
%%%%%%%%%%%%%%%%%

function nC = getUnitSurfNormC(xyzC)
% Input:
% 1) xyzC: point cloud, [NUM_Y, NUM_X, 3]
% Output:
% 1) nC: unit surface normal in {C}

[NUM_Y, NUM_X, ~] = size(xyzC);

x = xyzC(:,:,1);
y = xyzC(:,:,2);
z = xyzC(:,:,3);

z_i_jp1 = zeros(size(z));
z_i_jp1(1:NUM_Y, 1:NUM_X-1) = z(1:NUM_Y, 2:NUM_X);
z_ip1_j = zeros(size(z));
z_ip1_j(1:NUM_Y-1, 1:NUM_X) = z(2:NUM_Y, 1:NUM_X);
x_i_jp1 = zeros(size(x));
x_i_jp1(1:NUM_Y, 1:NUM_X-1) = x(1:NUM_Y, 2:NUM_X);
y_ip1_j = zeros(size(y));
y_ip1_j(1:NUM_Y-1, 1:NUM_X) = y(2:NUM_Y, 1:NUM_X);
p = (z_i_jp1 - z) ./ (x_i_jp1 - x);
q = (z_ip1_j - z) ./ (y_ip1_j - y);
p(:,NUM_X) = p(:,NUM_X-1);
q(NUM_Y,:) = q(NUM_Y-1,:);

% Compute unit surface normal, point outwards from the surface
nC = zeros(NUM_Y, NUM_X, 3);
nC(:,:,1) = p ./ sqrt(p.^2 + q.^2 + 1);
nC(:,:,2) = q ./ sqrt(p.^2 + q.^2 + 1);
nC(:,:,3) = -1 ./ sqrt(p.^2 + q.^2 + 1);

end