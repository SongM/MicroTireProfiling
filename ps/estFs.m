% Estimate fs
% by Boren Li
% 2016.09.03
function fs = estFs(imIrrad, invalidMap, dRGB, s_RGB)
% Inputs:
% 1) imIrrad: normalized image irradiance, NUM_Y x NUM_X x 3 x NUM_IMAGES
% 2) maskRegion: estimate pixel locations only within the mask (==1), NUM_Y x NUM_X
% 3) invalidMap: binary map to save invalid pixels (saturated OR
% shadowed), NUM_Y x NUM_X x NUM_IMAGES, == 1 means invalid
% 4) dRGB: unit diffuse color estimate, NUM_Y x NUM_X x 3
% 5) s_RGB: unit specular color, 3 x 1
% Output:
% 1) fs: specular geometrical scaling factor, NUM_Y x NUM_X x NUM_IMAGES
fprintf('Estimating fs...');

[NUM_Y, NUM_X, ~, NUM_IMAGES] = size(imIrrad);

%% Reshape and permute the matrices
imIrrad = reshape(imIrrad, [NUM_Y*NUM_X, 3, NUM_IMAGES]);
imIrrad = permute(imIrrad, [3,2,1]);    % NUM_IMAGES x 3 x NUM_Y*NUM_X
invalidMap = reshape(invalidMap, [NUM_Y*NUM_X, NUM_IMAGES])';   % NUM_IMAGES x NUM_Y*NUM_X
dRGB = reshape(dRGB, [NUM_Y*NUM_X, 3])';    % 3 x NUM_Y*NUM_X

%% Find pixel indices within the mask
% pixInd = find(maskRegion == 1); % find pixels within the mask
% NUM_PIX = size(pixInd, 1);  % number of pixels to process
NUM_PIX = NUM_X*NUM_Y;
%% Start the main loop
fs = NaN(NUM_IMAGES, NUM_Y*NUM_X);   % initialize
for ii = 1:NUM_PIX
    % Get normalized image irradiances in the RGB space
    tmp_imIrrad = imIrrad(:,:,ii);  % NUM_IMAGES x 3
    % Get invalid map
    tmp_invalidMap = invalidMap(:,ii);    % NUM_IMAGES x 1, shadows and saturations
    tmp_imIrrad = tmp_imIrrad(~tmp_invalidMap,:);
    % Get unit diffuse color
    tmp_dRGB = dRGB(:, ii);
    % Compute fs
    tmp_fs = (tmp_imIrrad * s_RGB - (tmp_imIrrad * tmp_dRGB) * (tmp_dRGB' * s_RGB)) / (1 - (tmp_dRGB' * s_RGB)^2);
    fs(~tmp_invalidMap, ii) = tmp_fs;
end

%% Reshape
fs = reshape(fs', [NUM_Y, NUM_X, NUM_IMAGES]);

fprintf('done!\n');
end
    
    
    
    
    
    
    
    