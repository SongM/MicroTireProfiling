% Compute unit light direction in {C}
% by Boren Li
% 2015.04.15
function [unitLightC, lightSurfDist] = getUnitLight(lightPosC, XYZC)
% Input:
% 1) lightPosC: light position in {C} in mm
% 2) XYZC: surface points in {C} in mm, NUM_Y x NUM_X x 3
% Output:
% unitLightC: unit light direction in {C}, 3 x NUM_PTS x NUM_IMAGE, NUM_PTS is number of
% lightSurfDist: light-surface distance, 1 x NUM_PTS x NUM_IMAGE
fprintf('Computing unit light direction AND light-surface distance in {C} No. ');
NUM_IMAGE = size(lightPosC,2);   % Number of images = number of lights
[NUM_Y, NUM_X, ~] = size(XYZC);
NUM_PTS = NUM_Y * NUM_X;  % Total number of surface patches on the object
surfPtsC = reshape(XYZC, [NUM_PTS, 3])';    % 3 x NUM_PTS
lightVecC = zeros(3, NUM_PTS, NUM_IMAGE);  % Light vector in {C}, a vector from surface point to the light source
unitLightC = zeros(3, NUM_PTS, NUM_IMAGE);
lightSurfDist = zeros(1, NUM_PTS, NUM_IMAGE);
for k = 1:NUM_IMAGE
    fprintf([num2str(k),'...']);
    % Light vector in {C}, a vector from surface point to the light source
    lightVecC(:,:,k) = repmat(lightPosC(:,k), [1,NUM_PTS]) - surfPtsC;
    % Light-surface distance
    lightSurfDist(:,:,k) = sqrt(lightVecC(1,:,k).^2 + lightVecC(2,:,k).^2 + lightVecC(3,:,k).^2);
    % Unit light direction in {C}
    unitLightC(:,:,k) = lightVecC(:,:,k) ./ repmat(lightSurfDist(:,:,k), [3,1]);
end
fprintf('done!\n');
end
