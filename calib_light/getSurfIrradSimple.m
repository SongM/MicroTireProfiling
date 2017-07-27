%%%%%%%%%%%%%%%%%%%%%%%%%
% Get surface irradiance (using vectors in {C})
% without cast shadow rendering (simple version)
% by Boren Li
% 2016.04.13
%%%%%%%%%%%%%%%%%%%%%%%%%
function [surfIrrad, shadowMap] = getSurfIrradSimple(nC, lC, flagSelfShadow)
% Inputs:
% 1) nC: unit surface normal in {C}, 3 x (NUM_Y*NUM_X)
% 2) lC: unit light direction in {C}, 3 x (NUM_Y*NUM_X) x NUM_IMAGE
% 3) flagSelfShadow: == 1 implies including self-shadow
% Outputs:
% 1) surfIrrad: surface irradiance, 1 x (NUM_Y*NUM_X) x NUM_IMAGE
% 2) shadowMap: .selfShadowMap, 1 x (NUM_Y*NUM_X) x NUM_IMAGE, 1 implies
% shadowed surface point

fprintf('Computing surface irradiance term (Lambert irradiance law) ...');
[~, ~, NUM_IMAGES] = size(lC);
%% Lambert irradiance term
surfIrrad = dot(repmat(nC, [1,1,NUM_IMAGES]), lC, 1);

%% Include self shadows
if flagSelfShadow == 1
    shadowMap.selfShadowMap = surfIrrad < 0;
    surfIrrad(shadowMap.selfShadowMap) = 0;
else
    shadowMap.selfShadowMap = false(size(surfIrrad));
end

fprintf('done!\n');

end
