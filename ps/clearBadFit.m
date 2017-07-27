% Clear bad fittings
% by Boren Li
% 2016.11.20

reconData = reconDataUV;
reconData.residual = nanmean(reconDataUV.residualMap, 3);

fprintf('Clearing bad fittings...');
% Reject boundary pixels
% SE = strel('square', 2);
% maskClear = imerode(maskRegion, SE);
maskClear = ones(NUM_Y,NUM_X);
absRes = abs(reconData.residual);
absRes = nanmean(absRes, 3);
absRes(maskClear(:)) = nan;

% Reject according to residual
badFitThres = quantile(absRes(~isnan(absRes(:))), (1-BAD_FIT_PERCENT));
badFitInd = absRes >= badFitThres;

% Reject according to gradient
nCRecon = reconData.n;
pCRecon = -nCRecon(:,:,1) ./ nCRecon(:,:,3);
qCRecon = -nCRecon(:,:,2) ./ nCRecon(:,:,3);
invalidInd = pCRecon >= tand(MAX_SLOPE) | pCRecon <= tand(-MAX_SLOPE) | qCRecon >= tand(MAX_SLOPE) | qCRecon <= tand(-MAX_SLOPE) | badFitInd;
pCRecon(invalidInd) = nan;
qCRecon(invalidInd) = nan;
nCRecon(:,:,1) = pCRecon ./ sqrt(pCRecon.^2 + qCRecon.^2 + 1);
nCRecon(:,:,2) = qCRecon ./ sqrt(pCRecon.^2 + qCRecon.^2 + 1);
nCRecon(:,:,3) = -1 ./ sqrt(pCRecon.^2 + qCRecon.^2 + 1);
clear pCRecon qCRecon

nCReconClear = reshape(nCRecon, [NUM_Y*NUM_X,3]);
nCReconClear(~maskClear(:),:) = nan;
nCReconClear = reshape(nCReconClear, [NUM_Y, NUM_X, 3]);

[maskRowInd, maskColInd] = find(ones(NUM_Y,NUM_X) == 1);

maskXMin = min(maskColInd);
maskXMax = max(maskColInd);
maskYMin = min(maskRowInd);
maskYMax = max(maskRowInd);
clear maskRowInd maskColInd

% Inpaint nans
pClear = -nCReconClear(:,:,1) ./ nCReconClear(:,:,3);
pClear(maskYMin:maskYMax, maskXMin:maskXMax) = inpaint_nans(pClear(maskYMin:maskYMax, maskXMin:maskXMax), 2);
qClear = -nCReconClear(:,:,2) ./ nCReconClear(:,:,3);
qClear(maskYMin:maskYMax, maskXMin:maskXMax) = inpaint_nans(qClear(maskYMin:maskYMax, maskXMin:maskXMax), 2);
nCReconClear(:,:,1) = pClear ./ sqrt(pClear.^2+qClear.^2+1);
nCReconClear(:,:,2) = qClear ./ sqrt(pClear.^2+qClear.^2+1);
nCReconClear(:,:,3) = -1 ./ sqrt(pClear.^2+qClear.^2+1);
clear pClear qClear

nCReconClear = reshape(nCReconClear, [NUM_Y*NUM_X, 3]);
% nCReconClear(~maskRegion(:),:) = nan;
nCReconClear = reshape(nCReconClear, [NUM_Y, NUM_X, 3]);

kdClear = reconData.rho ./ kappa_z;
kdClear(invalidInd) = nan;
kdClear(maskYMin:maskYMax, maskXMin:maskXMax) = inpaint_nans(kdClear(maskYMin:maskYMax, maskXMin:maskXMax), 2);

fprintf('done!\n');