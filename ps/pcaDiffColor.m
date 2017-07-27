%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the principal unit diffuse color for each pixel
% by Boren Li
% 2016.06.09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dRGB, outlierMap, mean_color_error] = pcaDiffColor(imIrrad, invalidMap, Td)
% Inputs:
% 1) imIrrad: normalized image irradiance, NUM_Y x NUM_X x 3 x NUM_IMAGES
% 2) maskRegion: estimate pixel locations only within the mask (==1), NUM_Y x NUM_X
% 3) invalidMap: binary map to save invalid pixels (saturated OR
% shadowed), NUM_Y x NUM_X x NUM_IMAGES, == 1 means invalid
% 4) Td: outlier rejection threshold for RPCA
% Outputs:
% 1) dRGB: unit diffuse color estimate, NUM_Y x NUM_X x 3
% 2) outlierMap: == 1 means outlier, NUM_Y x NUM_X x NUM_IMAGES


b_debug = 0;
    
fprintf('Estimating dRGB...');

[NUM_Y, NUM_X, ~, NUM_IMAGES] = size(imIrrad);

%% Reshape and permute the matrices
imIrrad = reshape(imIrrad, [NUM_Y*NUM_X, 3, NUM_IMAGES]);
imIrrad = permute(imIrrad, [3,2,1]);
invalidMap = reshape(invalidMap, [NUM_Y*NUM_X, NUM_IMAGES])';

%% Find pixel indices within the mask



%% Start the main loop
dRGB = NaN(NUM_Y*NUM_X, 3);   % initialize
outlierMap = false(NUM_Y*NUM_X, NUM_IMAGES);
mean_color_error = zeros(NUM_Y*NUM_X,1);
NUM_PIX = NUM_Y*NUM_X;
for ii = 1:NUM_PIX
    % Get normalized image irradiances in the RGB space
    
%     fprintf([num2str(ii),'..']);
    tmp_imIrrad = imIrrad(:,:,ii);  % NUM_IMAGES x 3
    
    % Get invalid map
    tmp_mask = invalidMap(:,ii);    % NUM_IMAGES x 1
    
    if ~isempty(find(tmp_mask == 0, 1)) % if there exists at least one valid pixel
        % Initialize
        tmp_flagTerminate = 0;  % termination flag in each iteration
        loop_count = 0;

        if(b_debug)
            fprintf([num2str(ii),', ']);
        end
        
        while tmp_flagTerminate == 0
            loop_count = loop_count + 1;
            ttmp_imIrrad = tmp_imIrrad(~tmp_mask,:); % choose the valid ones
            % Construct L = E^T * E (matrix for principle component analyses)
            ttmp_L = ttmp_imIrrad' * ttmp_imIrrad;
            % PCA
            [ttmp_V, ttmp_D] = eig(ttmp_L);    % eigenvalue
            [~, ttmp_ind] = max(diag(ttmp_D));   % find maximum absolute eigenvalue
            
            ttmp_dRGB = ttmp_V(:,ttmp_ind);   % Choose the corresponding column, remove the sign ambiguity
            
            %%
            %             ttmp_dRGB = abs(ttmp_dRGB);
            %%%%%%%%%%%%%%%%%%%
            % Outlier rejection using  residuals
            %%%%%%%%%%%%%%%%%%%
            % Project results on the principal diffuse color
            ttmp_fd = ttmp_imIrrad * ttmp_dRGB; % M x 1, M is number of valid pixels
            % Compute residuals
            ttmp_R = ttmp_fd * ttmp_dRGB' - ttmp_imIrrad;   % M x 3
            ttmp_r = sqrt(sum(ttmp_R.^2,2));    % M x 1
            %%
            if(b_debug)
                fprintf([num2str(mean(ttmp_r)),'. ']);
            end
            if (min(ttmp_dRGB)*max(ttmp_dRGB)<-1e-10)
                fprintf(['ii=',num2str(ii),', loop=',num2str(loop_count),', dRGB = ', num2str(ttmp_dRGB'),'\n']);
            end
            ttmp_dRGB = abs(ttmp_dRGB);
            %%
            
            
            % Detect outliers
            ttmp_outlier = true(NUM_IMAGES,1);  % initialize as all of the pixels are outliers
            if mean(ttmp_r) <= Td % residual already sufficiently small
                dRGB(ii,:) = ttmp_dRGB;
                outlierMap(ii,:) = tmp_mask';
                tmp_flagTerminate = 1;
                mean_color_error(ii,:) = mean(ttmp_r);
            else
%                 ttmp_r = ttmp_r - mean(ttmp_r); % mean-shift
                [~, ttmp_ind] = max(abs( ttmp_r-mean(ttmp_r) ));  % residual / std, exclude one outlier each time
                ttmp_validSeq = 1:NUM_IMAGES;
                ttmp_validSeq(tmp_mask == 1) = [];  % exclude past outliers
                ttmp_validSeq(ttmp_ind) = [];   % exclude current outlier
                ttmp_outlier(ttmp_validSeq) = 0;
                tmp_mask = ttmp_outlier;    % update
            end
        end
        if(b_debug)
            fprintf('\n');
        end
    end
end

%% Reshape
dRGB = reshape(dRGB, [NUM_Y, NUM_X, 3]);
outlierMap = reshape(outlierMap, [NUM_Y, NUM_X, NUM_IMAGES]);
mean_color_error = reshape(mean_color_error, [NUM_Y, NUM_X]);

fprintf('done!\n');
end
