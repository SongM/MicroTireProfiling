function [im, max_val, SAT_PERCENT] = getImFromPath(im_path)
    if ( strcmp(im_path(end-2:end),'mat') || ...
            strcmp(im_path(end-2:end),'MAT') )
        load(im_path);
        t_im = single(im);
        tt_im(:,:,1) = flipud(t_im(:,:,1)');
        tt_im(:,:,2) = flipud(t_im(:,:,2)');
        tt_im(:,:,3) = flipud(t_im(:,:,3)');
        im = tt_im;
        return;
    end
    if ( strcmp(im_path(end-2:end),'jpg') || ...
            strcmp(im_path(end-2:end),'JPG') )
        im = im2single(imread(im_path));
        
        t_im = single(im);
        tt_im(:,:,1) = flipud(t_im(:,:,1)');
        tt_im(:,:,2) = flipud(t_im(:,:,2)');
        tt_im(:,:,3) = flipud(t_im(:,:,3)');
        im = tt_im;

        max_val = 1;
        SAT_PERCENT = 0;
        return;
    end
    if ( strcmp(im_path(end-2:end),'dng') || ...
            strcmp(im_path(end-2:end),'DNG') )

        bayer = 'rggb'; % Bayer pattern of Nikon D810

        im_tiff = Tiff(im_path, 'r');
        offsets = getTag(im_tiff, 'SubIFD');
        setSubDirectory(im_tiff, offsets(1));
        im_raw = read(im_tiff); % Bayer CFA data
        close(im_tiff);
        meta_info = imfinfo(im_path);
        % Crop to only valid pixels
        x_origin = meta_info.SubIFDs{1}.ActiveArea(2)+1; % +1 due to MATLAB indexing
        width = meta_info.SubIFDs{1}.DefaultCropSize(1);
        y_origin = meta_info.SubIFDs{1}.ActiveArea(1)+1;
        height = meta_info.SubIFDs{1}.DefaultCropSize(2);
        im_raw = double(im_raw(y_origin:y_origin+height-1,x_origin:x_origin+width-1));        
        % Linearizing
        if isfield(meta_info.SubIFDs{1}, 'LinearizationTable')
            ltab=meta_info.SubIFDs{1}.LinearizationTable;
            ltab(1,(size(ltab,2)+1):(size(ltab,2)+100)) = ltab(1,size(ltab,2));
            im_raw = ltab(im_raw+1);
        end
        black = meta_info.SubIFDs{1}.BlackLevel(1);
        saturation = meta_info.SubIFDs{1}.WhiteLevel;
        lin_bayer = (im_raw-black)/(saturation-black);
        lin_bayer = max(0,min(lin_bayer,1));

        % White balancing
        wb_multipliers = (meta_info.AsShotNeutral).^-1;
        wb_multipliers = wb_multipliers/wb_multipliers(2);

        switch bayer
            case 'rggb'
                mask = wbmask(size(lin_bayer,1),size(lin_bayer,2),wb_multipliers, 'rggb');
            case 'gbrg'
                mask = wbmask(size(lin_bayer,1),size(lin_bayer,2),wb_multipliers, 'gbrg');
            case 'grbg'
                mask = wbmask(size(lin_bayer,1),size(lin_bayer,2),wb_multipliers, 'grbg');
            case 'bggr'
                mask = wbmask(size(lin_bayer,1),size(lin_bayer,2),wb_multipliers, 'bggr');
        end
        balanced_bayer = lin_bayer .* mask;
        % Demosaicing
        %     lin_rgb = double(demosaic(uint16(balanced_bayer/max(balanced_bayer(:))*2^16), bayer))/2^16;
        lin_rgb = double(demosaic(uint16(balanced_bayer/max(wb_multipliers)*2^16), bayer))/2^16;
        % Use linear RGB space directly
        im = im2single(max(0,min(lin_rgb,1)));
        
        t_im = single(im);
        tt_im(:,:,1) = flipud(t_im(:,:,1)');
        tt_im(:,:,2) = flipud(t_im(:,:,2)');
        tt_im(:,:,3) = flipud(t_im(:,:,3)');
        im = tt_im;

        max_val = 1;
        SAT_PERCENT = 0;
        return;
    end

    
    
    fprintf(2,'error path format.\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colormask = wbmask(m,n,wbmults,align)
% COLORMASK = wbmask(M,N,WBMULTS,ALIGN)
%
% Makes a white-balance multiplicative mask for an image of size m-by-n
% with RGB while balance multipliers WBMULTS = [R_scale G_scale B_scale].
% ALIGN is string indicating Bayer arrangement: ¡¯rggb¡¯,¡¯gbrg¡¯,¡¯grbg¡¯,¡¯bggr¡¯
colormask = wbmults(2) * ones(m,n); %Initialize to all green values
switch align
    case 'rggb'
        colormask(1:2:end,1:2:end) = wbmults(1); %r
        colormask(2:2:end,2:2:end) = wbmults(3); %b
    case 'bggr'
        colormask(2:2:end,2:2:end) = wbmults(1); %r
        colormask(1:2:end,1:2:end) = wbmults(3); %b
    case 'grbg'
        colormask(1:2:end,2:2:end) = wbmults(1); %r
        colormask(1:2:end,2:2:end) = wbmults(3); %b
    case 'gbrg'
        colormask(2:2:end,1:2:end) = wbmults(1); %r
        colormask(1:2:end,2:2:end) = wbmults(3); %b
end
end

