% clear all;close all
% SAT_PERCENT = 0.0001
function getImFromDNG(folder_name,SAT_PERCENT,b_display,figure_index,b_first_im)
    tic;
    N_figure = 1;
    if (~exist('b_display','var')),  b_display = 0;  end
    if (~exist('figure_index','var')),  figure_index = 1:N_figure;    end
    if ( (numel(figure_index)<N_figure) && (figure_index~=0) ),  
        disp(['N_figure = ', num2str(N_figure)]);
        figure_index = 1:N_figure;    
    end
    
    if (~exist('b_first_im','var')),  b_first_im = 0;  end
    

    disp('working folder = ');
    disp(folder_name);
    
    im_file_list = dir([folder_name,'\*dng']);
    N_image = numel(im_file_list);
    if (b_first_im),    N_image = 1;    b_display = 1;    end

    fprintf([num2str(N_image),' images need to be converted,'...,
        ' SAT_PERCENT = ', num2str(SAT_PERCENT),'\n']);
    
    
    for im_i = 1:N_image
        fileNames{1} = [folder_name,'\',im_file_list(im_i).name];
        fprintf(['working on image ', num2str(im_i),'/',num2str(N_image),', ']);
        %% Read raw images
%         disp('Reading raw images...');
        BAYER = 'rggb'; % Bayer pattern of Nikon D810
        % Read raw images
        rawIm = readRawImage(fileNames, BAYER);
        imraw = rawIm{1};

        max_val = quantile(imraw(:),1-SAT_PERCENT);
        im = imraw/max_val;
        im(im>1) = 1;
    
        imwrite(im, [folder_name,'_',num2str(im_i+1000),'_maxval=',num2str(max_val),'.jpg'], 'jpg');
        save([folder_name, '_',num2str(im_i+1000),'.mat'],'im','SAT_PERCENT','max_val');
        
        time = toc;
        fprintf(['max_val=',num2str(max_val),', done! time remain = ',num2str((N_image-im_i)*time/im_i),'s.\n']);
    
        if (b_display)
            figure(figure_index(1));
            clf;
            imshow(im);
            title(['im_i=',num2str(im_i),', max val=',num2str(max_val),', SAT PERCENT=',num2str(SAT_PERCENT)]);
            drawnow;
        end

    end
    toc;
end
