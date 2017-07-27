function_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\functions';
addpath(genpath(function_folder));
rmpath(genpath([function_folder,'\','New Folder']));
getStarted;

case_name = '20170724';
data_folder = 'C:\Users\s\Google Drive\projects\honda_20170118\laser_photo';

mat_file = checkCase(data_folder, case_name);load(mat_file);
load(mat_file);

%%



load(individual_mat_file.caliblightpos);
% load(mat_file);


% old_balls = balls;
% old_balls_im = balls_im;
task = ['2.C.extra.manuallyCorrectLightPos'];
displayTask(task,3);
figure_count = figure_count + 1;

prompt={'b_ck_x','b_ck_y','b_radius','b_sik_x','b_sik_y','b_sik_size','color_scale'};
dlg_title = 'Input';

color_scale = 2;
save(individual_mat_file.caliblightpos,'-v7.3');

% ind_ball_pool = 3:15;
ind_ball_pool= 7;
for ind_ball_i = 1:numel(ind_ball_pool) 
    ball_k = ind_ball_pool(ind_ball_i);
%     load(individual_mat_file.caliblightpos);
    ball = balls(ball_k);
    ball_im = balls_im{ball_k};
    ball_area = ball.b_ball_area;

    b_ck = ball.b_ck - ball_area(1:2) + 1;
    b_radius = ball.b_radius;

    for light_i = 1:param.N_lights

        ball_im_ind = ball_im.bright_spot_im_whole(:,:,light_i);
        b_sik_size_ind = ball_im.bright_spot_size(1,light_i);

        ball_im_ind_disp = ball_im_ind*color_scale;
        ball_im_ind_disp(ball_im_ind_disp>1) = 1;

        b_sik_ind = ball.b_sik(light_i,:) - ball_area(1:2) + 1;

        figure(figure_count);clf;set(figure_count,'OuterPosition',[10,200,800,600]);
        subplot(1,2,1);
        imshow(ball_im_ind_disp);
        viscircles(b_ck, b_radius, 'Color', 'r' ,'LineWidth', 1);
        viscircles(b_sik_ind, b_sik_size_ind, 'Color', 'g');
        title({['ball#',num2str(ball_k),', light#', num2str(light_i)],...
            [num2str(b_ck),', r=',num2str(b_radius)],...
            [num2str(b_sik_ind), ', r=',num2str(b_sik_size_ind)]});
        axis on;

        answer = 0;

        while(numel(answer)>0)
            defaultans = {num2str(b_ck(1)), num2str(b_ck(2)), num2str(b_radius),...
                num2str(b_sik_ind(1)), num2str(b_sik_ind(2)), num2str(b_sik_size_ind), ...
                num2str(color_scale)};

            answer = inputdlg(prompt,dlg_title,1,defaultans);
            if (numel(answer)>0)
                old_b_ck = b_ck;
                old_b_radius = b_radius;
                old_b_sik_ind = b_sik_ind;
                old_b_sik_size_ind = b_sik_size_ind;
                b_ck = [str2num(answer{1}),str2num(answer{2})];
                b_radius = str2num(answer{3});
                b_sik_ind = [str2num(answer{4}),str2num(answer{5})];
                b_sik_size_ind = str2num(answer{6});

                color_scale = str2num(answer{7});
                ball_im_ind_disp = ball_im_ind*color_scale;
                ball_im_ind_disp(ball_im_ind_disp>1) = 1;


                figure(figure_count);clf;set(figure_count,'OuterPosition',[10,200,800,600]);
                subplot(1,2,2);
                imshow(ball_im_ind_disp);
                viscircles(old_b_ck, old_b_radius, 'Color', 'r' ,'LineWidth', 1);
                viscircles(old_b_sik_ind, old_b_sik_size_ind, 'Color', 'g');
                title({'previous position', [num2str(old_b_ck),', r=',num2str(old_b_radius)],...
                    [num2str(old_b_sik_ind), ', r=',num2str(old_b_sik_size_ind)]});
                axis on;
                subplot(1,2,1);
                imshow(ball_im_ind_disp);
                viscircles(b_ck, b_radius, 'Color', 'r' ,'LineWidth', 1);
                viscircles(b_sik_ind, b_sik_size_ind, 'Color', 'g');
                title({['ball#',num2str(ball_k),', light#', num2str(light_i)],...
                    [num2str(b_ck),', r=',num2str(b_radius)], [num2str(b_sik_ind), ', r=',num2str(b_sik_size_ind)]});
                axis on;
            end
        end


        balls(ball_k).b_ck = b_ck + ball_area(1:2) - 1;
        balls(ball_k).b_radius = b_radius;
        balls(ball_k).b_sik(light_i,:) = b_sik_ind + ball_area(1:2) - 1;

        balls_im{ball_k}.bright_spot_size(1,light_i) = b_sik_size_ind;




    end
    
    clear ball_k ball ball_im ball_area b_ck b_radius
    clear light_i ball_im_ind b_sik_size_ind ball_im_ind_disp
    clear answer b_sik_ind defaultans old_b_ck old_b_radius old_b_sik_ind old_b_sik_size_ind

    save(individual_mat_file.caliblightpos);


end
    clear prompt dlg_title ball_k color_scale ind_ball_pool ind_ball_i

save(individual_mat_file.caliblightpos)


job_done = [job_done; task];
clear task;



task = '2.C.4.getLightPosFromCkAndSik';
displayTask(task,3);


balls_center_plane_eqn_3_variable = balls(1).ball_center_plane_eqn;

for ball_k = 1:numel(balls)

    ball = balls(ball_k);
    b_ck = ball.b_ck;
    b_sik = ball.b_sik;

    C_ck = get3DCoorFromPlaneEqnAnd2DCoor...
        (balls_center_plane_eqn_3_variable, b_ck', CamCalib);

    
    R = param.caliblightpos.ball_radius;
    fc = CamCalib.fc;
    cc = CamCalib.cc;
    kc = CamCalib.kc;
    alpha_c = CamCalib.alpha_c;

    i_x_undistort = normalize_pixel(b_sik', fc, cc, kc, alpha_c);
    l_s = [i_x_undistort;ones(1,size(i_x_undistort,2))];
    alpha = sum(l_s.^2);
    beta = C_ck'*l_s;
    gamma = sum(C_ck.^2)-R^2;
    C_sik_Z = (beta - sqrt(beta.^2-gamma*alpha))./alpha;
    C_sik = l_s.*repmat(C_sik_Z,3,1);

        
    RadiusFromImage = ball.b_radius/mean(CamCalib.fc)*C_ck(3);
    fprintf(['b_ck = [', num2str(b_ck), '], C_ck = [', num2str(C_ck'),...
        '], radius = ', num2str(RadiusFromImage),...
        ' (default_radius = ',num2str(param.caliblightpos.ball_radius),').\n']);
        
    balls(ball_k).C_ck = C_ck;
    balls(ball_k).C_sik = C_sik;
    balls(ball_k).C_Radius_FromImage = RadiusFromImage;



end



[balls, LightPos, error] = getLightPosFromCkAndSik(balls, error);

displayLightPosAndBalls(LightPos, balls, ref_board_balls, CamCalib , figure_count );



clear balls_center_plane_eqn_3_variable ball_k ball b_ck b_sik C_ck R fc cc kc alpha_c
clear i_x_undistort l_s alpha beta gamma C_sik_Z C_sik RadiusFromImage
save(individual_mat_file.caliblightpos);

clear old_balls old_balls_im balls balls_im balls_center_plane_eqn
job_done = [job_done;task];
clear task;
saveMatFile();

