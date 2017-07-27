function balls_im_struct = getRawBallIm(param,ball)

    ball_size_range = param.caliblightpos.ball_size_range;
    resize_scale = param.caliblightpos.resize_scale;
    ball_range_scale = param.caliblightpos.ball_range_scale;
    brightness_th = param.caliblightpos.brightness_th;
    ratio_radius_brightness_spot = param.caliblightpos.ratio_radius_brightness_spot;
    N_lights = param.N_lights;


    temp_im = 0;
    ims = [];
    for im_i = 1:N_lights
        ims(:,:,im_i) = rgb2gray( getImFromPath( ball.bright_spot_im_filepath{im_i} ) );
        temp_im = temp_im + ims(:,:,im_i);
    end
    temp_im = temp_im/N_lights;
end