
n_act = length(active_images);
ind_active = find(active_images);
if exist('center_optim'),
    center_optim = double(center_optim);
end;
if exist('est_alpha'),
    est_alpha = double(est_alpha);
end;
if exist('est_dist'),
    est_dist = double(est_dist);
end;
if exist('est_fc'),
    est_fc = double(est_fc);
end;
if exist('est_aspect_ratio'),
    est_aspect_ratio = double(est_aspect_ratio);
end;