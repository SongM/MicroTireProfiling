function ps = dropProcessingSizeForPS(ps_raw, DROP_RATIO)
    ps = ps_raw;
    ps.DROP_RATIO = DROP_RATIO;
    if (DROP_RATIO == 1)
        fprintf('original size');
    else
        fprintf(['DROP_RATIO = ', num2str(DROP_RATIO), ', ']);
        ps.region = [ceil( (ps_raw.region(1:2)-1)/DROP_RATIO ) + 1, floor( (ps_raw.region(3:4)-1)/DROP_RATIO ) + 1];
        shift_due_to_size_drop = (ps.region-1)*DROP_RATIO + 1 - ps_raw.region;
        start_point = shift_due_to_size_drop + 1;
        ps.region_due_to_size_drop = ps_raw.region_due_to_size_drop + shift_due_to_size_drop;
        
        
        x_ind = start_point(1):DROP_RATIO:size(ps_raw.im_raw,1);
        y_ind = start_point(2):DROP_RATIO:size(ps_raw.im_raw,2);

        ps.im_raw = ps_raw.im_raw(x_ind, y_ind, :, :);
%         ps.lightStrength = ps_raw.lightStrength(x_ind, y_ind, :, :);
        ps.xC = ps_raw.xC(x_ind, y_ind);
        ps.yC = ps_raw.yC(x_ind, y_ind);
        ps.zC_tire_original_est = ps_raw.zC_tire_original_est(x_ind, y_ind);
%         ps.zC_whiteboard = ps_raw.zC_whiteboard(x_ind, y_ind);
        disp('drop size down');
    end
end