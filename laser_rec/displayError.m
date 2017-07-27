function error = displayError(resid,b_display)
    if ~exist('b_display','var')
        b_display = 1 ;
    end
    error.resid = resid;
    error.std = std(abs(resid));
    error.max = max(abs(resid));
    if (b_display)
        fprintf(['std=',num2str(error.std),', max=', num2str(error.max),'...']);
    end

end

