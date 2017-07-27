function drawSphere(c, R, n_surface)
    if (~exist('n_surface','var')),  n_surface = 10;  end
    [x,y,z] = sphere(n_surface);
    hold on; surf(x*R + c(1), y*R + c(2), z*R + c(3));
end



