

function new_p = rotateAroundAxis(p,axis_p,axis_v,theta)
    x=p(1,:);
    y=p(2,:);
    z=p(3,:);
    a=axis_p(1,:);
    b=axis_p(2,:);
    c=axis_p(3,:);
    u=axis_v(1,:);
    v=axis_v(2,:);
    w=axis_v(3,:);
    
    new_x = ( a.*(v.^2+w.^2)-u.*(b.*v+c.*w-u.*x-v.*y-w.*z)).*(1-cos(theta) ) + x.*cos(theta) + (-c.*v+b.*w-w.*y+v.*z).*sin(theta);
    new_y = ( b.*(u.^2+w.^2)-v.*(a.*u+c.*w-u.*x-v.*y-w.*z)).*(1-cos(theta) ) + y.*cos(theta) + ( c.*u-a.*w+w.*x-u.*z).*sin(theta);
    new_z = ( c.*(u.^2+v.^2)-w.*(a.*u+b.*v-u.*x-v.*y-w.*z)).*(1-cos(theta) ) + z.*cos(theta) + (-b.*u+a.*v-v.*x+u.*y).*sin(theta);
    new_p(1,:) = new_x;
    new_p(2,:) = new_y;
    new_p(3,:) = new_z;
    
end



