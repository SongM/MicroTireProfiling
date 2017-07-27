function [d_theta, V, C0, figure_count, Z1, Z2, l1, l2] = getRotationAxisAndAngle(C0_est, V_est, Z1_est, Z2_est, Laser_X_c2, l1, l2, N_loop, b_debug, figure_count)
    if (~exist('b_debug','var')),  b_debug = 0;  end
    if (~exist('figure_count','var')),  figure_count = 1;  end
    if (~exist('N_loop','var')),  N_loop = 1;  end

    C0 = C0_est;    V = V_est;  Z1 = Z1_est;    Z2 = Z2_est;

    for loop_i = 1:N_loop
        % d_theta
        d_theta = mean( (l2(1,:).*Z2-l1(1,:).*Z1) ./ ((V(2)-V(3)*l1(2,:)).*Z1 + (V(3)*C0(2) - V(2)*C0(3))) );
        % Z1 & Z2
        for ii=1:size(l1,2)
            l1_ind = l1(:,ii);  l2_ind = l2(:,ii);
            A = [l1_ind+cross(V,l1_ind)*d_theta , -l2_ind];
            b = cross(V,C0)*d_theta;
            res = inv(A'*A)*A'*b;            
            Z1(ii) = res(1);    Z2(ii) = res(2);
        end

        ind = (Z1>max(Laser_X_c2(3,:))+20) | Z1<(min(Laser_X_c2(3,:))-20);
        l1(:,ind) = []; l2(:,ind) = []; Z1(ind) = [];   Z2(ind) = [];
        if(b_debug)
            figure_count = figure_count + 1;    figure(figure_count);
            plot3([l1(1,:).*Z1; l2(1,:).*Z2],  [l1(2,:).*Z1; l2(2,:).*Z2]  ,[Z1; Z2],'-b');hold on;
            scatter3(l1(1,:).*Z1,l1(2,:).*Z1,Z1, 5, 'og');
            scatter3(l2(1,:).*Z2,l2(2,:).*Z2,Z2, 5, '+r');
            scatter3(Laser_X_c2(1,:),Laser_X_c2(2,:),Laser_X_c2(3,:),5,'r','filled');
            
            scatter3(C0(1),C0(2),C0(3),10,'b','filled');
%             Axis = [C0-200*V,C0+200*V];
%             plot3(Axis(1,:),Axis(2,:),Axis(3,:),'b');
            quiver3(C0(1)+V(1)*150,C0(2)+V(2)*150,C0(3)+V(3)*150,-V(1)*300,-V(2)*300,-V(3)*300,'b');
            t_ind = ceil(size(l1,2)/2);
            P1 = l1(:,t_ind)*Z1(t_ind);
            P2 = l2(:,t_ind)*Z2(t_ind);
            H = (P1-C0)'*V;
            P1_new = C0 + H*V;
            plot3([P1_new(1),P1(1)],[P1_new(2),P1(2)],[P1_new(3),P1(3)],'b');
            plot3([P1_new(1),P2(1)],[P1_new(2),P2(2)],[P1_new(3),P2(3)],'b');
 

            title({[num2str(loop_i),',',num2str(d_theta),','],num2str(C0'),num2str(V')});
            view(0,0);  axis equal
        end

        % V & C0
        % (V x l1)*Z1 + (-V x C0) = (l2*Z2 -l1*Z1)/d_theta;
        A = []; b = [];
        for ii=1:size(l1,2)
            l1_ind = l1(:,ii);  l2_ind = l2(:,ii);  
            Z1_ind = Z1(ii);    Z2_ind = Z2(ii);
            temp_A = [];
            temp_A(1,:) = [0,           Z1_ind,     -l1_ind(2,:)*Z1_ind,    1, 0, 0];
            temp_A(2,:) = [-Z1_ind,     0,          l1_ind(1,:)*Z1_ind,     0, 1, 0];
            temp_A(3,:) = [l1_ind(2,:)*Z1_ind, -l1_ind(1,:)*Z1_ind, 0,      0, 0, 1];
            A = [A;temp_A];
            b = [b;(l2_ind*Z2_ind - l1_ind*Z1_ind)/d_theta];
        end
%             A(2:3:end,:) = [];
%             b(2:3:end) = [];
        res = inv(A'*A)*A'*b;            
        % V = res(1:3) % direction;
        % -V x C0 = res(4:5) 
        V = res(1:3)/sqrt(sum(res(1:3).^2));
        A = [0,-res(2); -res(3),res(1); res(2),0];
        b = res(4:6);
        res = inv(A'*A)*A'*b;            
        C0 = [res(1);0;res(2)];
    end


end