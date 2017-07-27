function calculatedZ = getZfromPQ(p_measured,q_measured,NUM_X,NUM_Y,delta,h1,h2,h3,Z_measured,Z_measured_ind)

% [NUM_Y,NUM_X] = size(zC);
% 
% NUM_X = 5;
% NUM_Y = 3;
% delta = 0.1;

% Z_extra = rand(NUM_Y+1,NUM_X+1);
% Z_extra = Z_extra - Z_extra(1);
% Z = Z_extra(1:NUM_Y,1:NUM_X);
% 
% 
% p = (Z_extra(1:NUM_Y,2:(NUM_X+1))-Z_extra(1:NUM_Y,1:NUM_X)) / delta;
% q = (Z_extra(2:(NUM_Y+1),1:NUM_X)-Z_extra(1:NUM_Y,1:NUM_X)) / delta;
% p_measured = p;
% q_measured = q;


p_measured = inpaint_nans(p_measured,3);
q_measured = inpaint_nans(q_measured,3);
if (nargin<5)
    fprintf(2,'not enough input');
    return;
end
if (nargin==5)
    h1 = 0.1;
    h2 = 1;
    h3 = h2;
    Z_measured_ind = zeros(NUM_Y,NUM_X);
    Z_measured_ind(1,1) = 1;

    Z_measured = zeros(NUM_Y,NUM_X);
    Z_measured(1,1) = 0;
end

D_i_minus1 = -ones(NUM_Y,NUM_X) * h3;
D_i_minus1(1,:) = 0;
D_i_plus1 = -ones(NUM_Y,NUM_X) * h3;
D_i_plus1(end,:) = 0;
D_j_minus1 = -ones(NUM_Y,NUM_X) * h2;
D_j_minus1(:,1) = 0;
D_j_plus1 = -ones(NUM_Y,NUM_X) * h2;
D_j_plus1(:,end) = 0;
D_main = -(D_i_minus1+D_i_plus1+D_j_minus1+D_j_plus1) + delta^2*h1*Z_measured_ind;

D_all = [D_main(:), D_i_minus1(:), D_i_plus1(:), D_j_minus1(:), D_j_plus1(:)];
D_pos = [0, 1, -1, NUM_Y, -NUM_Y];
M=spdiags(D_all, D_pos, NUM_X*NUM_Y, NUM_X*NUM_Y);
% A=full(M);




p_j_minus1 = [zeros(NUM_Y,1), p_measured];
p_j_minus1(:,end) = [];
p_ij = p_measured;
p_ij(:,end) = 0;
q_i_minus1 = [zeros(1,NUM_X); q_measured];
q_i_minus1(end,:) = [];
q_ij = q_measured;
q_ij(end,:) = 0;
v = h2*p_ij - h2*p_j_minus1 + h3*q_ij - h3*q_i_minus1 - delta*h1*Z_measured_ind.*Z_measured;
V = -delta * v(:);
temp_Z = M \ V;
calculatedZ = reshape(temp_Z,NUM_Y,NUM_X);




end