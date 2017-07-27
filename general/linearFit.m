function [params, err] = linearFit(data)
    % assume data is n data with dimention of m X'=(x1,x2,...,xm) 
    % data is in col vector, data size = m*n 
    % then the linear fitting eqn can be write as p(1)*X(1) + 
    % p(2)*X(2) + ... + p(m-1)*X(m-1) + p(m)*X(m) = 1;
    % or P'*X = X'*P = 1
    % with n data => data'*P = ones(n,1)
    % params is an m*1-dimension col vector.
    
    A = data';
    b = ones(size(data,2),1);
    params = A\b;
    resid = (A*params - b)/sqrt(sum(params.^2));
    err.resid = resid;
%     err.mse = sqrt(mean(resid.^2));
    err.std = std(resid);
    err.max = max(abs(resid));
%     err.min = min(abs(resid));
end
