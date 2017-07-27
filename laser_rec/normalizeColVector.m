function [normalizedV,length] = normalizeColVector(V)
    length = sqrt(sum(V.^2));
    normalizedV = V./repmat(length,3,1);

end