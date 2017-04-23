function [idx] = SPARSEfun3Db(p, thresholddB)

threshold =  2*10^-5 * 10^(thresholddB/20);

idx = p;
idx(idx > threshold) = 1;
idx = smooth3(idx);
idx = find(idx > threshold);
% idx(idx > threshold) = 1;
% idx = floor(idx);
end