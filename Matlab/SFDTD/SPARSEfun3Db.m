function [idx] = SPARSEfun3Db(p, thresholddB)

threshold =  2*10^-5 * 10^(thresholddB/20);

idx = zeros(size(p)+2);
idx(2:end-1,2:end-1,2:end-1) = p;
idx(idx > threshold) = 1;
idx = smooth3(idx);
idx = find(idx(2:end-1,2:end-1,2:end-1) > threshold);
% idx(idx > threshold) = 1;
% idx = floor(idx);
end