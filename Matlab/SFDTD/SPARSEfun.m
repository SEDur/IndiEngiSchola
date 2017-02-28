function [pidxRow, pidxCol, uxidx, uyidx] = SPARSEfun(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/20);

[pidxRow pidxCol] = find(abs(p)>threshold);

uxidx = 0;
uyidx = 0;

if length(pidxCol) > 1
surf(p(pidxRow, pidxCol))
view(2);
end
end