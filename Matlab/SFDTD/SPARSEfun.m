function [pidxRow, pidxCol, uxidx, uyidx] = SPARSEfun(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/20);

interpP = interp2(p);

[pidxRow, pidxCol] = find(abs(interpP)>threshold);

uxidx = 0;
uyidx = 0;

if length(pidxCol) > 1
surf(p(pidxRow, pidxCol))
view(2);
end
end