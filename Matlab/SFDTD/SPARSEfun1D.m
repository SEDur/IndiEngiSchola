function [idx] = SPARSEfun1D(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/10);

p(end+1) = 0;

temp = abs(p) ./ threshold;

temp = floor(temp);

temp(temp > 1) = 1;

temp2 = decimate(temp, 2);

temp2(temp2 > threshold) = 1.0;
temp3 = round(temp2);
temp4 = interp(temp3,2);
% plot(temp4);
% drawnow;
idx = temp4;
end