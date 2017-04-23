function [idx] = SPARSEfun2D(p, thresholddB, p0)
% Convert threshold from dB to Pa
threshold = p0 * 10^(thresholddB/20);
% Pad edge of p with 0s to accomodate truncation
p(end+1,1:end) = 0;
p(1:end,end+1) = 0;

% Decimate matrix to operate on fewer points, and to smooth
% Decimate p in x direction
for i = 1 : size(p, 1)
temp(i,:) = decimate(p(i,:), 2);
end
% Decimate p in y direction
for i = 1 : size(temp, 2)
temp2(:,i) = decimate(temp(:,i), 2);
end
% Normalise array by threshold
temp3 = abs(temp2) ./ threshold;
% Cut out low levels
temp3 = floor(temp3);
% Bring index of interest to 1
temp3(temp3 > 1) = 1;
% Interp to complete smoothing and bring back array scale
temp4 = ceil(interp2(temp3));
% Bring back to size of p
idx = temp4(1:end-1, 1:end-1);
end