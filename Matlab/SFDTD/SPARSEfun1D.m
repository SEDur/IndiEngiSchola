function [idx] = SPARSEfun1D(p, thresholddB)
%Creates an indexing array that can be used to limit the number of
%cells to solve

%Convert threshold to pressure value
threshold = 10^-12 * 10^(thresholddB/10);
%Append a 0 to accomodate for later bin loss due to decimation and
%interpolation
p(end+1) = 0;
%Normalise domain by threshold
temp = abs(p) ./ threshold;
%Floor noise to 0
temp = floor(temp);
%round big numbers to 1 [COULD BE DEPRECATED?]
temp(temp > 1) = 1;
%decimate the numbers of points 
temp2 = decimate(temp, 2);
%round up the smoothed points
temp2(temp2 > 0.1) = 1.0;
%supress noise
temp3 = floor(temp2);
%Interpolate to return to correct number of points
temp4 = interp(temp3,2);
%round up the smoothed points
temp2(temp2 > 0.1) = 1.0;
%supress noise
temp3 = floor(temp2);
idx = temp4;
end