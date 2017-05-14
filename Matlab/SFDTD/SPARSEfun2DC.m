function [idx] = SPARSEfun2DC(p, thresholddB, p0)
%Transform threshold to Pa
threshold = p0 * 10^(thresholddB/20);
%Calculate blurring filter
PSF = fspecial('gaussian',7,10);
%Create copy of p that is normalised by threshold
temp = abs(p) ./ threshold;
%Set the 'quiet' regions to 0
temp2 = floor(temp);
%Reduce values above 1 to 1, so that blurring isnt too strong in some areas
temp2(temp2 > 1) = 1;
%Convolve 2d matrix with gaussian blurring filter to smooth
idx = imfilter(temp2,PSF,'symmetric','conv');
end
