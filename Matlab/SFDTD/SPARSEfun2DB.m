function [idx] = SPARSEfun2DB(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/10);

p(end+1,1:end) = 0;
p(1:end,end+1) = 0;

for i = 1 : size(p, 1)
temp(i,:) = decimate(p(i,:), 2);
end
for i = 1 : size(temp, 2)
temp2(:,i) = decimate(temp(:,i), 2);
end

temp3 = abs(temp2) ./ threshold;

temp3 = floor(temp3);

temp3(temp3 > 1) = 1;

% temp3(temp3 > 0.1) = 1.0;

% temp3 = floor(temp3);

temp4 = ceil(interp2(temp3));
% figure(2);
% mesh(temp4);
% drawnow;
idx = temp4(1:end-1, 1:end-1);
end