function [idx] = SPARSEfun(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/10);

p(end+1,1:end) = 0;
p(1:end,end+1) = 0;

temp = abs(p) ./ threshold;

temp = floor(temp);

temp(temp > 1) = 1;

for i = 1 : size(temp, 1)
temp2(i,:) = decimate(temp(i,:), 2);
end
for i = 1 : size(temp2, 2)
temp3(:,i) = decimate(temp2(:,i), 2);
end


temp3(temp3 > 0.1) = 1.0;
temp3 = round(temp3);

% temp4 = interp2(temp3);
temp4 = interp2(temp3);
figure(2);
mesh(temp4);
idx = temp4(1:end-1, 1:end-1);
end