function [idx] = SPARSEfun(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/20);

temp = p ./ threshold;

temp = floor(temp);

temp(temp > 1) = 1;

for i = 1 : size(temp, 1)
temp2(i,:) = decimate(temp(i,:), 2);
end
for i = 1 : size(temp2, 2)
temp3(:,i) = decimate(temp2(:,i), 2);
end
temp3 = ceil(temp3);

temp4 = interp2(temp3);
figure(2);
mesh(temp4);
idx = temp4;
end