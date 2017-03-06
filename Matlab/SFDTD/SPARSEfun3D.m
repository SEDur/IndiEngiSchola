function [idx] = SPARSEfun3D(p, thresholddB)

threshold = 10^-12 * 10^(thresholddB/10);

p(end+1,1:end,1:end) = 0;
p(1:end,end+1,1:end) = 0;
p(1:end,1:end,end+1) = 0;

temp = abs(p) ./ threshold;

temp = floor(temp);

temp(temp > 1) = 1;

for i = 1 : size(temp, 1)
    for i1 = 1 : size(temp,3)
        temp2(i,:,i1) = decimate(temp(i,:,i1), 2);
    end
end
for i = 1 : size(temp2, 2)
    for i1 = 1 : size(temp,3)
        temp3(:,i,i1) = decimate(temp2(:,i,i1), 2);
    end
end
for i = 1 : size(temp3, 1)
    for i1 = 1 : size(temp3,2)
        temp4(i,i1,:) = decimate(temp3(i,i1,:), 2);
    end
end


temp4(temp4 > 0.1) = 1.0;
temp4 = round(temp4);

temp5 = interp3(temp4);
figure(2);
spy(temp5(:,:,ceil(size(temp5,2)/2)));
drawnow;
idx = temp5();
end