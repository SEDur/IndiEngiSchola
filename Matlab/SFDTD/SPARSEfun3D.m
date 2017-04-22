function [idx] = SPARSEfun3D(p, thresholddB)

threshold =  2*10^-5 * 10^(thresholddB/20);

p(end+1,1:end,1:end) = 0;
p(1:end,end+1,1:end) = 0;
p(1:end,1:end,end+1) = 0;
for i2 = 1 : size(p, 3)
    for i = 1 : size(p, 1)
        temp(i,:,i2) = decimate(p(i,:,i2), 2);
    end
    for i = 1 : size(temp, 2)
        temp2(:,i) = decimate(temp(:,i), 2);
    end
    temp3(:,:,i2) = temp2;
end
for i = 1 : size(temp3, 1)
    for i2 = 1 : size(temp3, 2)
        temp4(i,i2,:) = decimate(temp3(i,i2,:), 2);
    end
end
temp5 = abs(temp4) ./ threshold;

temp5 = floor(temp5);

temp5(temp5 > 1) = 1;

temp6 = ceil(interp3(temp5));

idx = temp6(1:end-1, 1:end, 1:end-1);

end