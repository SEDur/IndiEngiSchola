function [p, ux, uy, uz] = SFDTD3Dfun(p, pCx, pCy, pCz, ux, uy, uz, uCx, uCy, uCz, Rx, Ry, Rz, ZL, ZR, ZT, ZB, ZF, ZG, idx)

if(size(idx(idx > 0),1) > 10)
    %  mat(col, row)
    for i = 1 : size(ux,1)-1
        for i1 = 2 : size(ux,2)-2
            for i2 = 1 : size(ux,3)-1
                if(idx(i, i1, i2) > 0)
                    ux(i,i1,i2) = ux(i,i1,i2) - uCx*(p(i,i1,i2)-p(i,i1-1,i2));
                end
            end
        end
    end
    
    for i = 2 : size(uy,1)-2
        for i1 = 1 : size(uy,2)-1
            for i2 = 1 : size(uy,3)-1
                if(idx(i, i1, i2) > 0)
                    uy(i,i1,i2) = uy(i,i1,i2) - uCx*(p(i,i1, i2)-p(i-1,i1, i2));
                end
            end
        end
    end
    for i = 1 : size(uz,1)-2
        for i1 = 1 : size(uz,2)-1
            for i2 = 2 : size(uz,3)-1
                if(idx(i, i1, i2) > 0)
                    uz(i,i1,i2) = uz(i,i1,i2) - uCx*(p(i,i1, i2)-p(i,i1, i2-1));
                end
            end
        end
    end
else
    % update the non-boundary condition nodes for velocity
    ux(:, 2:end-1, :) = ux(:, 2:end-1,:) - uCx*(p(:, 2:end,:) - p(:, 1:end-1, :));
    uy(2:end-1, :, :) = uy(2:end-1, :, :) - uCy*(p(2:end, :, :) - p(1:end-1, :, :));
    uz(:, :, 2:end-1) = uz(:, :, 2:end-1) - uCz*(p(:, :, 2:end) - p(:, :, 1:end-1));
end

% update the velocity at the right wall
ux(:, end, :) = ((Rx - ZR)/(Rx + ZR))*ux(:, end, :) ...
    + (2/(Rx + ZR))*p(:, end, :);

%update the velocity at the left wall
ux(:, 1, :) = ((Rx - ZL)/(Rx + ZL))*ux(:, 1, :) - (2/(Rx + ZL))*p(:, 1, :);

%update the velocity at the top wall
uy(end, :, :) = ((Ry - ZF)/(Ry + ZF))*uy(end, :, :) ...
    + (2/(Ry + ZT))*p(end, :, :);

%update the velocity at the bottom wall
uy(1, :, :) = ((Ry - ZB)/(Ry + ZB))*uy(1, :, :) - (2/(Ry + ZB))*p(1, :, :);

%update the velocity at the ceiling
uz(:, :, end) = ((Rz - ZT)/(Rz + ZT))*uz(:, :, end) + ...
    (2/(Rz + ZT))*p(:, :, end);

%update the velocity at the floor
uz(:, :, 1) = ((Rz - ZG)/(Rz + ZG))*uz(:, :, 1) - ...
    (2/(Rz + ZG))*p(:, :, 1);

if(size(idx(idx > 0),1) > 10)
    for i = 1 : size(p,1)
        for i1 = 1 : size(p,2)
            for i2 = 1 : size(p,3)
                if(idx(i, i1, i2) > 0)
                    p(i,i1,i2) = p(i,i1,i2) - pCx*(ux(i, i1 + 1, i2) - ux(i, i1, i2))...
                        - pCy*(uy(i+1, i1, i2) - uy(i, i1, i2))...
                        - pCz*(uz(i,i1, i2 + 1) - uz(i, i1, i2));
                end
            end
        end
    end
else
    p = p - pCx*(ux(:, 2:end, :) - ux(:, 1:end-1, :))...
        - pCy*(uy(2:end, :, :) - uy(1:end-1, :, :))...
        - pCz*(uz(:, :, 2:end) - uz(:, :, 1:end-1));
end
end