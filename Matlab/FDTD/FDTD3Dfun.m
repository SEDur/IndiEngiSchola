function [p, ux, uy, uz] = FDTD3Dfun(p, pCx, pCy, pCz, ux, uy, uz, uCx, uCy, uCz, Rx, Ry, Rz, ZL, ZR, ZF, ZB, ZT, ZG)

    ux(:, 2:end-1, :) = ux(:, 2:end-1,:) - uCx*(p(:, 2:end,:) - p(:, 1:end-1, :));
    uy(2:end-1, :, :) = uy(2:end-1, :, :) - uCy*(p(2:end, :, :) - p(1:end-1, :, :));
    uz(:, :, 2:end-1) = uz(:, :, 2:end-1) - uCz*(p(:, :, 2:end) - p(:, :, 1:end-1));

    % update the velocity at the right wall
    ux(:, end, :) = ((Rx - ZR)/(Rx + ZR))*ux(:, end, :) ...
        + (2/(Rx + ZR))*p(:, end, :);

    %update the velocity at the left wall
    ux(:, 1, :) = ((Rx - ZL)/(Rx + ZL))*ux(:, 1, :)...
        - (2/(Rx + ZL))*p(:, 1, :);

    %update the velocity at the top wall
    uy(end, :, :) = ((Ry - ZF)/(Ry + ZF))*uy(end, :, :) ...
        + (2/(Ry + ZT))*p(end, :, :);

    %update the velocity at the bottom wall
    uy(1, :, :) = ((Ry - ZB)/(Ry + ZB))*uy(1, :, :)...
        - (2/(Ry + ZB))*p(1, :, :);
    
    %update the velocity at the ceiling
    uz(:, :, end) = ((Rz - ZT)/(Rz + ZT))*uz(:, :, end) + ...
    (2/(Rz + ZT))*p(:, :, end);

    %update the velocity at the floor
    uz(:, :, 1) = ((Rz - ZG)/(Rz + ZG))*uz(:, :, 1) - ...
    (2/(Rz + ZG))*p(:, :, 1);

    % update the pressure at all nodes
    p = p - pCx*(ux(:, 2:end, :) - ux(:, 1:end-1, :))...
        - pCy*(uy(2:end, :, :) - uy(1:end-1, :, :))...
        - pCz*(uz(:, :, 2:end) - uz(:, :, 1:end-1));
end