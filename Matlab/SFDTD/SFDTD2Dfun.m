function [p, ux, uy] = SFDTD2Dfun(p, pCx, pCy, ux, uy, uCx, uCy, Rx, Ry, ZL, ZR, ZT, ZB, idx)

    % update the non-boundary condition nodes for velocity
    ux(:, 2:end-1) = ux(:, 2:end-1) - uCx*(p(:, 2:end) - p(:, 1:end-1));
    uy(2:end-1, :) = uy(2:end-1, :) - uCy*(p(2:end, :) - p(1:end-1, :));

%     ux(pidxRow, pidxCol+1) = ux(pidxRow, pidxCol+1) - uCx*(p(pidxRow, pidxCol+1) - p(pidxRow, pidxCol));
%     uy(pidxRow+1, pidxCol) = uy(pidxRow+1, pidxCol) - uCy*(p(pidxRow+1, pidxCol) - p(pidxRow, pidxCol));


    % update the velocity at the right wall
    ux(:, end) = ((Rx - ZR)/(Rx + ZR))*ux(:, end) ...
        + (2/(Rx + ZR))*p(:, end);

    %update the velocity at the left wall
    ux(:, 1) = ((Rx - ZL)/(Rx + ZL))*ux(:, 1) - (2/(Rx + ZL))*p(:, 1);

    %update the velocity at the top wall
    uy(end, :) = ((Ry - ZT)/(Ry + ZT))*uy(end, :) ...
        + (2/(Ry + ZT))*p(end, :);

    %update the velocity at the bottom wall
    uy(1, :) = ((Ry - ZB)/(Ry + ZB))*uy(1, :) - (2/(Ry + ZB))*p(1, :);

    % update the pressure at all nodes
%     p(pidxRow, pidxCol) = p(pidxRow, pidxCol) - pCx*(ux(pidxRow, pidxCol+1) - ux(pidxRow, pidxCol))...
%         - pCy*(uy(pidxRow+1, pidxCol) - uy(pidxRow, pidxCol));
   p = p - pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
        - pCy*(uy(2:end, :) - uy(1:end-1, :));
end