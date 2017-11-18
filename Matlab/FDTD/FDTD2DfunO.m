%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD3DfunO.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A function that will execute one time step of the FDTD method in 2D,
% including forcing the velocity at some point to 0, thus acting like an
% obstacle to sound waves.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p, ux, uy] = FDTD2DfunO(p, pCx, pCy, ux, uy, uCx, uCy, Rx, Ry, ZL, ZR, ZT, ZB)

    % update the non-boundary condition nodes for velocity
    ux(:, 2:end-1) = ux(:, 2:end-1) - uCx*(p(:, 2:end) - p(:, 1:end-1));
    uy(2:end-1, :) = uy(2:end-1, :) - uCy*(p(2:end, :) - p(1:end-1, :));
    
    dimsx = size(ux)-1;
    dimsy = size(uy)-1;
    
    ux(ceil(dimsx(1)/4):ceil(dimsx(1)/2),...
        ceil(dimsx(1)/8)-1:ceil(dimsx(1)/4)-1) =...
        ux(ceil(dimsx(1)/4):ceil(dimsx(1)/2),...
        ceil(dimsx(1)/8)-1:ceil(dimsx(1)/4)-1) .* 0.0002;
    uy(ceil(dimsy(1)/4)-1:ceil(dimsy(1)/2)-1,...
        ceil(dimsx(1)/8):ceil(dimsx(1)/4)) =...
        uy(ceil(dimsy(1)/4)-1:ceil(dimsy(1)/2)-1,...
        ceil(dimsx(1)/8):ceil(dimsx(1)/4)) .* 0.0002;
    
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
    p = p - pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
        - pCy*(uy(2:end, :) - uy(1:end-1, :));
end