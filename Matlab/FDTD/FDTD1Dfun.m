%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDTD1Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pressure, velocityX] = FDTD1Dfun(pressure, pCx, velocityX, uCx, Rx, ZL, ZR)
    % update the non-boundary condition nodes for velocity
    velocityX(2:end-1) = velocityX(2:end-1) - uCx*(pressure(2:end) - pressure(1:end-1));

    % update the velocity at the right wall
    velocityX(end) = ((Rx - ZR)/(Rx + ZR))*velocityX(end) + (2/(Rx + ZR))*pressure(end);

    %update the velocity at the left wall
    velocityX(1) = ((Rx - ZL)/(Rx + ZL))*velocityX(1) - (2/(Rx + ZL))*pressure(1);

    % update the pressure at all nodes
    pressure = pressure - pCx*(velocityX(2:end) - velocityX(1:end-1));
end