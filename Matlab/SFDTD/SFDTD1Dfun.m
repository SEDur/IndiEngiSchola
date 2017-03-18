%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SFDTD1Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pressure, ux] = SFDTD1Dfun(pressure, pCx, ux, uCx, Rx, ZL, ZR, idx, minWindowSize)

if(length(idx(idx > 0)) > minWindowSize)
    %  mat(col, row)
    for i = 2 : length(ux)-1
        if(idx(i) > 0)
            ux(i) = ux(i) - uCx*(pressure(i)-pressure(i-1));
        end
    end
else
    % update the non-boundary condition nodes for velocity
    ux(2:end-1) = ux(2:end-1) - uCx*(pressure(2:end) - pressure(1:end-1));
end

% update the velocity at the right wall
ux(end) = ((Rx - ZR)/(Rx + ZR))*ux(end) + (2/(Rx + ZR))*pressure(end);

%update the velocity at the left wall
ux(1) = ((Rx - ZL)/(Rx + ZL))*ux(1) - (2/(Rx + ZL))*pressure(1);

% update the pressure at all nodes
%     p(idx > 0) = p(idx > 0) - pCx*(ux(idx3 > 0) - ux(idx > 0))...
%         - pCy*(uy(idx2 > 0) - uy(idx > 0));
if(length(idx(idx > 0)) > minWindowSize)
    for i = 1 : length(pressure)
            if(idx(i) > 0)
                pressure(i) = pressure(i) - pCx*(ux(i+1) - ux(i));
            end
    end
else
    pressure = pressure - pCx*(ux(2:end) - ux(1:end-1));
end
end