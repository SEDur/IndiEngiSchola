%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD1Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pressure, velocity] = PSTD1DfunB(pressure, velocity, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst)
    %% Function solves using the PSTD method for a pressure vector,...
    %  velocity vector and differentiation impulse response in 1 dimension
    %  and returns the solved pressure and velocity vectors
    
    alphaL = 1.0;
    alphaR = 1.0;
    alphaF = 1.0;
    alphaB = 1.0;
    ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
    ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
    ZT = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
    ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));
 
    phat = fft(pressure);
    temp = phat .* diffmatrix;
    pdiffhat = ifft(temp);
    velocity = velocity .* PMLdiff  - PMLalphau  .* (pdiffhat ./PMLconst); 
    
        % update the velocity at the right wall
    ud(:, end) = ((Rx - ZR)/(Rx + ZR))*ud(:, end) ...
        + (2/(Rx + ZR))*p(:, end);

    %update the velocity at the left wall
    ud(:, 1) = ((Rx - ZL)/(Rx + ZL))*ud(:, 1) - (2/(Rx + ZL))*p(:, 1);

    %update the velocity at the top wall
    ud(end, :) = ((Ry - ZT)/(Ry + ZT))*ud(end, :) ...
        + (2/(Ry + ZT))*p(end, :);

    %update the velocity at the bottom wall
    ud(1, :) = ((Ry - ZB)/(Ry + ZB))*ud(1, :) - (2/(Ry + ZB))*p(1, :);

    uhat = fft(velocity);
    temp = uhat .* diffmatrix;
    udiffhat = ifft(temp);
    pressure  = pressure  .* PMLdiff  - PMLalphap  .* (udiffhat ./PMLconst);

    
end