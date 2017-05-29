%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD2Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pd, udx, udy] = PSTD2Dfun(pd, udx, udy, diffmatrixX, diffmatrixY,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst)
    %% Function solves using the PSTD method for a pressure matrix...
    %  velocity matrix and differentiation impulse response in 2D
    %  and returns the solved pressure and velocity vectors
    
%% Velocity for y and x terms from 2D pressure matrix
    phat = fft2(pd);
    temp1 = phat .* diffmatrixX;
    temp2 = phat .* diffmatrixY;
    pdiffhatx = ifft2(temp1,'symmetric');
    pdiffhaty = ifft2(temp2,'symmetric');

%% Total Velocity
    udx = udx .* PMLdiff - PMLalphau .* (pdiffhatx.*PMLconst);
    udy = udy .* PMLdiff - PMLalphau .* (pdiffhaty.*PMLconst);
%% velocity differentiation in x direction
    uhat = fft2(udx);
    temp = uhat .* diffmatrixX;
    udiffhatx = ifft2(temp,'symmetric');
    
%% velocity differentiation in y direction
    uhat = fft2(udy);
    temp = uhat .* diffmatrixY;
    udiffhaty = ifft2(temp,'symmetric');

%% Total Pressure for new time step
    pd = pd .* PMLdiff -(PMLalphap .* (udiffhatx.*PMLconst))...
        - (PMLalphap .* (udiffhaty.*PMLconst));


end