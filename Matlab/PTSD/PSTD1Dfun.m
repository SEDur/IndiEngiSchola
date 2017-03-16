%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD1Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pd, ud] = PSTD1Dfun(pd, ud, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst)
    %% Function solves using the PSTD method for a pressure vector,...
    %  velocity vector and differentiation impulse response in 1 dimension
    %  and returns the solved pressure and velocity vectors
    
    phat = fft(pd);
    temp = phat .* diffmatrix;
    pdiffhat = ifft(temp);
    ud = ud .* PMLdiff - PMLalphau .* (pdiffhat./PMLconst); 
    uhat = fft(ud);
    temp = uhat .* diffmatrix;
    udiffhat = ifft(temp);
    pd = pd .* PMLdiff - PMLalphap .* (udiffhat./PMLconst);
    
end