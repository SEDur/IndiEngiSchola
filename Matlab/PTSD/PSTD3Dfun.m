%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD3Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pd, udx, udy, udz] = PSTD3Dfun(pd, udx, udy, udz,...
    diffmatrixX, diffmatrixY, diffmatrixZ,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst)
    %% Function solves using the PSTD method for a pressure vector,...
    %  velocity vector and differentiation impulse response in 1 dimension
    %  and returns the solved pressure and velocity vectors
    
%% Velocity in 3d
    phat = fftn(pd);
    temp1 = phat .* diffmatrixX;
    temp2 = phat .* diffmatrixY;
    temp3 = phat .* diffmatrixZ;
    
    pdiffhatx = ifftn(temp1);
    pdiffhaty = ifftn(temp2);
    pdiffhatz = ifftn(temp3);
    
%% Total Velocity
    udx = udx .* PMLdiff - PMLalphau .* (pdiffhatx./PMLconst);
    udy = udy .* PMLdiff - PMLalphau .* (pdiffhaty./PMLconst);
    udz = udz .* PMLdiff - PMLalphau .* (pdiffhatz./PMLconst);
    
%% Pressure in 3d
    uhat = fftn(udx);
    temp = uhat .* diffmatrixX;
    udiffhatx = ifftn(temp);
    
%% Pressure in 3d
    uhat = fftn(udy);
    temp = uhat .* diffmatrixY;
    udiffhaty = ifftn(temp);
    
%% Pressure in 3d
    uhat = fftn(udz);
    temp = uhat .* diffmatrixY;
    udiffhatz = ifftn(temp);

%% Total Pressure
    pd = pd .* PMLdiff -(PMLalphap .* (udiffhatx./(PMLconst./2)))...
        - (PMLalphap .* (udiffhaty./(PMLconst./2)))...
        - (PMLalphap .* (udiffhatz./(PMLconst./2)));

end