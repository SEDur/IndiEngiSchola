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
    %% Function solves using the PSTD method, takes pressure and velocity
    % matrices in 3D, and PML matrices and returns new velocity and
    % pressure matrices. 
    
%% Pressure differentials in X, Y and Z directions due to linear wave equations
    phat = fftn(pd);
    temp1 = phat .* diffmatrixX;
    temp2 = phat .* diffmatrixY;
    temp3 = phat .* diffmatrixZ;
    
    pdiffhatx = ifftn(temp1,'symmetric');
    pdiffhaty = ifftn(temp2,'symmetric');
    pdiffhatz = ifftn(temp3,'symmetric');
    
%% Total New Velocity
    udx = udx .* PMLdiff - PMLalphau .* (pdiffhatx./PMLconst);
    udy = udy .* PMLdiff - PMLalphau .* (pdiffhaty./PMLconst);
    udz = udz .* PMLdiff - PMLalphau .* (pdiffhatz./PMLconst);
    
%% Velocity differential in x dimension
    uhat = fftn(udx);
    temp = uhat .* diffmatrixX;
    udiffhatx = ifftn(temp,'symmetric');
    
%% Velocity differential in y dimension
    uhat = fftn(udy);
    temp = uhat .* diffmatrixY;
    udiffhaty = ifftn(temp,'symmetric');
    
%% Velocity differential in z dimension
    uhat = fftn(udz);
    temp = uhat .* diffmatrixZ;
    udiffhatz = ifftn(temp,'symmetric');

%% Total Pressure
    pd = pd .* PMLdiff -(PMLalphap .* (udiffhatx./(PMLconst)))...
        - (PMLalphap .* (udiffhaty./(PMLconst)))...
        - (PMLalphap .* (udiffhatz./(PMLconst)));

end