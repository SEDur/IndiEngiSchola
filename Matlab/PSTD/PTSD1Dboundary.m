%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD1Dboundary.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% This functtion takes arguments of pressure array, velocty array and 
% impedances for boundaries as a xi variable. This function is written to
% implement partially absorbing boundary conditions in 1D PSTD acoustic
% simulations.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pd, ud] = PTSD1Dboundary(pd, ud, PMLdepth, xiL, xiR)
% This function is written to implement partially absorbing boundary 
% conditions in 1D PSTD acoustic simulations. pd and ud are the
% pressure and velocity arrays respectively, and the xi values are the
% normalised impedances at each boundary of the domain.

if xiL <= 1
    pd(PMLdepth+1) = pd(PMLdepth+1) * xiL;
else
    ud(PMLdepth+1) = ud(PMLdepth+1) * (1/xiL);
end

if xiR <= 1
    pd(end-PMLdepth) = pd(end-PMLdepth) * xiR;
else
    ud(end-PMLdepth) = ud(end-PMLdepth) * (1/xiR);
end

end
