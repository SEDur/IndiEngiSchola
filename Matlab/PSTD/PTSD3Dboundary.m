%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD3Dboundary.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% This functtion takes arguments of pressure matrix, velocty matrices and 
% impedances for boundaries as a xi variable. This function is written to
% implement partially absorbing boundary conditions in 3D PSTD acoustic
% simulations.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pd, udx, udy, udz] = PTSD3Dboundary(pd, udx, udy, udz,...
    PMLdepth, xiXn, xiXp, xiYn, xiYp, xiZn, xiZp)
% This function is written to implement partially absorbing boundary 
% conditions in 3D PSTD acoustic simulations. pd, udx, and udy are the
% pressure and velocity matrices respectively, and the xi values are the
% normalised impedances are each boundary of a rectangular room.

if xiXn <= 1
    pd(PMLdepth+1,:,:) = pd(PMLdepth+1,:,:) .* xiXn;
else
    udy(PMLdepth+1,:,:) = udy(PMLdepth+1,:,:) .* (1/xiXn);
end

if xiXp <= 1
    pd(end-PMLdepth,:,:) = pd(end-PMLdepth,:,:) .* xiXp;
else
    udy(end-PMLdepth,:,:) = udy(end-PMLdepth,:,:) .* (1/xiXp);
end

if xiYn <= 1
    pd(:,PMLdepth+1,:) = pd(:,PMLdepth+1,:) .* xiYn;
else
    udx(:,PMLdepth+1,:) = udx(:,PMLdepth+1,:) .* (1/xiYn);
end

if xiYp <= 1
    pd(:,end-PMLdepth,:) = pd(:,end-PMLdepth,:) .* xiYp;
else
    udx(:,end-PMLdepth,:) = udx(:,end-PMLdepth,:) .* (1/xiYp);
end

if xiZn <= 1
    pd(:,:,PMLdepth+1) = pd(:,:,PMLdepth+1) .* xiZn;
else
    udz(:,:,PMLdepth+1) = udz(:,:,PMLdepth+1) .* (1/xiZn);
end

if xiZp <= 1
    pd(:,:,end-PMLdepth) = pd(:,:,end-PMLdepth) .* xiZp;
else
    udz(:,:,end-PMLdepth) = udz(:,:,end-PMLdepth) .* (1/xiZp);
end

end
