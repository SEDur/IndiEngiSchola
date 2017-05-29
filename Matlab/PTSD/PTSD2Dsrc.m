%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD2Dsrc.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% This functtion takes arguments of pressure matrix, source term and source
% location, and returns the pressure matrix with a source term added.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pd = PTSD2Dsrc(pd, src, srcloc)
% This function applies the source term to a point in a pressure matrix,
% based on the values in the 'srcloc' source location array, and applied
% the inverse of the source term src to the pressure matrix.
pd(srcloc) = pd(srcloc) - src;

end