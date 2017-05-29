%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD3Dsrc.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% This function is written to wrap adding the source term of a 3D pstd
% simulation into a simple function.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pd = PTSD3Dsrc(pd, src, srcloc)

pd(srcloc(1),srcloc(2),srcloc(3)) = pd(srcloc(1),srcloc(2),srcloc(3))+src;

end