%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD3Dsources.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A function to wrap the source terms of an FDTD simulation
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = FDTD3Dsources(p,sourcelocations ,sourcesamples , sourcetypes)

p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3))...
    = p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3))...
    - sourcesamples(1,1);

end