%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD3Dplotdomain.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A function that wraps plotting of a 3D fDTD simulation domain.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FDTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0)

        xslice = xcells-1; 
        yslice = ycells-1; 
        zslice = (ceil(zcells-1)/2);
        slice(p,xslice,yslice,zslice)
        shading interp;
        title(sprintf('Time = %.6f s Max P = %.3f dB',n*dt,20*log10(real(max(max(max(abs(p(:,:,:))))))/p0)),...
            'Color',[0 0 0],'FontSize', 14);
        shading('interp');
        drawnow;

end