%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD3Dplotdomain.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% This function is used to consolidate plotting of a 3D PSTD simulations
% pressure matrix, down into an easy to manage function with a restricted
% namespace.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PSTD3Dplotdomain(p, n, dt, p0, roundtime, PMLdepth)
% This function is written to consolidate the plotting of slices of a 3D
% pressure domain into a simple function.

        xslice = ceil((size(p,2)-(2*PMLdepth))/2); 
        yslice = ceil((size(p,1)-(2*PMLdepth))/2);  
        zslice = ceil((size(p,3)-(2*PMLdepth))/3);
        
        slice(p(PMLdepth:end-PMLdepth,...
            PMLdepth:end-PMLdepth,...
            PMLdepth:end-PMLdepth),...
            xslice, yslice, zslice)
        
%           mesh(p(:,:,ceil(size(p,3)/2)));
%         slice(xvec,yvec,zvec,abs(p),xslice,yslice,zslice)
%           view(2);
%         zslice = (ceil(zcells-1)/2);
%         figure(1);
%         slice(abs(p),[xcells-1 ceil(xcells/2)],ycells-1,zslice) 
%         shading interp;
        title(sprintf('Time = %.6f s, Max P = %.3f dB, ExecTime = %.4f',n*dt,20*log10(real(max(max(max(abs(p(:,:,:))))))/p0),roundtime),...
            'Color',[0 0 0],'FontSize', 14);
%         xlabel('Width (meters)', 'Color', [0 0 0]);
%         ylabel('Length (meters)', 'Color', [0 0 0]);
        shading('interp');
        drawnow;

end