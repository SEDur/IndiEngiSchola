function PSTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0, roundtime, PMLdepth)

%         figure(1);
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