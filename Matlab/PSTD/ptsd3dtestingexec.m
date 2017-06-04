%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PTSD3Dtestingexec.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that was written to test the execution speed of the PSTD method
% for 3 different domain sizes.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initz
clc;
clear all;
% close all;
addpath(genpath('../mls/mls/'));
%% Make Variables
p0 = 2*10^-5;

%alpha 
alphaXn = 0.45;
alphaXp = 0.45;
alphaYn = 0.45;
alphaYp = 0.45;
alphaZn = 0.45;
alphaZp = 0.45;

%define FS
fs = 3000;

%define density
rho = 1.21;
%define speed of sound
c = 3430; %<< works for 2k
%define total time
T = 1.0;

%Target stability number 
St = 2/(pi*sqrt(3));
dx = (c / fs);

%define timestep
dt = ((1/c)*dx)/2;

%calculate pconst
pconst = rho * c^2 * (dt/dx);

%calculate uconst
uconst = (1/rho)*(dt/dx);

%set PML depth
PMLdepth = 30;

%calc time steps
timestep = abs(T/dt);

%Calc source
source1 = GenerateMLSSequence(2,11,0).*((2*10^-5)*10^(100/20));
T = length(source1)*dt;

%set widths
widths = [5 10 20 40 60];

for cntr = 1 : 5
%define grid width
gridWidthX = widths(cntr);
gridWidthY = gridWidthX;
gridWidthZ = gridWidthX;

%calc grid size
Nx = ceil(abs(gridWidthX/dx)+2*PMLdepth);
Ny = ceil(abs(gridWidthY/dx)+2*PMLdepth);
Nz = ceil(abs(gridWidthZ/dx)+2*PMLdepth);

%set source location
srcloc = [PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx)];

%preallocate differentiation matrices
tempdiffmatrixX = zeros(1,Nx);
tempdiffmatrixY = zeros(1,Ny);
tempdiffmatrixZ = zeros(1,Nz);

%preallocate domain matrices
pd =  zeros(Nx,Ny,Nz);
udx = zeros(Nx,Ny,Nz);
udy = zeros(Nx,Ny,Nz);
udz = zeros(Nx,Ny,Nz);

    for i2 = 1 : Nx-1
        if i2 <  ceil(Nx/2)
            tempdiffmatrixX(i2) =  (i2-1);
        end
        if i2 ==  ceil(Nx/2)
            tempdiffmatrixX(i2) = 0;
        end
        if i2 >  ceil(Nx/2)
            tempdiffmatrixX(i2) = (i2 - 1 - Nx) ;
        end
    end
    
    for i2 = 1 : Ny-1
        if i2 <  ceil(Ny/2)
            tempdiffmatrixY(i2) =  (i2-1);
        end
        if i2 ==  ceil(Ny/2)
            tempdiffmatrixY(i2) = 0;
        end
        if i2 >  ceil(Ny/2)
            tempdiffmatrixY(i2) = (i2 - 1 - Ny) ;
        end
    end
    for i2 = 1 : Nz-1
        if i2 <  ceil(Nz/2)
            tempdiffmatrixZ(i2) =  (i2-1);
        end
        if i2 ==  ceil(Nz/2)
            tempdiffmatrixZ(i2) = 0;
        end
        if i2 >  ceil(Nz/2)
            tempdiffmatrixZ(i2) = (i2 - 1 - Nz) ;
        end
    end

mgx = reshape(tempdiffmatrixX,[Nx,1,1]);
mgy = reshape(tempdiffmatrixY,[1,Ny,1]);
mgz = reshape(tempdiffmatrixZ,[1,1,Nz]);
diffmatrixX =  1i.*mgx ;
diffmatrixY =  1i.*mgy ;
diffmatrixZ =  1i.*mgz ;

PMLconst = ones(Nx,Ny,Nz);
PMLconst = PMLconst * (pi*sqrt(Nx^2 + Ny^2 + Nz^2));
PMLdiff = zeros(Nx,Ny,Nz);
PMLdiff2 = zeros(Nx,Ny,Nz);
PMLdiff4 = zeros(Nx,Ny,Nz);
    for i = 1 : PMLdepth
        PMLdiff(i,:,:) = i;
        PMLdiff(i,:,:) = (1.0/3.0).*(((PMLdepth-PMLdiff(i,:,:))./PMLdepth).^3);
    end
PMLdiff(Nx-PMLdepth+1:end,:,:) = fliplr(PMLdiff(PMLdepth:-1:1,:,:));

    for i = 1 : PMLdepth
        PMLdiff2(:,i,:) = i;
        PMLdiff2(:,i,:) = (1.0/3.0).*(((PMLdepth-PMLdiff2(:,i,:))./PMLdepth).^3);
    end
PMLdiff2(:,end:-1:Ny-PMLdepth+1,:) = fliplr(PMLdiff2(:,PMLdepth:-1:1,:));
    for i = 1 : PMLdepth
        PMLdiff4(:,:,i) = i;
        PMLdiff4(:,:,i) = (1.0/3.0).*(((PMLdepth-PMLdiff4(:,:,i))./PMLdepth).^3);
    end
PMLdiff4(:,:,Nz-PMLdepth+1:end) = fliplr(PMLdiff4(:,:,PMLdepth:-1:1));

PMLdiff5 = ((PMLdiff) + (PMLdiff2) + (PMLdiff4))./3;
clear('PMLdiff');
clear('PMLdiff2');
clear('PMLdiff3');
clear('PMLdiff4');

PMLalphau = uconst*(1./(1+PMLdiff5));
PMLalphap = pconst*(1./(1+PMLdiff5));
PMLdiff5 = ((1-PMLdiff5)./(1+PMLdiff5));

%set up values for partially absorbing boundaries
S = c*(dt/dx);
Rxn = 1 - alphaXn;
Rxp = 1 - alphaXp;
Ryn = 1 - alphaYn;
Ryp = 1 - alphaYp;
Rzn = 1 - alphaZn;
Rzp = 1 - alphaZp;
xiXn = (1 + Rxn)/(1 + Rxn - 2 * S * Rxn);
xiXp = (1 + Rxp)/(1 + Rxn - 2 * S * Rxp);
xiYn = (1 + Ryn)/(1 + Ryn - 2 * S * Ryn);
xiYp = (1 + Ryp)/(1 + Ryp - 2 * S * Ryp);
xiZn = (1 + Rzn)/(1 + Rzn - 2 * S * Rzn);
xiZp = (1 + Rzp)/(1 + Rzp - 2 * S * Rzp);
xcells = 0 : dx : (Nx-1-PMLdepth)*dx;
ycells = 0 : dx : (Ny-1-PMLdepth)*dx;
zcells = 0 : dx : (Nz-1-PMLdepth)*dx;

%% solve for some time

for i = 1 : T/dt
    tic();
    [pd, udx, udy, udz] = PTSD3Dboundary(pd, udx, udy, udz, PMLdepth,...
        xiXn, xiXp, xiYn, xiYp, xiZn, xiZp);
    [pd, udx, udy, udz] = PSTD3Dfun(pd, udx, udy, udz,...
        diffmatrixX, diffmatrixY , diffmatrixZ,...
     PMLdiff5, PMLalphau, PMLalphap, PMLconst);
    pd = PTSD3Dsrc(pd, source1(i), srcloc);
    reciever(i,1) = pd(ceil(Nx/4), ceil(Ny/4),ceil(Nz/2));
    reciever(i,2) = pd(ceil(Nx/2)+ceil(Nx/4), ceil(Ny/4),ceil(Nz/2));    
    reciever(i,3) = pd(ceil(Nx/2), ceil(Ny/2)+ceil(Nx/4),ceil(Nz/2));
    reciever(i,4) = pd(ceil(Nx/2+ceil(Nx/4)), ceil(Ny/2)+ceil(Nx/4),ceil(Nz/2));
    reciever(i,5) = pd(ceil(Nx/2), ceil(Ny/2),ceil(Nz/2));
    srcnorm(i) = pd(srcloc(1,1),srcloc(1,2),srcloc(1,3));
    exectime(i) = toc();
    T - (i*dt)
    PSTD3Dplotdomain(pd, i, dt, p0, exectime(i), PMLdepth)

end
%% Some really minor postprocessing
% Hd = postprocessingDCfilter;
norec = reciever ./ max(abs(reciever));
recanal = AnalyseMLSSequence(reciever(:,1),0,2,5,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(ceil(length(norec)/2)),[],ceil(length(norec)/2),fs);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(ceil(length(srcnrm)/2)),[],ceil(length(srcnrm)/2),fs);
filename = strcat('xwidth',num2str(gridWidthX),'.mat');
save(filename,'exectime', 'norec', 'recanal', 'lpsd', 'lf', 'srcnrm', 'spsd', 'sf');
end
%% Display the results
subplot(5,1,1);
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold on;
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold off;
axis('tight')
legend('reciever','source');
title('raw input and output');
subplot(5,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm)
hold off;
axis('tight')
legend('reciever','source');
title('normalised input and output');
subplot(5,1,3);
plot(lf, db(lpsd));
hold on;
plot(sf, db(spsd));
hold off;
legend('reciever','source');
grid('on');
title('power spectral density of input and output');
subplot(5,1,4);
plot(0:dt:((length(reciever)-1)*dt),roundtime)
axis('tight')
ttlstr = sprintf('computation time per cycle, total time is %d',sum(roundtime));
title(ttlstr);
subplot(5,1,5);
plot(0:dt:((length(recanal)-1)*dt),recanal);
title('MLS Analysed');


