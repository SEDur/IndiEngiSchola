%% PTSD3D testing script

%% Initz
clc;
clear all;
% close all;

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
% fs = 44100.0;
% fs = 10000;
fs = 1000;

%define density
rho = 1.21;
%define speed of sound
% c = 343;
c = 3430; %<< works for 2k
% c = 34300;
%define total time
T = 1.0;

%Target stability number 
St = 2/(pi*sqrt(3));
dx = (c / fs);
%define timestep
% dt = (1/fs);
dt = ((1/c)*dx)/2;
%dfine grid spacing
% dx = c * dt / St;
%calculate pconst
% pconst = pi * 100 * (rho * c^2 * dt/dx);
% pconst = rho * c^2 * dt;
% pconst = rho * c^2 * (dt/dx);
% pconst = dt/(dx*rho);
pconst = rho * c^2 * (dt/dx);
%calculate uconst
uconst = (1/rho)*(dt/dx);
% pconst = c^2*rho*dt/dx;
%calculate uconst
% uconst = dt/(dx*rho);
% uconst = (1/rho)*(dt/dx);
% uconst = dt/(dx*rho);
% uconst = c^2*rho*dt/dx;
%define pml depth 
% PMLdepth = ceil(abs(gridWidth/dx)/(2*(fs/c)));
% PMLdepth = 30;
PMLdepth = 30;

%calc time steps
timestep = abs(T/dt);

%Calc source
% w1 = window(@gausswin,0.4/(dt),2.5); 
precursor = -(sin(2*pi*0.5*(0:0.001:2))).*0.01;
chirp = dsp.Chirp(...
    'SweepDirection', 'Unidirectional', ...
    'TargetFrequency', ceil((fs/4)/2), ...
    'InitialFrequency', 10,...
    'TargetTime', 0.4, ...
    'SweepTime', 0.4, ...
    'SamplesPerFrame', 0.4/dt, ...
    'SampleRate', 1/dt);
% % plot(chirp());
% % source1 = fir(chirp());
% source1 = chirp();
% source1(1:76) = 0;
source1 = GenerateMLSSequence(2,9,0).*((2*10^-5)*10^(100/20));
T = length(source1)*dt;
% source1 = source1(1:(T/dt));
w1 = window(@gausswin,length(source1),2.5); 
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
% temp = zeros(N, N);

% source1 = source1 .* w1;
% source1 = [precursor'; source1];
% source1 = [source1; zeros((T/dt) - length(source1),1)].*((2*10^-5)*10^(100/20));

srcloc = [PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx)];
tempdiffmatrixX = zeros(1,Nx);
tempdiffmatrixY = zeros(1,Ny);
tempdiffmatrixZ = zeros(1,Nz);
% spin = -180 :0.005 : 180;     
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

% [mgx mgy mgz] = meshgrid(tempdiffmatrixY,tempdiffmatrixX,tempdiffmatrixZ);
mgx = reshape(tempdiffmatrixX,[Nx,1,1]);
mgy = reshape(tempdiffmatrixY,[1,Ny,1]);
mgz = reshape(tempdiffmatrixZ,[1,1,Nz]);
diffmatrixX =  1i.*mgx ;
diffmatrixY =  1i.*mgy ;
diffmatrixZ =  1i.*mgz ;

PMLconst = ones(Nx,Ny,Nz);
% PMLconst = PMLconst .* (pi*max([Ny Nx Nz]));
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
% PMLdiff2 = permute(PMLdiff, [2 1 3]);
% PMLdiff3 = sqrt(PMLdiff.^2 + PMLdiff2.^2);
% PMLdiff4 = permute(PMLdiff, [3 2 1]);
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
    reciever(i) = real(pd(ceil(Ny/2), ceil(Nx/2), ceil(Nz/2)));
    srcnorm(i) = pd(srcloc(1,1),srcloc(1,2),srcloc(1,3));
    roundtime(i) = toc();
    T - (i*dt)
    PSTD3Dplotdomain(pd, xcells, ycells, zcells, i, dt, p0, roundtime(i), PMLdepth);

end
%% Some really minor postprocessing
% Hd = postprocessingDCfilter;
norec = reciever ./ max(abs(reciever));
recanal = AnalyseMLSSequence(reciever',0,2,9,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(200),[],200,fs);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(200),[],200,fs);
filename = strcat('xwidth',num2str(gridWidthX),'.mat');
save(filename,'roundtime', 'norec', 'recanal', 'lpsd', 'lf', 'srcnrm', 'spsd', 'sf');
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


