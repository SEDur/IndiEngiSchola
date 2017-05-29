%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PTSD2Dtesting.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that will execute the PSTD method for a 2D simulation up to 10kHz
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initz
clc;
clear all;
% close all;

%% Make Variables


%alpha 
alphaXn = 0.45;
alphaXp = 0.45;
alphaYn = 0.45;
alphaYp = 0.45;

%define FS
fs = 10000;
% fmax = 5000;
%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 0.1;
%define grid width
gridWidthX = 34.0;
gridWidthY = 34.0;
%Define Stability Condition
St = 2/(pi * sqrt(2));
%define timestep
% dx = (c / (pi*fmax));
% dt = ((1/c)*dx)/2;
% dt = 1 / fs;
dt = 1/(2*fs);
% dx = c * sqrt(2) * dt;
% dt = (1/fs);
% dt = 1/(2*fs);
% % %dfine grid spacing
% dx = c * dt / St;
% dx = c * sqrt(2) * dt;
dx = 2 * dt * c;
% assert(isequal((c*dt/dx),St));
%calculate pconst
pconst = rho * c^2 * (dt/dx);
%calculate uconst
uconst = (1/rho)*(dt/dx);
% pconst = c^2 * rho * dt/dx;
% pconst = rho * c^2 * (dt/dx) * dt * c;
%calculate uconst
% uconst = dt / (dx*rho);
% uconst = (1/rho) * (dt/dx) * dt * c;
%define pml depth 
% PMLdepth = ceil(abs(gridWidth/dx)/(2*(fs/c)));
PMLdepth = 30;
%calc time steps
timesteps = abs(T/dt);
%calc grid size
Nx = ceil(abs(gridWidthX/dx)+2*PMLdepth);
Ny = ceil(abs(gridWidthY/dx)+2*PMLdepth);
%calculate differentiation matrix
tempdiffmatrixX = zeros(1,Nx);
tempdiffmatrixY = zeros(1,Ny);

% temp = zeros(N, N);
%Calc source
% sStart = 44100 * 40;
% src = zeros(1,ceil(T/dt)+1);
% srctime = 0 : dt/2 : 0.2;
% srcf = 10;
% src(10:10+length(srctime)-1) = (10^-12)*10^(50/20) * sin(2*pi*srcf.*srctime);
% win = kaiser(length(srctime) + 20,2.0);
% src(1:length(srctime) + 20) = src(1:length(srctime) + 20) .* win';
% clear('win');

% tnum = ceil(T/dt);
% src = zeros(1,tnum);
% fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 10;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% src = exp(-dt^2*(n-n0).^2/(2*sigma^2));
% scr = ((2*10^-5)*10^(100/20))*src;

% music = audioread('track.mp3');
% src = (10^-12)*10^(50/20) .* music(sStart:sStart + length(src));
srcloc = [PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx)];
source1 = GenerateMLSSequence(2,11,0).*((2*10^-5)*10^(100/20));
w1 = window(@gausswin,length(source1),2.5); 
src = source1;
T = length(src)*dt;
% source1 = source1(1:(T/dt));

% tnum = ceil(T/dt);
% fc = 0.25;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 30;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% src=exp(-dt^2*(n-n0).^2/(2*sigma^2));
% for n = 37 : length(src)
%     if(src(n) < 0)
%        src(n) = 0; 
%     end
% end

% spin = -180 :0.005 : 180;
        
% alpha = 0;
% calculate geometry matricies
% phat = zeros(N,N);
% uhat = zeros(N,N);
% pdiffhat = zeros(N,N);
% udiffhat = zeros(N,N);
% pd = zeros(Nx,Ny);
udx = zeros(Nx,Ny);
udy = zeros(Nx,Ny);
pd = zeros(Nx,Ny);

% linex = 0 : dx : dx * (Nx-1-(2*PMLdepth));
% liney = 0 : dx : dx * (Ny-1-(2*PMLdepth));

linex = 0 : dx : dx * (Nx-1);
liney = 0 : dx : dx * (Ny-1);

% udy = zeros(N,N);
    for i2 = 1 : Nx-1
        if i2 <  ceil((Nx-2)/2)
            tempdiffmatrixX(i2) =  (i2-1);
        end
        if i2 ==  ceil((Nx-1)/2)
            tempdiffmatrixX(i2) = 0 * (1+0j);
        end
        if i2 >  ceil((Nx-1)/2)
            tempdiffmatrixX(i2) = (i2 - 1 - Nx) ;
        end
    end
    
     for i2 = 1 : Ny-1
        if i2 <  ceil((Ny-2)/2)
            tempdiffmatrixY(i2) =  (i2-1);
        end
        if i2 ==  ceil((Ny-1)/2)
            tempdiffmatrixY(i2) = 0 * (1+0j);
        end
        if i2 >  ceil((Ny-1)/2)
            tempdiffmatrixY(i2) = (i2 - 1 - Ny) ;
        end
    end
[mgx, mgy] = meshgrid(tempdiffmatrixY, tempdiffmatrixX);
diffmatrixX =  1i.*mgx ;
diffmatrixY =  1i.*mgy ;

PMLconst = ones(Nx,Ny);
PMLconst = 1./(PMLconst .* (pi*sqrt(Nx^2 + Ny^2)));
PMLdiff = zeros(Nx,Ny);
PMLdiff2 = zeros(Nx,Ny);
for i = 1 : Nx
PMLdiff(i,1:PMLdepth) = 1:PMLdepth;
PMLdiff(i,1:PMLdepth) = (1.0/3.0).*(((PMLdepth-PMLdiff(i,1:PMLdepth))./PMLdepth).^3);
end
PMLdiff(:,Ny-PMLdepth+1:end) = fliplr(PMLdiff(:,1:PMLdepth));
for i = 1 : Ny
PMLdiff2(1:PMLdepth,i) = 1:PMLdepth;
PMLdiff2(1:PMLdepth,i) = (1.0/3.0).*(((PMLdepth-PMLdiff2(1:PMLdepth,i))./PMLdepth).^3);
end
PMLdiff2(end:-1:Nx-PMLdepth+1,:) = fliplr(PMLdiff2(1:PMLdepth,:));
% PMLdiff2 = PMLdiff';
% mesh(PMLdiff2);

PMLdiff = sqrt(PMLdiff.^2 + PMLdiff2.^2);
%%
PMLdiffmax = max(max(PMLdiff));
PMLdiffsetmax = 0.3011;
PMLdiff(PMLdiff > PMLdiffsetmax) = PMLdiffsetmax;
% mesh(PMLdiff);

% PMLdiff((N-PMLdepth+1):end) = 1 : PMLdepth;
% PMLdiff(1:PMLdepth) = (1/3.0).*(((PMLdepth-PMLdiff(1:PMLdepth))./PMLdepth).^3);
% PMLdiff((N-PMLdepth+1):end) = (1/3.0).*(PMLdiff((N-PMLdepth+1):end)./PMLdepth).^3;
PMLalphau = uconst*(1./(1+PMLdiff));
PMLalphap = pconst*(1./(1+PMLdiff));
PMLdiff = ((1-PMLdiff)./(1+PMLdiff));

S = c*(dt/dx);
Rxn = 1 - alphaXn;
Rxp = 1 - alphaXp;
Ryn = 1 - alphaYn;
Ryp = 1 - alphaYp;
xiXn = (1 + Rxn)/(1 + Rxn - 2 * S * Rxn);
xiXp = (1 + Rxp)/(1 + Rxn - 2 * S * Rxp);
xiYn = (1 + Ryn)/(1 + Ryn - 2 * S * Ryn);
xiYp = (1 + Ryp)/(1 + Ryp - 2 * S * Ryp);

%% solve for some time
% linkdata on;
% tic();
for i = 1 : T/dt+1
    tic();
    [pd, udx, udy] = PTSD2Dboundary(pd, udx, udy, PMLdepth,...
        xiXn, xiXp, xiYn, xiYp);
    [pd, udx, udy] = PSTD2Dfun(pd, udx, udy, diffmatrixX,diffmatrixY,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst);
    pd = PTSD2Dsrc(pd, src(i), srcloc);
    srcnorm(i,1) = pd(srcloc(1), srcloc(2));
    reciever(i,1) = pd(ceil(Nx/4), ceil(Ny/4));
    reciever(i,2) = pd(ceil(Nx/2)+ceil(Nx/4), ceil(Ny/4));    
    reciever(i,3) = pd(ceil(Nx/2), ceil(Ny/2)+ceil(Nx/4));
    reciever(i,4) = pd(ceil(Nx/2+ceil(Nx/4)), ceil(Ny/2)+ceil(Nx/4));
    reciever(i,5) = pd(ceil(Nx/2), ceil(Ny/2));
    exTime(i) = toc();
%     if mod(i, 100) < 1
%     mesh(liney, linex, real(pd(PMLdepth:end-PMLdepth-1,...
%         PMLdepth:end-PMLdepth-1)));
% mesh(liney, linex, abs(pd));
mesh(liney, linex, pd);

    
%     zlim([-0.4 0.4]);
     zlim([-1 1]);
%     set(gca,'zlim',[-0.04 0.04]);
    caxis([-0.04 0.04]);
    shading interp;
    title(sprintf('Time = %.6f s,ExecTime = %.4f',dt*(i-1),exTime(i)));
%     view([spin(i) 13]);
    view(2)
    axis tight;
    zlim([-1 1]);
    drawnow;
%     end
end
%% Some really minor postprocessing
% Hd = postprocessingDCfilter;
norec = reciever ./ max(abs(reciever));
% recanal = AnalyseMLSSequence(reciever(:,1)',0,2,11,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(2000),[],2000,fs);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(2000),[],2000,fs);
for i = 1 : size(norec,2)
lags(i) = getlag(norec(:,i),srcnrm);
norec(:,i) = circshift(norec(:,i),lags(i));
end

%% Display the results
subplot(4,1,1);
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold on;
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold off;
axis('tight')
legend('source','topleft','topright','bottomleft','bottomright','center');
title('Raw Input And Output');
subplot(4,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm)
hold off;
axis('tight')
legend('topleft','topright','bottomleft','bottomright','center','source');
title('Normalised Input And Output');
subplot(4,1,3);
plot(lf, db(lpsd),'--','Linewidth',2.0);
hold on;
plot(sf, db(spsd));
hold off;
legend('topleft','topright','bottomleft','bottomright','center','source');
grid('on');
title('Power Spectral Density of Input and Output');
subplot(4,1,4);
plot(0:dt:((length(reciever)-1)*dt),exTime)
axis('tight')
ttlstr = sprintf('Computation Time Per Cycle, Total Time: %i',sum(exTime));
title(ttlstr);
% 
%% Display the results

