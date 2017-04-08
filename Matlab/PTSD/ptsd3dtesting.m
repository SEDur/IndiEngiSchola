%% PTSD2D testing script

%% Initz
clc;
clear all;
% close all;

%% Make Variables


%alpha 
alphaXn = 1.0;
alphaXp = 1.0;
alphaYn = 1.0;
alphaYp = 1.0;
alphaZn = 1.0;
alphaZp = 1.0;

%define FS
fs = 1000.0;
%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 10.0;
%define grid width
gridWidth = 10.0;
%Target stability number 
St = 2/(pi*sqrt(3));
%define timestep
dt = (1/fs);
%dfine grid spacing
dx = c * dt * 1/St;
%calculate pconst
pconst = rho * c^2 * (dt/dx);
%calculate uconst
uconst = dt/(dx*rho);
%define pml depth 
% PMLdepth = ceil(abs(gridWidth/dx)/(2*(fs/c)));
% PMLdepth = 30;
PMLdepth = 60;

%calc time steps
timestep = abs(T/dt);
%calc grid size
N = ceil(abs(gridWidth/dx)+2*PMLdepth);
% temp = zeros(N, N);
%Calc source
% sStart = 44100 * 40;
src = zeros(1,ceil(T/dt)+1);
srctime = 0 : dt : 1.0;
srcf = 20;
src(10:10+length(srctime)-1) = (10^-12)*10^(50/20) * sin(2*pi*srcf.*srctime);
win = kaiser(length(srctime) + 20,2.0);
src(1:length(srctime) + 20) = src(1:length(srctime) + 20) .* win';
clear('win');
% music = audioread('track.mp3');
% src = (10^-12)*10^(50/20) .* music(sStart:sStart + length(src));
srcloc = ceil(N/2);
tempdiffmatrix = zeros(1,N);
spin = -180 :0.005 : 180;     
pd = zeros(N,N,N);
udx = zeros(N,N,N);
udy = zeros(N,N,N);
udz = zeros(N,N,N);
    for i2 = 1 : N
        if i2 <  ceil(N/2)
            tempdiffmatrix(i2) =  (i2-1);
        end
        if i2 ==  ceil(N/2)
            tempdiffmatrix(i2) = 0;
        end
        if i2 >  ceil(N/2)
            tempdiffmatrix(i2) = (i2 - 1 - N) ;
        end
    end

[mgx mgy mgz] = meshgrid(tempdiffmatrix);
diffmatrixX =  1i.* mgx ;
diffmatrixY =  1i.* mgy ;
diffmatrixZ =  1i.* mgz ;

PMLconst = ones(N,N,N);
PMLconst = PMLconst .* (pi*N);
PMLdiff = zeros(N,N,N);

    for i = 1 : PMLdepth
        PMLdiff(i,:,:) = i;
        PMLdiff(i,:,:) = (1.0/3.0).*(((PMLdepth-PMLdiff(i,:,:))./PMLdepth).^3);
    end
PMLdiff(N-PMLdepth+1:end,:,:) = fliplr(PMLdiff(PMLdepth:-1:1,:,:));
PMLdiff2 = permute(PMLdiff, [2 1 3]);
% PMLdiff3 = sqrt(PMLdiff.^2 + PMLdiff2.^2);
PMLdiff4 = permute(PMLdiff, [3 2 1]);
PMLdiff5 = ((PMLdiff) + (PMLdiff2) + (PMLdiff4))./3;
% PMLdiff5 = nthroot((PMLdiff.^3 + PMLdiff2.^3 + PMLdiff4.^3),3);


% PMLdiffmax = max(max(PMLdiff5));
% PMLdiffsetmax = 0.3011;
% PMLdiff5(PMLdiff5 > PMLdiffsetmax) = PMLdiffsetmax;
% mesh(PMLdiff);

% PMLdiff((N-PMLdepth+1):end) = 1 : PMLdepth;
% PMLdiff(1:PMLdepth) = (1/3.0).*(((PMLdepth-PMLdiff(1:PMLdepth))./PMLdepth).^3);
% PMLdiff((N-PMLdepth+1):end) = (1/3.0).*(PMLdiff((N-PMLdepth+1):end)./PMLdepth).^3;
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

%% solve for some time
% linkdata on;
% tic();
for i = 1 : T/dt
%     [pd, udx, udy] = PTSD3Dboundary(pd, udx, udy, udz, PMLdepth,...
%         xiXn, xiXp, xiYn, xiYp, xiZn, xiZp);
    [pd, udx, udy, udz] = PSTD3Dfun(pd, udx, udy, udz,...
        diffmatrixX, diffmatrixY , diffmatrixZ,...
     PMLdiff5, PMLalphau, PMLalphap, PMLconst);
    pd = PTSD3Dsrc(pd, src(i), srcloc);
    reciever(i) = pd(ceil(N/2), ceil(N/2), ceil(N/2));
    if mod(i, 100) < 1
    mesh(abs(pd(:, :, ceil(N/2))));
    
%     zlim([-10^-10 10^-10]);
%     set(gca,'zlim',[-10^-12 10^-12]);
    caxis([-10^-9 10^-9])
    shading interp;
    title(sprintf('Time = %.6f s',dt*(i-1)));
%     view([spin(i) 13]);
%      view(2);
    drawnow;
    end
end
% toc();

%% Display the results

