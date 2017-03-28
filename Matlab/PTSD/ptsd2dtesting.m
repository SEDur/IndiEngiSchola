%% PTSD1D testing script

%% Initz
clc;
clear all;
close all;

%% Make Variables

%alpha 
a = 1.0;
%define FS
fs = 10000.0;
%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 4.0;
%define grid width
gridWidth = 10.0;
%define timestep
dt = 1/fs;
%dfine grid spacing
dx = c * sqrt(2) * dt;
%calculate pconst
pconst = rho * c^2 * (dt/dx);
%calculate uconst
uconst = (1/rho)*(dt/dx);
%define pml depth 
% PMLdepth = ceil(abs(gridWidth/dx)/(10*(fs/c)));
PMLdepth = 100;

%calc time steps
timestep = abs(T/dt);
%calc grid size
N = ceil(abs(gridWidth/dx)+2*PMLdepth);
%calculate differentiation matrix
tempdiffmatrix = zeros(1,N);
% temp = zeros(N, N);
%Calc source
src = zeros(1,ceil(T/dt)+1);
% src(10:2010) = (10^-12)*30*10^(50/20) * sin(2*(pi/2000)*(1:2001));
srcloc = ceil(N/2);
% 
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


src(10:4010) = (10^-12*10^(50/20)) * sin(2*(pi/8010)*(1:4001));
        
% alpha = 0;
% calculate geometry matricies
% phat = zeros(N,N);
% uhat = zeros(N,N);
% pdiffhat = zeros(N,N);
% udiffhat = zeros(N,N);
pd = zeros(N,N);
ud = zeros(N,N);
    for i2 = 1 : N-1
        if i2 <  ceil((N-2)/2)
            tempdiffmatrix(i2) =  (i2-1);
        end
        if i2 ==  ceil((N-1)/2)
            tempdiffmatrix(i2) = 0 * (1+0j);
        end
        if i2 >  ceil((N-1)/2)
            tempdiffmatrix(i2) = (i2 - 1 - N) ;
        end
    end
diffmatrix = 1i * tempdiffmatrix;
PMLconst = ones(N,N);
PMLconst = PMLconst .* (3.142*N);
PMLdiff = zeros(N,N);
for i = 1 : N
PMLdiff(i,1:PMLdepth) = 1:PMLdepth;
PMLdiff(i,1:PMLdepth) = (a/3.0).*(((PMLdepth-PMLdiff(i,1:PMLdepth))./PMLdepth).^3);
end
PMLdiff(:,N-PMLdepth+1:end) = fliplr(PMLdiff(:,1:PMLdepth));
% mesh(PMLdiff);
% for i = PMLdepth : N-PMLdepth+1
%     PMLdiff(1:PMLdepth,i) = 1:PMLdepth';
% PMLdiff((N-PMLdepth+1):end, i) = 1 : PMLdepth';
% end
PMLdiff2 = PMLdiff';
% mesh(PMLdiff2);

PMLdiff = sqrt(PMLdiff.^2 + PMLdiff2.^2);
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

linex = 0 : dx : (size(pd, 1)-1)*dx;
%% solve for some time
% linkdata on;
% tic();
for i = 1 : T/dt
   [pd, ud] = PSTD2Dfun(pd, ud, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst, N, PMLdepth);
    pd = PTSD2Dsrc(pd, src(i), srcloc);
    mesh(linex, linex, pd);
%     zlim([-4e-10 4e-10]);
% %     set(gca,'zlim',[-10^-12 10^-12]);
% %     caxis([-10^-12 10^-12])
%     shading interp;
%     title(sprintf('Time = %.6f s',dt*i));
    drawnow;
i
end
% toc();

%% Display the results

