%% PTSD2D testing script

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
fs = 5000.0;

%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 0.1;
%define grid width
gridWidthX = 5.0;
gridWidthY = 4.0;
gridWidthZ = 3.0;
%Target stability number 
St = 2/(pi*sqrt(3));
%define timestep
dt = (1/fs);
%dfine grid spacing
dx = c * dt / St;
%calculate pconst
pconst = pi * 100 * (rho * c^2 * dt/dx);
%calculate uconst
uconst = dt/(dx*rho);
%define pml depth 
% PMLdepth = ceil(abs(gridWidth/dx)/(2*(fs/c)));
% PMLdepth = 30;
PMLdepth = 30;

%calc time steps
timestep = abs(T/dt);
%calc grid size
Nx = ceil(abs(gridWidthX/dx)+2*PMLdepth);
Ny = ceil(abs(gridWidthY/dx)+2*PMLdepth);
Nz = ceil(abs(gridWidthZ/dx)+2*PMLdepth);
% temp = zeros(N, N);
%Calc source
% sStart = 44100 * 40;
% src = zeros(1,ceil(T/dt)+1);
% srctime = 0 : dt : 1.0;
% srcf = 20;
% src(10:10+length(srctime)-1) = (10^-12)*10^(50/20) * sin(2*pi*srcf.*srctime);
% win = kaiser(length(srctime) + 20,2.0);
% src(1:length(srctime) + 20) = src(1:length(srctime) + 20) .* win';
% clear('win');
% music = audioread('track.mp3');
% src = (10^-12)*10^(50/20) .* music(sStart:sStart + length(src));
tnum = ceil(T/dt);
fc = 0.25;     % Cutoff frequency (normalised 0.5=nyquist)
n0 = 10;        % Initial delay (samples)
sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
n=0:tnum;
source1=exp(-dt^2*(n-n0).^2/(2*sigma^2));
source1 = -(((2*10^-5)*10^(100/20))*source1);
source1(source1 > (2*10^-5)) = (2*10^-5);
% source1 = -source1;
% srcloc = PMLdepth + ceil(1/dx);
% srcloc = [ceil(Ny/2) ceil(Nx/2) ceil(Nz/2)];
srcloc = [PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx) PMLdepth+ceil(1/dx)];
tempdiffmatrixX = zeros(1,Nx);
tempdiffmatrixY = zeros(1,Ny);
tempdiffmatrixZ = zeros(1,Nz);
% spin = -180 :0.005 : 180;     
pd = ones(Nx,Ny,Nz).*p0;
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
xcells = 0 : dx : (Nx-1-PMLdepth)*dx;
ycells = 0 : dx : (Ny-1-PMLdepth)*dx;
zcells = 0 : dx : (Nz-1-PMLdepth)*dx;
% xcells = ceil(1/dx);
% ycells = ceil(1/dx);
% zcells = ceil(1/dx);
%% solve for some time
% linkdata on;
% tic();
for i = 1 : T/dt
    tic();
    [pd, udx, udy, udz] = PTSD3Dboundary(pd, udx, udy, udz, PMLdepth,...
        xiXn, xiXp, xiYn, xiYp, xiZn, xiZp);
    [pd, udx, udy, udz] = PSTD3Dfun(pd, udx, udy, udz,...
        diffmatrixX, diffmatrixY , diffmatrixZ,...
     PMLdiff5, PMLalphau, PMLalphap, PMLconst);
    pd = PTSD3Dsrc(pd, source1(i), srcloc);
    reciever(i) = pd(ceil(Ny/2), ceil(Nx/2), ceil(Nz/2));
    roundtime(i) = toc();
    PSTD3Dplotdomain(pd, xcells, ycells, zcells, i, dt, p0, roundtime(i), PMLdepth);

%     zlim([-10^-10 10^-10]);
%     set(gca,'zlim',[-10^-12 10^-12]);
%     caxis([-10^-9 10^-9])
%     shading interp;
%     title(sprintf('Time = %.6f s',dt*(i-1)));
% %     view([spin(i) 13]);
% %      view(2);
%     drawnow;

end
subplot(3,1,1);
plot(0:dt:((length(reciever)-1)*dt),reciever)
axis('tight')
subplot(3,1,2);
plot(0:dt:((length(reciever))*dt),source1)
% xlim([0 0.01])
subplot(3,1,3);
plot(0:dt:((length(reciever)-1)*dt),roundtime)
axis('tight')
% toc();

%% Display the results

