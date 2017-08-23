%% Make Variables
%Create Figure
figure(1);
set(1,'Windowstyle','docked');
%Define pressure mimimum
p0 = 2*10^-5;
%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 0.12;
%define FS
fs = 44100;
sfdtdcstab = 2/(pi*sqrt(2));
cstab = 2/(pi*sqrt(2));
%define timestep
pstddt = (1/fs);
fdtddt = (1/fs)/3;
sfdtddt = (1/fs)/3;
%grid size
fdtdgx = c * fdtddt / cstab;
fdtdgy = c * fdtddt / cstab;
fdtdgz = c * fdtddt / cstab;
%grid size
sfdtdgx = c * sfdtddt / sfdtdcstab;
sfdtdgy = c * sfdtddt / sfdtdcstab;
sfdtdgz = c * sfdtddt / sfdtdcstab;
%define grid width
dims = [32 64 128 256 512 1024 2048 4096 8192];
gridWidthX = dims;
gridWidthY = dims;
gridWidthZ = dims;
lx = dims;
ly = dims;
lz = dims;
fdtdxcells = ceil(lx/fdtdgx);
fdtdycells = ceil(ly/fdtdgy);
fdtdzcells = ceil(lz/fdtdgz);
sfdtdxcells = ceil(lx/sfdtdgx);
sfdtdycells = ceil(lx/sfdtdgy);
sfdtdzcells = ceil(lx/sfdtdgz);
%Target stability number 
pstdSt = 2/(pi*sqrt(3));
%dfine grid spacing
pstddx = c * pstddt / pstdSt;
%define pml depth 
PMLdepth = 30;
%calc grid size
Nx = ceil(abs(gridWidthX/pstddx)+2*PMLdepth);
Ny = ceil(abs(gridWidthY/pstddx)+2*PMLdepth);
Nz = ceil(abs(gridWidthZ/pstddx)+2*PMLdepth);
pstdncells = Nx.*Ny.*Nz;
pstd2dncells = Nx.*Ny;
fdtdncells = fdtdxcells .* fdtdycells .* fdtdzcells;
fdtd2dncells = fdtdxcells .* fdtdycells;
sfdtdncells = sfdtdxcells .* sfdtdycells .* sfdtdzcells;
sfdtd2dncells = sfdtdxcells .* sfdtdycells;

semilogy(dims,[fdtd2dncells; pstd2dncells; sfdtd2dncells])
legend('fdtd','pstd','sfdtd');
axis('tight');
set(gca,'xscale','log');
set(gcf,'Color','White');
xlabel('N Dimension size');
ylabel('Number of Cells in Domain');
title('Number of cells Compared to Domain Size for FDTD, SFDTD and PSTD simulations');
grid('on');