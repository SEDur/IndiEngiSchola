f = 10000;
c = 343;
lambda = c/f;
dx = lambda/6;
spacex = 30 / dx;
spacey = 60 / dx;
spacez = 40 / dx;
volume = spacex * spacey * spacez;
dt = 1 / c * sqrt((1/(dx^2))+(1/(dx^2))+(1/(dx^2)));
numTimeSteps = 1/dt;
myArray = ones(ceil(spacex), ceil(spacey), ceil(spacez));
sizeOfArray = numel(myArray);
numBytes = numel(myArray) * 8;
numMegaBytesSingleMatrix = numBytes / 1e6;
minReqMemoryMegaBytes = numMegaBytesSingleMatrix * 4;
%% 

f = 500;
c = 343;
lambda = c/f;
dx = lambda/6;
spacex = 30 / dx;
spacey = 60 / dx;
spacez = 40 / dx;
volume = spacex * spacey * spacez;
dt = 1 / c * sqrt((1/(dx))^2+(1/(dx))^2+(1/(dx)^2));
numTimeSteps = 1/dt;
myArray = ones(ceil(spacex), ceil(spacey), ceil(spacez));
sizeOfArray = numel(myArray);
numBytes = numel(myArray) * 8;
numMegaBytesSingleMatrix = numBytes / 1e6;
minReqMemoryMegaBytes = numMegaBytesSingleMatrix * 4;
%% 

f = 1000 : 1000 : 20000;
c = 343;
lambda = c./f;
dx = lambda./6;
spacex = 30 * 1./dx;
spacey = 60 * 1./dx;
spacez = 40 * 1./dx;
volume = spacex .* spacey .* spacez;
dt = 1 / c .* (sqrt((1./(dx.^2))+(1./(dx.^2))+(1./(dx.^2))));
numTimeSteps = 1./dt;
myArraySize = ceil(spacex).* ceil(spacey) .* ceil(spacez);

scatter(f / 1000, myArraySize, '*')
set(gcf, 'Color', 'White');
grid;
xlabel('Frequency [kHz]')
ylabel('Array Size [Elements]');
ylabel('Matrix Size [Elements]');
title('Matrix Size Vs Maximum Valid Simulation Data Frequency');
% numberOfElements = numel(myArray);
% numBytes = numel(myArray) * 4;

%% 
clear all;

c = 200 : 1 : 400;
f = 500;
lambda = c./f;
dx = lambda./6;
dy = lambda./6;
dz = lambda./6;

dt = (1 ./ c) .* sqrt((1./(dx)).^2+(1./(dy)).^2+(1./(dz).^2));

plot(c, dt)

The fundamental problem here is that we have to deal with large constructs of memory, which means we need a large amount of memory available and hardware needs to be able to work with large amount of memory quickly. This is not trivial, and as we will see shortly, a fundamental problem.

%% 
clear all;

f = 500;
c = 343;
lambda = c/f;

dx = lambda/6;
dy = lambda/6;
dz = lambda/6;

dimSize = 30 : 1 : 500;
spacex = dimSize ./ dx;
spacey = dimSize ./ dy;
spacez = dimSize ./ dz;
volume = spacex .* spacey .* spacez;

sizeOfArray = ceil(volume);
numBytes = sizeOfArray * 8;
numMegaBytesSingleMatrix = numBytes / 1e6;
minReqMemoryMegaBytes = numMegaBytesSingleMatrix * 4;

scatter(dimSize, sizeOfArray, '*')
set(gcf, 'Color', 'White');
grid;
xlabel('Domain Size (m^3)')
ylabel('Matrix Size [Elements]');
title('Domain Size Vs Matrix Volume for a Rectangular 500Hz FDTD simulation');
axis('tight')

%% 
clear all;

f = 500;
c = 343;
lambda = c/f;

dx = lambda/6;
dy = lambda/6;
dz = lambda/6;

dimSize = 30 : 1 : 500;
spacex = dimSize ./ dx;
spacey = dimSize ./ dy;
spacez = dimSize ./ dz;
volume = spacex .* spacey .* spacez;

sizeOfArray = ceil(volume);
numBytes = sizeOfArray * 8;
numMegaBytesSingleMatrix = numBytes / 1e6;
minReqMemoryMegaBytes = numMegaBytesSingleMatrix * 4;

scatter(dimSize, numMegaBytesSingleMatrix, '*')
set(gcf, 'Color', 'White');
grid;
xlabel('Domain Size (m^3)')
ylabel('Matrix Size [MegaBytes]');
title('Domain Size Vs Matrix Volume for a Rectangular 500Hz FDTD simulation');
axis('tight')