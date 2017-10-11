f = 500;
c = 343;
lambda = c/f;
dx = lambda/6;
spacex = 30 * 1/dx;
spacey = 60 * 1/dx;
spacez = 40 * 1/dx;
volume = spacex * spacey * spacez;
dt = 1 / c * sqrt((1/(dx^2))+(1/(dx^2))+(1/(dx^2)));
numTimeSteps = 1/dt;
myArray = ones(ceil(spacex), ceil(spacey), ceil(spacez));
sizeOfArray = numel(myArray);
numBytes = numel(myArray) * 4;
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


% numberOfElements = numel(myArray);
% numBytes = numel(myArray) * 4;