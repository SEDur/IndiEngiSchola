[p, ux, uy, meta] = setupSim(6000, 0.1);
meta = setupDomain(meta, p, ux, uy);
meta = makeMap(meta, p, ux, uy);
yRec = ceil(25/meta.gx);
xRec = ceil(5/meta.gy);
ct = 1;
t = 0;

while t < meta.T
    [p, ux, uy] = fun2do(p, ux, uy, meta, ct);
    %     [p, ux, uy] = fun2d(p, ux, uy, meta, ct);
    t = t + meta.dt;
    ct = ct + 1;
    myRecieverData(ct) = p(yRec, xRec);
    stepPlot(p, meta, t);
end

myTimeVector = 0 : meta.dt : (length(myRecieverData) - 1) * meta.dt;
plot(myTimeVector, myRecieverData);
axis tight
grid

function [p, ux, uy, meta] = setupSim(fmax, T)

%Units

%%Distance
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
%Time
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;

%% constants
c     = 343 * meters / seconds; %Speed of sound m/s
rho    = 1.21; %Density of air kg/m^3
p0 = 2*10^-5;
% cstab = sqrt(1/2);
cstab = 2/(pi*sqrt(2));
% cstab = 1;

%%
%%Hard Code Variables
%Maximum calculation frequency
% fmax = 10000 * hertz;

%grid size
gx = c * (1/fmax) / cstab;
gy = c * (1/fmax) / cstab;

meta.gx = gx;
meta.gy = gy;

%Dims
%Dim Size (m)
lx = 30*meters;
ly = 30*meters;
xcells = ceil(lx/gx);
ycells = ceil(lx/gy);

%Boundary Absorption Coefs (0 to 1)
alphaL = 1.0;
alphaR = 1.0;
alphaF = 1.0;
alphaB = 1.0;

alphaPositiveX = 0.001;
alphaNegativeX = 0.001;
alphaNegativeY = 0.001;
alphaPositiveY = 0.001;

%number of sources
snum = 2;

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/2)];

%Time of sim
dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))));
meta.dt = dt;
meta.T = T;

% generate the source(s) & determine number of time steps needed
tnum = ceil(T/dt);
source1 = zeros(1,tnum);
source2 = zeros(1,tnum);

stimulus = create_stimulus('tone', 1e-4, 60, T, 0, [], 1/dt);
windowWidth = ceil(length(stimulus)/8);
stimulus = stimulus .* hannrectwindow(1, windowWidth, length(stimulus) - windowWidth, length(stimulus));
gradDelay = ceil(4.17e-3 / dt);
meta.signals = [stimulus,...
    -stimulus*0.5];
    
%     zeros(length(stimulus),1)];

% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);

pObs = zeros(ycells - 1, xcells - 1);
uxObs = zeros(ycells - 1, xcells);
uyObs = zeros(ycells, xcells - 1);

% set up the multiplication constants for the update equations
meta.uCx = ones(size(p)) .* dt/(gx*rho);
meta.uCy = ones(size(p)) .* dt/(gy*rho);
meta.pCx = c^2*rho*dt/gx;
meta.pCy = c^2*rho*dt/gy;

% characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
ZT = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));

ZPositiveX = rho*c*(1 + sqrt(1 - alphaPositiveX))/(1 - sqrt(1 - alphaPositiveX));
ZNegativeX = rho*c*(1 + sqrt(1 - alphaNegativeX))/(1 - sqrt(1 - alphaNegativeX));
ZNegativeY = rho*c*(1 + sqrt(1 - alphaNegativeY))/(1 - sqrt(1 - alphaNegativeY));
ZPositiveY = rho*c*(1 + sqrt(1 - alphaPositiveY))/(1 - sqrt(1 - alphaPositiveY));


% calulcate the coefficients used for the boundary conditions
Rx = rho*gx/dt;
Ry = rho*gy/dt;

% plot vectors
meta.linex = linspace(0, lx - gx, xcells-1);
meta.liney = linspace(0, ly - gx, ycells-1);

meta.xOneuImpedance = ((Rx - ZL)/(Rx + ZL));
meta.xOnepImpedance = (2/(Rx + ZL));
meta.xEnduImpedance = ((Rx - ZR)/(Rx + ZR));
meta.xEndpImpedance = (2/(Rx + ZR));
meta.yOneuImpedance = ((Ry - ZB)/(Ry + ZB));
meta.yOnepImpedance = (2/(Ry + ZB));
meta.yEnduImpedance = ((Ry - ZT)/(Ry + ZT));
meta.yEndpImpedance = (2/(Ry + ZT));

meta.xPositiveuImpedance = ((Rx - ZPositiveX)/(Rx + ZPositiveX));
meta.xPositivepImpedance = (2/(Rx + ZPositiveX));
meta.xNegativeuImpedance = ((Rx - ZNegativeX)/(Rx + ZNegativeX));
meta.xNegativepImpedance = (2/(Rx + ZNegativeX));
meta.yPositiveuImpedance = ((Ry - ZPositiveY)/(Ry + ZPositiveY));
meta.yPositivepImpedance = (2/(Ry + ZPositiveY));
meta.yNegativeuImpedance = ((Ry - ZNegativeY)/(Ry + ZNegativeY));
meta.yNegativepImpedance = (2/(Ry + ZNegativeY));

dimsx = size(ux)-1;
dimsy = size(uy)-1;

end

function meta = setupDomain(meta, p, ux, uy)

referencePointX = ceil(12.5/meta.gx);
referencePointY = ceil(12.5/meta.gy);

%source locations
s1loc = [referencePointY+ceil((5/meta.gy)), referencePointX+ceil((2.5/meta.gx))];
s2loc = [referencePointY+ceil((5/meta.gy))-1, referencePointX+ceil((2.5/meta.gx))];
sources = [s1loc; s2loc];
meta.sources = sub2ind(size(p), sources(:,1), sources(:,2));

% Row % Col
wallOne = [referencePointY : referencePointY + ceil(5/meta.gy);...
    ones(1, length(referencePointY : referencePointY + ceil(5/meta.gy))) * referencePointX];

wallTwo = [ones(1, length(referencePointX : referencePointX + ceil(5/meta.gx))) * referencePointY;...
    referencePointX : referencePointX + ceil(5/meta.gx)];
% wallTwo = [1; 1];

wallThree = [referencePointY : referencePointY + ceil(5/meta.gy);...
    ones(1, length(referencePointY : referencePointY + ceil(5/meta.gy))) * (referencePointX + ceil(5/meta.gx))];

wallFour = [ones(1, length(referencePointX : referencePointX + ceil(5/meta.gx))) * (referencePointY+ ceil(5/meta.gy));...
    referencePointX : referencePointX + ceil(5/meta.gx)];
wallFour = [1;1];
% wallFoursection1 = [ones(1, length(referencePointX : referencePointX + ceil(2/meta.gx))) * (referencePointY+ ceil(5/meta.gy));...
%     referencePointX : referencePointX + ceil(2/meta.gx)];
% wallFoursection2 = [ones(1, length(referencePointX : referencePointX + ceil(2/meta.gx))) * (referencePointY+ ceil(5/meta.gy));...
%     referencePointX+ ceil(3/meta.gx) : referencePointX + ceil(5/meta.gx)];
% wallFour = [wallFoursection1 wallFoursection2];

meta.negativeXWalls = sub2ind(size(p), wallOne(1,:), wallOne(2,:));
meta.negativeYWalls = sub2ind(size(p), wallTwo(1,:), wallTwo(2,:));
meta.positiveXWalls = sub2ind(size(p), wallThree(1,:), wallThree(2,:));
meta.positiveYWalls = sub2ind(size(p), wallFour(1,:), wallFour(2,:));

%     meta.deadSpace = poly2mask([meta.negativeXWalls meta.positiveXWalls], [meta.negativeYWalls meta.positiveYWalls], size(p,1), size(p,2));

meta.walls = [meta.negativeXWalls meta.negativeYWalls meta.positiveXWalls meta.positiveYWalls];

meta.pMask = ones(size(p));
meta.pMask(meta.walls) = zeros;

meta.uxMask = zeros(size(ux));
meta.uxMask(:, 1:end-1) = abs(meta.pMask - 1);
% meta.uyMask = ones(size(uy,2), size(uy,1));
meta.uyMask = zeros(size(uy));
meta.uyMask(1:end-1, :) = abs(meta.pMask - 1);


end

function meta = makeMap(meta, p, ux, uy)

referencePointX = ceil(5/meta.gx);
referencePointY = ceil(5/meta.gy);

velocityRectangle= drawRectangle([referencePointX , referencePointY], [5, 5]);
meta.pressureRectangle = drawRectangle([referencePointX , referencePointY], [5, 5]);
meta.xWallsNeg = sub2ind(size(ux), velocityRectangle(4,:),velocityRectangle(3,:));
meta.yWallsNeg = sub2ind(size(uy), velocityRectangle(2,:),velocityRectangle(1,:));
pxWallsNeg = sub2ind(size(p), meta.pressureRectangle(4,:),meta.pressureRectangle(3,:));
pyWallsNeg = sub2ind(size(p), meta.pressureRectangle(2,:),meta.pressureRectangle(1,:));

mask = zeros(size(p));
mask(pxWallsNeg) = ones;
mask(pyWallsNeg) = ones;

uxNcounter = 1;
uxPcounter = 1;
uyNcounter = 1;
uyPcounter = 1;
uxDiffCounter = 1;
uyDiffCounter = 1;

for iYcoord = 1 : size(mask, 1) - 1
    for iXcoord = 1 : size(mask, 2) - 1
        % Check to see if this value is 1, if so, handle the north and east
        % values
        if mask(iYcoord, iXcoord) >= 1
            if ((mask(iYcoord, iXcoord + 1) < 1) && ((iXcoord + 1) < size(mask,2)))
                uxPosBoundariesCoord(uxPcounter, :) = [iYcoord, iXcoord + 1];
                uxPosBoundariesPCoord(uxPcounter, :) = [iYcoord, iXcoord + 1];
                uxPosBoundariesImpedace(uxPcounter) = mask(iYcoord, iXcoord);
                uxPcounter = uxPcounter + 1;
            end
            if ((mask(iYcoord + 1, iXcoord) < 1) && ((iYcoord + 1) < size(mask,1)))
                uyPosBoundariesCoord(uyPcounter, :) = [iYcoord + 1, iXcoord];
                uyPosBoundariesPCoord(uyPcounter, :) = [iYcoord + 1, iXcoord];
                uyPosBoundariesImpedace(uyPcounter) = mask(iYcoord, iXcoord);
                uyPcounter = uyPcounter + 1;
            end
        else
            % This value is 0, check east north west and south
            % Look west
            if mask(iYcoord, iXcoord + 1) >= 1
                uxNegBoundariesCoord(uxNcounter, :) = [iYcoord, iXcoord + 1];
                uxNegBoundariesPCoord(uxNcounter, :) = [iYcoord, iXcoord];
                uxNegBoundariesImpedace(uxNcounter) = mask(iYcoord, iXcoord + 1);
                uxNcounter = uxNcounter + 1;
            elseif  ((mask(iYcoord, iXcoord + 1) == 0) && ((iXcoord + 1) < size(mask,2)))
                uxDiffCoord(uxDiffCounter, :) = [iYcoord, iXcoord + 1];
                uxDiffPCoords(uxDiffCounter, :) = [[iYcoord, iXcoord] [iYcoord, iXcoord + 1]];
                uxDiffCounter = uxDiffCounter + 1;
            end
            % Look North
            if mask(iYcoord + 1, iXcoord) >= 1
                uyNegBoundariesCoord(uyNcounter, :) = [iYcoord + 1, iXcoord];
                uyNegBoundariesPCoord(uyNcounter, :) = [iYcoord, iXcoord];
                uyNegBoundariesImpedace(uyNcounter) = mask(iYcoord+ 1, iXcoord);
                uyNcounter = uyNcounter + 1;
            elseif  ((mask(iYcoord + 1, iXcoord) == 0) && ((iYcoord + 1) < size(mask, 1)))
                uyDiffCoord(uyDiffCounter, :) = [iYcoord, iXcoord + 1];
                uyDiffPCoords(uyDiffCounter, :) = [[iYcoord, iXcoord] [iYcoord + 1, iXcoord]];
                uyDiffCounter = uyDiffCounter + 1;
            end
        end
    end
end

                meta.uxPosBoundariesCoord = uxPosBoundariesCoord;
                meta.uxPosBoundariesPCoord = uxPosBoundariesPCoord;
                meta.uxPosBoundariesImpedace = uxPosBoundariesImpedace;
                meta.uyPosBoundariesCoord  = uyPosBoundariesCoord;
                meta.uyPosBoundariesPCoord = uyPosBoundariesPCoord;
                meta.uyPosBoundariesImpedace = uyPosBoundariesImpedace;
                meta.uxNegBoundariesCoord = uxNegBoundariesCoord;
                meta.uxNegBoundariesPCoord = uxNegBoundariesPCoord;
                meta.uxNegBoundariesImpedace = uxNegBoundariesImpedace;
                meta.uxDiffCoord = uxDiffCoord;
                meta.uxDiffPCoords = uxDiffPCoords;
                meta.uyNegBoundariesCoord = uyNegBoundariesCoord;
                meta.uyNegBoundariesPCoord = uyNegBoundariesPCoord;
                meta.uyNegBoundariesImpedace = uyNegBoundariesImpedace ;
                meta.uyDiffCoord = uyDiffCoord;
                meta.uyDiffPCoords = uyDiffPCoords;

end

function newCoords = translate(shape, p)

end

function [p, ux, uy] = fun2do(p, ux, uy, meta, ct)
% update the non-boundary condition nodes for velocity
% ux(:, 2:end-1) = ux(:, 2:end-1) - meta.uCx*(p(:, 2:end) - p(:, 1:end-1));
% uy(2:end-1, :) = uy(2:end-1, :) - meta.uCy*(p(2:end, :) - p(1:end-1, :));

ux(:, 2:end-1) = ux(:, 2:end-1) - meta.uCx(:, 2:end) .*(p(:, 2:end) - p(:, 1:end-1));
uy(2:end-1, :) = uy(2:end-1, :) - meta.uCy(2:end, :) .*(p(2:end, :) - p(1:end-1, :));

% update the velocity at the right wall
ux(:, end) = meta.xEnduImpedance * ux(:, end) ...
    + meta.xEndpImpedance * p(:, end);

%update the velocity at the left wall
ux(:, 1) = meta.xOneuImpedance * ux(:, 1) ...
    - meta.xOnepImpedance * p(:, 1);

%update the velocity at the top wall
uy(end, :) = meta.yEnduImpedance * uy(end, :) ...
    + meta.yEndpImpedance * p(end, :);

%update the velocity at the bottom wall
uy(1, :) = meta.yOneuImpedance * uy(1, :) ...
    - meta.yOnepImpedance * p(1, :);

% update the pressure at all nodes
p = meta.pMask .* (p - meta.pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
    - meta.pCy*(uy(2:end, :) - uy(1:end-1, :)));

% p(meta.s1loc(1),meta.s1loc(2)) = p(meta.s1loc(1),meta.s1loc(2)) - meta.source1(ct);
% p(meta.s2loc(1),meta.s2loc(2)) = p(meta.s2loc(1),meta.s2loc(2)) - meta.source2(ct);
p(meta.sources) = p(meta.sources) - meta.signals(ct, :)';

end

function [p, ux, uy] = fun2d(p, ux, uy, meta, ct)
% update the non-boundary condition nodes for velocity
ux(:, 2:end-1) = ux(:, 2:end-1) - meta.uCx*(p(:, 2:end) - p(:, 1:end-1));
uy(2:end-1, :) = uy(2:end-1, :) - meta.uCy*(p(2:end, :) - p(1:end-1, :));

% update the velocity at the right wall
ux(:, end) = meta.xEnduImpedance * ux(:, end) ...
    + meta.xEndpImpedance * p(:, end);

%update the velocity at the left wall
ux(:, 1) = meta.xOneuImpedance * ux(:, 1) ...
    - meta.xOnepImpedance * p(:, 1);

%update the velocity at the top wall
uy(end, :) = meta.yEnduImpedance * uy(end, :) ...
    + meta.yEndpImpedance * p(end, :);

%update the velocity at the bottom wall
uy(1, :) = meta.yOneuImpedance * uy(1, :) ...
    - meta.yOnepImpedance * p(1, :);

% update the pressure at all nodes
p = p - meta.pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
    - meta.pCy*(uy(2:end, :) - uy(1:end-1, :));

p(meta.s1loc(1),meta.s1loc(2)) = p(meta.s1loc(1),meta.s1loc(2)) - meta.source1(ct);
p(meta.s2loc(1),meta.s2loc(2)) = p(meta.s2loc(1),meta.s2loc(2)) - meta.source2(ct);
end

function stepPlot(p, meta, t)
surf(meta.linex, meta.liney, p);
colormap('winter');
shading interp;
title(sprintf('Time = %.6f s',t),...
    'Color',[0 0 0],'FontSize', 14);
xlabel('Width (meters)', 'Color', [0 0 0]);
ylabel('Length (meters)', 'Color', [0 0 0]);
colorbar;
caxis([-1e-6 1e-6]);
zlim([-1e-6 1e-6]);
axis equal;
view(2);
drawnow();
end

function [deadspace] = findDeadSpace(walls, domain)
% essentially performs a floodfill to find dead space  and set its
% value to zero
mask = zeros(size(domain));
mask(walls) = 0;

for  arrayYIndex = 1 : size(domain, 1)
    for arrayXIndex = 1 : size(domain, 2)
        % get the next index
        newIndex = mask(arrayYIndex, arrayXIndex);
        
    end
end

end

function circle = drawCircle(referenceCoordinate, diameter)

end

function rectangle = drawRectangle(referenceCoordinate, edgeLengths)
% rectangle = drawRectangle(referenceCoordinate, edgeLengths)
% rectangle = drawRectangle([1,1], [5,5])
% Draw a rectangle from the bottom left corner, going clockwise. This
% returns a two column array of [yCoordinates xCoordinates] in the
% following order:  [westWall northWall eastWall southWall]

referencePointY = referenceCoordinate(1);
referencePointX = referenceCoordinate(2);
westEdgeLength = edgeLengths(1);
northEdgeLength = edgeLengths(2);
eastEdgeLength = edgeLengths(1);
southEdgeLength = edgeLengths(2);

westWall = [referencePointY : referencePointY + westEdgeLength;...
    ones(1, length(referencePointY : referencePointY + westEdgeLength)) * referencePointX];

northWall = [ones(1, length(referencePointX : referencePointX + northEdgeLength)) * referencePointY;...
    referencePointX : referencePointX + northEdgeLength];

eastWall = [referencePointY : referencePointY + eastEdgeLength;...
    ones(1, length(referencePointY : referencePointY + eastEdgeLength)) * (referencePointX + northEdgeLength)];

southWall = [ones(1, length(referencePointX : referencePointX + southEdgeLength)) * (referencePointY + eastEdgeLength);...
    referencePointX : referencePointX + southEdgeLength];

rectangle = [[westWall eastWall]; [northWall southWall]];

end

function lineSegment = drawLineSegment(startCoordinate, endCoordinate)

% wallOne = [referencePointY : referencePointY + ceil(5/meta.gy);...
%     ones(1, length(referencePointY : referencePointY + ceil(5/meta.gy))) * referencePointX];

% Ycoords; Xcoords
lineSegment = [startCoordinate(1) : endCoordinate(1);...
    startCoordinate(2) : endCoordinate(2)];

end