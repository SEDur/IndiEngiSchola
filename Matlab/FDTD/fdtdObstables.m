[p, ux, uy, meta] = setupSim();
meta = setupObstacles(meta, p, ux, uy);

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

function [p, ux, uy, meta] = setupSim()

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
fmax = 6000 * hertz;

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
%source locations
% meta.s1loc = [ceil((15/gx)) ceil((15/gx)/2)];
% meta.s2loc = [(meta.s1loc(1) - ceil(1.05/gx)) (meta.s1loc(2) + ceil(1.05/gx))];%source frequency
meta.s1loc = [ceil((15/gx)) ceil((15/gy))];
% meta.s2loc = [(meta.s1loc(1)) (meta.s1loc(2) + ceil(1.05/gy))];%source frequency
meta.s2loc = [(meta.s1loc(1)) (meta.s1loc(2) + ceil(sqrt(1.05^2 + 1.05^2)/gy))];%source frequency

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/2)];

%Time of sim
dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))));
meta.dt = dt;
T = 0.2;
meta.T = T;

% generate the source(s) & determine number of time steps needed
tnum = ceil(T/dt);
source1 = zeros(1,tnum);
source2 = zeros(1,tnum);

stimulus = create_stimulus('tone', 1e-4, 100, T, 0, [], 1/dt);
stimulus = stimulus .* hannrectwindow(1, 200 , length(stimulus) - 200, length(stimulus));
gradDelay = ceil(4.17e-3 / dt);
meta.source1 = stimulus;
meta.source2 = zeros(length(stimulus),1);
% meta.source2 = [zeros(gradDelay,1); -stimulus];

% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);

pObs = zeros(ycells - 1, xcells - 1);
uxObs = zeros(ycells - 1, xcells);
uyObs = zeros(ycells, xcells - 1);

% set up the multiplication constants for the update equations
meta.uCx = dt/(gx*rho);
meta.uCy = dt/(gy*rho);
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

function meta = setupObstacles(meta, p, ux, uy)

referencePointX = ceil(20/meta.gx);
referencePointY = ceil(5/meta.gy);

% Row % Col
wallOne = [referencePointY : referencePointY + ceil(5/meta.gy);...
    ones(1, length(referencePointY : referencePointY + ceil(5/meta.gy))) * referencePointX];

wallTwo = [ones(1, length(referencePointX : referencePointX + ceil(5/meta.gx))) * referencePointY;...
    referencePointX : referencePointX + ceil(5/meta.gx)];

wallThree = [referencePointY : referencePointY + ceil(5/meta.gy);...
    ones(1, length(referencePointY : referencePointY + ceil(5/meta.gy))) * (referencePointX + ceil(5/meta.gx))];

wallFour = [ones(1, length(referencePointX : referencePointX + ceil(5/meta.gx))) * (referencePointY+ ceil(5/meta.gy));...
    referencePointX : referencePointX + ceil(5/meta.gx)];

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

% template = zeros(size(ux));
% template(:, 2:end-1) = ones;
% template(meta.walls) = zeros;
% % meta.uXBiasTemplate = find(template(:, 2:end-1));
% meta.uXBiasTemplate = find(template);
% 
% 
% template = zeros(size(ux));
% template(:, 2:end) = ones;
% template(meta.walls) = zeros;
% % meta.uXNegBiasTemplate = find(template(:, 2:end));
% meta.uXNegBiasTemplate = find(template);
% 
% template = zeros(size(ux));
% template(:, 1:end-1) = ones;
% template(meta.walls) = zeros;
% % meta.uXPosBiasTemplate = find(template(:, 1:end-1));
% meta.uXPosBiasTemplate = find(template);
% 
% template = zeros(size(uy));
% template(2:end-1, :) = ones;
% template(meta.walls) = zeros;
% % meta.uYBiasTemplate = find(template(2:end-1, :));
% meta.uYBiasTemplate = find(template);
% 
% template = zeros(size(uy));
% template(2:end, :) = ones;
% template(meta.walls) = zeros;
% % meta.uYNegBiasTemplate = find(template(2:end, :));
% meta.uYNegBiasTemplate = find(template);
% 
% template = zeros(size(uy));
% template(1:end-1, :) = ones;
% template(meta.walls) = zeros;
% % meta.uYPosBiasTemplate = find(template(1:end-1, :));
% meta.uYPosBiasTemplate = find(template);
% 
% template = ones(size(p));
% template(meta.walls) = zeros;
% meta.pTemplate = find(template);
% 
% template = zeros(size(p));
% template(2:end, :) = ones;
% template(meta.walls) = zeros;
% % meta.PuYNegBiasTemplate = find(template(2:end, :));
% meta.PuYNegBiasTemplate = find(template);
% 
% template = zeros(size(p));
% template(1:end-1, :) = ones;
% template(meta.walls) = zeros;
% % meta.PuYPosBiasTemplate = find(template(1:end-1, :));
% meta.PuYPosBiasTemplate = find(template);
% 
% template = zeros(size(p));
% template(:, 2:end) = ones;
% template(meta.walls) = zeros;
% % meta.PuXNegBiasTemplate = find(template(:, 2:end));
% meta.PuXNegBiasTemplate = find(template);
% 
% template = zeros(size(p));
% template(:, 1:end-1) = ones;
% template(meta.walls) = zeros;
% % meta.PuXPosBiasTemplate = find(template(:, 1:end-1));
% meta.PuXPosBiasTemplate = find(template);



%     boxBottomLeftCorner1 = [ceil(5/gx) ceil(5/gx)];
%     boxWidth1 = ceil(20/gx);
%     boxDepth1 = ceil(2/gx);
%
%     boxBottomLeftCorner2 = [ceil(5/gx) ceil(25/gx)];
%     boxWidth2 = ceil(2/gx);
%     boxDepth2 = ceil(20/gx);
%
%     meta.uxIndexes1 = [boxBottomLeftCorner1(1), boxBottomLeftCorner1(1) + boxDepth1;...
%         boxBottomLeftCorner1(2), (boxBottomLeftCorner1(2) + boxWidth1)-1];
%     meta.uyIndexes1 = [boxBottomLeftCorner1(1), (boxBottomLeftCorner1(1) + boxDepth1)-1;...
%         boxBottomLeftCorner1(2), boxBottomLeftCorner1(2) + boxWidth1];
%
%      meta.uxIndexes2 = [boxBottomLeftCorner2(1), boxBottomLeftCorner2(1) + boxDepth2;...
%         boxBottomLeftCorner2(2), (boxBottomLeftCorner2(2) + boxWidth2)-1];
%     meta.uyIndexes2 = [boxBottomLeftCorner2(1), (boxBottomLeftCorner2(1) + boxDepth2)-1;...
%         boxBottomLeftCorner2(2), boxBottomLeftCorner2(2) + boxWidth2];

end

function [p, ux, uy] = fun2do(p, ux, uy, meta, ct)
% update the non-boundary condition nodes for velocity
    ux(:, 2:end-1) = ux(:, 2:end-1) - meta.uCx*(p(:, 2:end) - p(:, 1:end-1));
    uy(2:end-1, :) = uy(2:end-1, :) - meta.uCy*(p(2:end, :) - p(1:end-1, :));
%     ux(:, 2:end-1) = meta.uxMask(:, 2:end-1) .* ux(:, 2:end-1) - meta.uCx*(p(:, 2:end) - p(:, 1:end-1));
%     uy(2:end-1, :) = meta.uyMask(2:end-1, :) .* uy(2:end-1, :) - meta.uCy*(p(2:end, :) - p(1:end-1, :));
% ux(meta.uXBiasTemplate) = ux(meta.uXBiasTemplate)...
%     - meta.uCx*(p(meta.PuXNegBiasTemplate) - p(meta.PuXPosBiasTemplate));
% uy(meta.uYBiasTemplate) = uy(meta.uYBiasTemplate)...
%     - meta.uCy*(p(meta.PuYNegBiasTemplate) - p(meta.PuYPosBiasTemplate));

% update the velocity at the right wall
ux(:, end) = meta.xEnduImpedance * ux(:, end) ...
    + meta.xEndpImpedance * p(:, end);

% ux(meta.negativeXWalls) = meta.xNegativeuImpedance * ux(meta.negativeXWalls) ...
%     + meta.xNegativepImpedance * p(meta.negativeXWalls);
%     ux(meta.walls) = meta.xNegativeuImpedance * ux(meta.walls) ...
%         + meta.xNegativepImpedance * p(meta.walls);

%update the velocity at the left wall
ux(:, 1) = meta.xOneuImpedance * ux(:, 1) ...
    - meta.xOnepImpedance * p(:, 1);

% ux(meta.positiveXWalls) = meta.xPositiveuImpedance * ux(meta.positiveXWalls) ...
%     + meta.xPositivepImpedance * p(meta.positiveXWalls);
%     x(meta.walls) = meta.xPositiveuImpedance * ux(meta.walls) ...
%         + meta.xPositivepImpedance * p(meta.walls);

%update the velocity at the top wall
uy(end, :) = meta.yEnduImpedance * uy(end, :) ...
    + meta.yEndpImpedance * p(end, :);

% uy(meta.negativeYWalls) = meta.yNegativeuImpedance * ux(meta.negativeYWalls) ...
%     + meta.yNegativepImpedance * p(meta.negativeYWalls);
%     uy(meta.walls) = meta.yNegativeuImpedance * ux(meta.walls) ...
%         + meta.yNegativepImpedance * p(meta.walls);

%update the velocity at the bottom wall
uy(1, :) = meta.yOneuImpedance * uy(1, :) ...
    - meta.yOnepImpedance * p(1, :);

% uy(meta.positiveYWalls) = meta.yPositiveuImpedance * uy(meta.positiveYWalls) ...
%     - meta.yPositivepImpedance * p(meta.positiveYWalls);
%     uy(meta.walls) = meta.yPositiveuImpedance * uy(meta.walls) ...
%         - meta.yPositivepImpedance * p(meta.walls);

% 
% ux(meta.walls) = zeros;
% 
% uy(meta.walls) = zeros;

% update the pressure at all nodes
p = meta.pMask .* (p - meta.pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
    - meta.pCy*(uy(2:end, :) - uy(1:end-1, :)));
% p = meta.pMask .* p - meta.pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
%     - meta.pCy*(uy(2:end, :) - uy(1:end-1, :));

% p(meta.pTemplate) = p(meta.pTemplate) - meta.pCx*(ux(meta.uXNegBiasTemplate) - ux(meta.uXPosBiasTemplate))...
%     - meta.pCy*(uy(meta.uYNegBiasTemplate) - uy(meta.uYPosBiasTemplate));

% meta.PuXNegBiasTemplate = find(template(:, 2:end));
% meta.PuXPosBiasTemplate = find(template(:, 1:end-1));
% meta.PuYNegBiasTemplate = find(template(2:end, :));
% meta.PuYPosBiasTemplate = find(template(1:end-1, :));

p(meta.s1loc(1),meta.s1loc(2)) = p(meta.s1loc(1),meta.s1loc(2)) - meta.source1(ct);
p(meta.s2loc(1),meta.s2loc(2)) = p(meta.s2loc(1),meta.s2loc(2)) - meta.source2(ct);

%     ux(meta.uxIndexes1(1,1):meta.uxIndexes1(1,2),...
%         meta.uxIndexes1(2,1):meta.uxIndexes1(2,2)) =...
%         ux(meta.uxIndexes1(1,1):meta.uxIndexes1(1,2),...
%         meta.uxIndexes1(2,1):meta.uxIndexes1(2,2)) .* 0.0002;
%
%     uy(meta.uyIndexes1(1,1):meta.uyIndexes1(1,2),...
%         meta.uyIndexes1(2,1):meta.uyIndexes1(2,2)) =...
%         uy(meta.uyIndexes1(1,1):meta.uyIndexes1(1,2),...
%         meta.uyIndexes1(2,1):meta.uyIndexes1(2,2)) .* 0.0002;
%
%     ux(meta.uxIndexes2(1,1):meta.uxIndexes2(1,2),...
%         meta.uxIndexes2(2,1):meta.uxIndexes2(2,2)) =...
%         ux(meta.uxIndexes2(1,1):meta.uxIndexes2(1,2),...
%         meta.uxIndexes2(2,1):meta.uxIndexes2(2,2)) .* 0.0002;
%
%     uy(meta.uyIndexes2(1,1):meta.uyIndexes2(1,2),...
%         meta.uyIndexes2(2,1):meta.uyIndexes2(2,2)) =...
%         uy(meta.uxIndexes2(1,1):meta.uyIndexes2(1,2),...
%         meta.uyIndexes2(2,1):meta.uyIndexes2(2,2)) .* 0.0002;
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