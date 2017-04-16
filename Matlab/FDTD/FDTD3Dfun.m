function [p, ux, uy, uz] = FDTD3Dfun(p, pCx, pCy, pCz, ux, uy, uz, uCx,...
    uCy, uCz, Rx, Ry, Rz, ZxN, ZxP, ZyN, ZyP, ZzN, ZzP)
% Function that performs one timestep of FDTD method for acoustic simulation.
% 
% This function performs central finite difference calculations on
% matricies that represent pressure and velocity. This function assumes
% that a linear acoustic wave equation is being solved, and so assumes that
% the velocity terms are orthoganal and there are no cross-terms. This
% function solves empirical semi-absorbing boundary conditions, using the
% acoustic impedance of the boundary based on a normalised aproximation of
% absorption coefficient.
% 
% Takes the following arguments:
% p = N:N:N matrix of pressure values
% ux = N:N+1:N matrix of velocity values
% uy = N+1:N:N matrix of velocity values
% uz = N:N:N+1 matrix of velocity values
% pCx = constant related to pressure calculation in x direction
% pCy = constant related to pressure calculation in y direction
% pCz = constant related to pressure calculation in z direction
% uCx = constant related to velocity calculation in x direction
% uCy = constant related to velocity calculation in y direction
% uCz = constant related to velocity calculation in z direction
% Rx = (rho0*dx)/(0.5*dt) Constant related to field constants
% Ry = (rho0*dy)/(0.5*dt) Constant related to field constants
% Rz = (rho0*dz)/(0.5*dt) Constant related to field constants
% ZxN = acoutsitc impedance term at boundary in -x direction
% ZxP = acoutsitc impedance term at boundary in +x direction
% ZyN = acoutsitc impedance term at boundary in -y direction
% ZyP = acoutsitc impedance term at boundary in +y direction
% ZzN = acoutsitc impedance term at boundary in -z direction
% ZzP = acoutsitc impedance term at boundary in +z direction
%
% This functions returns the pressure and velocity field matricies
%

    % Calculate central difference aproximation to velocity field 
    % Velocity in a direction at current timestep excluding the boundarys 
    % = velocity 1 time step ago - constants * pressure 
    % differential half a time step ago in that direction
    ux(:, 2:end-1, :) = ux(:, 2:end-1,:) - uCx*(p(:, 2:end,:) - p(:, 1:end-1, :));
    uy(2:end-1, :, :) = uy(2:end-1, :, :) - uCy*(p(2:end, :, :) - p(1:end-1, :, :));
    uz(:, :, 2:end-1) = uz(:, :, 2:end-1) - uCz*(p(:, :, 2:end) - p(:, :, 1:end-1));

    % update the velocity at the negative x boundary
    % Velocity at this boundary for all of y and z = time and space step
    % normalised by the lovel impedance condition * current velocity values
    % - 2 / time and space discretization * local pressure value 
    ux(:, 1, :) = ((Rx - ZxN)/(Rx + ZxN))*ux(:, 1, :)...
        - (2/(Rx + ZxN))*p(:, 1, :);
    
    % update the velocity at the positive x boundary
    ux(:, end, :) = ((Rx - ZxP)/(Rx + ZxP))*ux(:, end, :) ...
        + (2/(Rx + ZxP))*p(:, end, :);

    % update the velocity at the negative y boundary
    uy(1, :, :) = ((Ry - ZyP)/(Ry + ZyP))*uy(1, :, :)...
        - (2/(Ry + ZyP))*p(1, :, :);
    
    % update the velocity at the positive y boundary
    uy(end, :, :) = ((Ry - ZyN)/(Ry + ZyN))*uy(end, :, :) ...
        + (2/(Ry + ZzN))*p(end, :, :);
    
    % update the velocity at the negative z boundary
    uz(:, :, 1) = ((Rz - ZzP)/(Rz + ZzP))*uz(:, :, 1) - ...
    (2/(Rz + ZzP))*p(:, :, 1);

    % update the velocity at the positive z boundary
    uz(:, :, end) = ((Rz - ZzN)/(Rz + ZzN))*uz(:, :, end) + ...
    (2/(Rz + ZzN))*p(:, :, end);

    % update the pressure at all nodes
    % new pressure across domain = pressure across domain 1 time step ago - 
    % (space,time and wave speed constant) * central difference of 
    % velocities half a time step  ago in all three dimensions
    p = p - pCx*(ux(:, 2:end, :) - ux(:, 1:end-1, :))...
        - pCy*(uy(2:end, :, :) - uy(1:end-1, :, :))...
        - pCz*(uz(:, :, 2:end) - uz(:, :, 1:end-1));
end