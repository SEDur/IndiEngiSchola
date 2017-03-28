% Create an object of generic parameters
classdef Consts < handle
    %class def goes here
    properties
        %Distance
        meters
        centimeters
        millimeters
        inches
        feet
        %Time
        seconds
        hertz
        kilohertz
        megahertz
        gigahertz
        % constants
        c %Speed of sound m/s
        rho %Density of air kg/m^3
        p0 %Lower limit of pressure audability
        cstab %Stability margine
    end
    methods
        function obj = Consts()
            %Distance
            obj.meters      = 1;
            obj.centimeters = 1e-2 * meters;
            obj.millimeters = 1e-3 * meters;
            obj.inches      = 2.54 * centimeters;
            obj.feet        = 12 * inches;
            
            %Time
            obj.seconds     = 1;
            obj.hertz       = 1/seconds;
            obj.kilohertz   = 1e3 * hertz;
            obj.megahertz   = 1e6 * hertz;
            obj.gigahertz   = 1e9 * hertz;
            
            % constants
            obj.c     = 343 * consts.meters / consts.seconds; %Speed of sound m/s
            obj.rho    = 1.21; %Density of air kg/m^3
            obj.p0 = 10^-12; %Lower limit of pressure audability
            obj.cstab = sqrt(1/3); %Stability margine
        end
    end
end