classdef fdtdDomain
    properties
        p
        ux
        uy
        uz
        dt
    end
    methods
        function obj = fdtdDomain(pdims, uxdims, uydims, uzdims, dt)
            obj.p = zeros(pdims);
            obj.ux = zeros(uxdims);
            obj.uy = zeros(uydims);
            obj.uz = zeros(uzdims);
            obj.dt = dt;
        end
    end
end