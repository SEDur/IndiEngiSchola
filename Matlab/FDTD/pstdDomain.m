classdef pstdDomain
    properties
        pd
        ux
        uy
        uz
        dt
    end
    methods
        function obj = pstdDomain(pd, ux, uy, uz, dt)
            obj.pd = pd;
            obj.ux = ux;
            obj.uy = uy;
            obj.uz = uz;
            obj.dt = dt;
        end
    end
end