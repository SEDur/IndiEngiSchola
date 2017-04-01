classdef Simulation
    %Class def does here
    properties
        %Sim Properties
        simMethod
        noDims
        problemSize
        fmax
        lx
        ly
        lz
        T
        Consts
        %Sim Boundary
        alphaXP
        alphaXN
        alphaYP
        alphaYN
        alphaZP
        alphaZN
        zXP
        zXN
        zYP
        zYN
        zZP
        zZN
        %Sim Grids
        p
        ux
        uy
        uz
        %Sim Consts
        rX
        rY
        rZ
        pCx
        pCy
        pCz
        uCx
        uCy
        uCz
        %Source Properties
        noSources
        sourceLoc
        sourceFreq
        sourcePhase
        sourceSignal
        sourceSignalType
        %Reciever Properties
        noRecievers
        recieverLocs
        recievers
        %Sim engine
        Engine
        
        
    end

    methods
        function obj = Simulation (simMethod, noDims, problemSize, fmax)
            obj.Consts = Consts;
            obj.lx = problemSize(1);
            if length(problemSize) > 1
                obj.ly = problemSize(2);
            end
            if length(problemSize) > 2
                obj.lz = problemSize(3);
            end
            obj.fmax = fmax;
            obj.sampleRate = 2*fmax;
            obj.simMethod = simMethod;
            
            switch simMethod
                case 'fdtd'
                    switch noDims
                        case 1
                            obj.Engine = FDTD1D;
                        case 2
                            obj.Engine = FDTD2D;
                        case 3
                            obj.Engine = FDTD3D;
                    end
                case 'pstd'
                    switch noDims
                        case 1
                            obj.Engine = PSTD1D;
                        case 2
                            obj.Engine = PSTD2D;
                        case 3
                            obj.Engine = PSTD3D;
                    end
                case 'sfdtd'
                    switch noDims
                        case 1
                            obj.Engine = SFDTD1D;
                        case 2
                            obj.Engine = SFDTD2D;
                        case 3
                            obj.Engine = SFDTD3D;
                    end
            end
        end
        function run(obj)
            obj.Engine.uCacl;
            obj.Engine.uBoundary;
            obj.Engine.pCalc;
            obj.Engine.pBoundary;
            obj.Engine.source;
            obj.Engine.recieve;
        end
    end
        properties (Dependent)
        domainVolume
        sampleRate
        end
        methods
            function value = get.domainVolume(obj)
                switch obj.noDims
                    case 1
                        value = obj.lx;
                    case 2
                        value = obj.lx * obj.ly;
                    case 3
                        value = obj.lx * obj.ly * obj.lz;
                end
            end
            function value = get.sampleRate(obj)
                value = 2 * obj.fmax;
            end
        end
end