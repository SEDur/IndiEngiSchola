classdef Initz
    %Class def does here
    properties
        %Sim Properties
        simMethod
        noDims
        problemSize
        fmax
        dx
        dy
        dz
        T
        dt
        Consts
        sampleRate
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
        function obj = Initz (simMethod, noDims, problemSize, fmax)
            obj.Consts = Consts;
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
    end
end