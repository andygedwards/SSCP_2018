function AdaptProtocol(q,tCycle,fac_q,fac_tCycle,nAdaptCycle)

global P

% Number of full adaptation cycles
for i=1:nAdaptCycle
    
    % bring to exercise state
    P.Adapt.FunctionName = 'Adapt0P';
    P.General.q0 = fac_q * q;
    P.General.tCycle = fac_tCycle * tCycle;
    P.General.DtSimulation =100 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    P.Adapt.Fast = 1;
    save P P
    
    CircAdaptP;
    
    % adaption at exercise
    P.Adapt.FunctionName = 'AdaptExcP';
    P.General.DtSimulation = 100 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    P.Adapt.Fast = 0;
    save P P
    
    CircAdaptP;
    
    % bring to resting state
    P.Adapt.FunctionName = 'Adapt0P';
    P.General.q0 = q;
    P.General.tCycle = tCycle;
    P.General.DtSimulation =100 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    P.Adapt.Fast = 1;
    save P P
    
    CircAdaptP;

    % Run to steady state (rest)
    
    % adaption at rest
    P.Adapt.FunctionName = 'AdaptRestP';
    P.General.DtSimulation = 100 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    P.Adapt.Fast = 0;
    save P P
    
    CircAdaptP;
    
end

P.Adapt.FunctionName = 'Adapt0P';
P.General.DtSimulation = 5 * tCycle;
P.General.tEnd = P.t( end ) + P.General.DtSimulation;
P.Adapt.Fast = 0;
save P P
    
CircAdaptP;

load PTemp
save P P

