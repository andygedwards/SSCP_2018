% Number of iterations
numIterations = 10;

% Adapts the currently loaded P file.

tCycleAtRest = P.General.tCycle;
flowAtRest = P.General.q0;

for it = 1 : numIterations
    % Go to exercise
    P.General.tCycle = 0.5 * tCycleAtRest;
    P.General.q0     = 3 * flowAtRest;
    P.General.DtSimulation = 30;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    P.Adapt.Fast = 1;
    CircAdaptP;
    P.Adapt.Fast = 0;
    
    % Adapt at exercise
    P.Adapt.FunctionName = 'AdaptExcP';
    P.General.DtSimulation = 100 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    CircAdaptP;
    
    % Go to rest
    P.Adapt.FunctionName = 'Adapt0P';
    P.General.tCycle = tCycleAtRest;
    P.General.q0     = flowAtRest;
    P.General.DtSimulation = 30;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    P.Adapt.Fast = 1;
    CircAdaptP;
    P.Adapt.Fast = 0;
    
    % Adapt at rest
    P.Adapt.FunctionName = 'AdaptRestP';
    P.General.DtSimulation = 50 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    CircAdaptP;
    
    % Return to normal
    P.Adapt.FunctionName = 'Adapt0P';
    P.General.DtSimulation = 1.5 * P.General.tCycle;
    P.General.tEnd = P.t( end ) + P.General.DtSimulation;
    CircAdaptP;
    
end

save PAdapted P;
