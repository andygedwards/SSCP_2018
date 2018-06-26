function [ COthresh, numIter, COvec, pVec ] = ExerciseCapacity( config )
% EXERCISECAPACITY Determines exercise capacity represented by the maximum CO with
%   mean pressure below a threshold pressure at a given location.

%   The Newton-Raphson method is used to step towards the desired cardiac
%   output at which the critical pressure occurs.
%
%   We aim to solve p( CO )- pCrit = 0, where p( CO ) is the pressure at
%   the desired location.
%   Equivalently, we solve log( p( CO ) ) - log( pCrit ) = 0 to improve
%   convergence of the algorithm.
%
%   Options can be set in the config structure, if supplied:
%   config.pCrit - critical pressure at which peak exercise is reached in
%   mmHg (default 30 mmHg)
%
%   config.location - Location where pressure is measured - must be a node
%   in P.Node.Name (Default 'PuVen')
%
%   config.tolerance - tolerance within which the simulated pressure must
%   be of pCrit (default 1 mmHg)
%
%   config.COInit - Starting point for algorithm in l/min - ideally as
%   close to the critical value as possible (default 15 l/min)
%
%   config.maxIter - maximum number of iterations to perform (default 5).
%
%   config.initStep - initial step size for derivative calculation in l/min
%   (default 1 l/min).
%
%   config.minStepSize - restricts the minimum step sized used when
%   calculating derivatives in l/min (default 0.1 l/min).
%
%   config.stabCrit - tightens the stability tolerances used to determine
%   whether circadapt is stable. Can be useful for preventing lack of
%   convergence when changes in mean pressure become small and unreliable, 
%   leading to negative derivatives that cause divergence from the solution.
%   Default 1, smaller values  = stricter tolerances. Will increase running
%   time dramatically. 
%
%   21 August 2017, Joost Lumens and John Walmsley

global P

% Load default values:
pCrit = 30;
maxIter = 5;
COInit = 15;
initStep = 1;
location = 'PuVen';
tolerance = 1;
minStepSize = 0.1;
stabCrit = 1;

% Remove previous P files to prevent erroneous reloading of files.
if exist( [ pwd '\P.mat' ], 'file' )
    delete P.mat;
end
if exist( [ pwd '\PBackup.mat' ], 'file' )
    delete PBackup.mat;
end
if exist( [ pwd '\PTemp.mat' ], 'file' )
    delete PTemp.mat;
end

% if config is set, overwrite default values:
if nargin == 1
    if isfield( config, 'pCrit' )
        pCrit = config.pCrit;
    end
    if isfield( config, 'maxIter' )
        maxIter = config.maxIter;
    end
    if isfield( config, 'COInit' )
        COInit = config.COInit;
    end
    if isfield( config, 'initStep' )
        initStep = config.initStep;
    end
    if isfield( config, 'location' )
        location = config.location;
    end
    if isfield( config, 'tolerance' )
        tolerance = config.tolerance;
    end
    if isfield( config, 'minStepSize' )
        minStepSize = config.minStepSize;
    end
    if isfield( config, 'stabCrit' )
        stabCrit = config.stabCrit;
    end
end

% In case the simulation crashes at any point, we will need a backup of the file to reload.
save PBackup P;

% initialise the iteration counter
numIter = 1;

% Get the 'resting' pressure flow relation and pCrit (assumes the current P
% structure is viable).
CORest = 60 * 1e3 * P.General.q0; % convert to l/min

% Calculate initial value - take exercise adaptation set point for healthy
% simulation as our beginning value (15 l/min)
[ q0Init, tCycleInit ] = Exercise_COHR_Relation( COInit );
pInit = CalculatePCrit( q0Init, tCycleInit, location, pCrit, stabCrit );

% This point may already be unviable for the P structure (pCrit = NaN)
if isnan( pInit )
    while( isnan( pInit ) ) % pCrit cannot be measured because of simulation crashing
        COInit = CORest + (1/2)*( COInit - CORest ); % reduce the current CO value
        disp( COInit )
        [ q0Init, tCycleInit ] = Exercise_COHR_Relation( COInit ); % Get input values (pressure, flow) from CO.
        pInit = CalculatePCrit( q0Init, tCycleInit, location, pCrit ); % Calculate new pCritInit
    end
end

% Initialise the iterative problem
p_n = pInit;
CO_n = COInit;
diff_n = initStep; % initial step size for derivatives - 10% of distance from current CO to rest v0v.
COvec = COInit;
pVec = p_n;
logP_n = log( p_n ); % logarithms are used to render the pressure - flow relationship pseudolinear, as in most cases it has an exponential form. This improves the convergence of the method.

% Iterate to solution
while ( abs( pCrit - p_n ) > tolerance && numIter < maxIter )
    % Calculate step for derivative:
    if pCrit > p_n % if current pressure too low
        COStep = CO_n + diff_n; % forward difference - use monotonicity here.
        % calculate derivative
        [ q0Step, tCycleStep ] = Exercise_COHR_Relation( COStep ); % Get input values (pressure, flow) from CO.
        pStep = CalculatePCrit( q0Step, tCycleStep, location, pCrit, stabCrit ); % Calculate new pCritInit
        logP_step = log( pStep );
        derivative = ( logP_step - logP_n ) / diff_n;
        
    else % current pressure too high
        COStep = CO_n - diff_n; % backwards difference
        % calculate derivative
        [ q0Step, tCycleStep ] = Exercise_COHR_Relation( COStep ); % Get input values (pressure, flow) from CO.
        pStep = CalculatePCrit( q0Step, tCycleStep, location, pCrit, stabCrit ); % Calculate new pCritInit
        logP_step = log( pStep );
        derivative = ( logP_n - logP_step ) / diff_n;
    end
    % Calculate Newton-Raphson Step:
    DCO = - ( logP_n - log( pCrit ) ) / derivative;
    % update CO_n
    CO_n = CO_n + DCO;
    % Calculate p_n at new CO_n
    [ q0_n, tCycle_n ] = Exercise_COHR_Relation( CO_n );
    p_n = CalculatePCrit( q0_n, tCycle_n, location, pCrit, stabCrit );
    
    % It can happen that the CO is too high and the simulation crashes...
    if isnan( p_n )
        it = 0;
        while( isnan( p_n ) && it < 10 ) % pCrit cannot be measured because of simulation crashing
            disp( 'NaN while calculating p_n' )
            load PBackup % Get the last known good configuration
            CO_n = CO_n - (1/2)* DCO; % reduce the current CO value to a step half the size
            [ q0_n, tCycle_n ] = Exercise_COHR_Relation( CO_n ); % Get input values (pressure, flow) from CO.
            p_n = CalculatePCrit( q0_n, tCycle_n, location, pCrit, stabCrit ); % Calculate new pCritInit
            DCO = 0.5* DCO; % ensures we never get back to the last good point, which is too small
            it = it + 1;
        end
    end
    if isnan( p_n )
        disp( 'Convergence failure. Return NaN' )
        break
    else
        save PBackup P;
    end
    % Update the log fo P_n
    logP_n = log( p_n );
    % update difference - need to make sure it is smaller than the step
    % size or derivative calculation will become unreliable.
    diff_n = min( diff_n, 0.5 * abs( DCO ) );
    % Prevent derivative step getting too small and preventing CircAdapt
    % from differentiating between two cardiac outputs.
    if diff_n < minStepSize
        diff_n = minStepSize;
    end
    % update iterations:
    numIter = numIter + 1;
    COvec = [COvec CO_n ]; % Note: shifted by 1 due to saving initial condition
    pVec = [ pVec p_n ]; % Note: shifted by 1 due to saving initial condition
    
end

% Finalise CO for return
COthresh = CO_n;

% Display message if convergence has not been achieved
if numIter == maxIter
    disp( 'Max. Iterations reached! Convergence may not be adequate' )
    disp( [ 'Error: ' num2str( pCrit - p_n ) ' mmHg' ] )
end

end

