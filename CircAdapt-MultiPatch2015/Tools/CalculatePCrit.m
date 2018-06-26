function [ pCrit ] = CalculatePCrit( q0, tCycle, location, pCritThresh, stabCrit )
%CALCULATEPCRIT Calculates pressures at a given location
%   For a flow q0 (m^3/s) and a cycle time tCycle (s), this function
%   attempts to calculate the average pressure over a beat at the node
%   location.
%
%   If the simulation crashes, NaN is returned.

global P

try
    it = 0;
    stationary = 2 / stabCrit;
    P.General.q0 = q0;
    P.General.tCycle = tCycle;
    % Try to reduce the duration of simulations by not running further once
    % stationarity = 0.
    while ( it < 20 && stationary > 0 )
        P.General.DtSimulation = 5;
        P.Adapt.FunctionName = 'Adapt0P'; %No adaptation
        P.General.PressFlowContr = 1;
        if stationary > 2 / stabCrit % Guarantees that the first iteration runs fully
            P.Adapt.Fast = 1;
        else 
            P.Adapt.Fast = 0;
        end
        P.General.tEnd = P.t(end) + P.General.DtSimulation;
        save P P
        
        evalc( 'CircAdaptP' );
        % Calculate stationarity as in Adapt0P
        ErrVec = 1000 * log( P.Adapt.Out( end, : ) ./ P.Adapt.In( end, : ) );
        stationary = round( ( 1 / stabCrit ) * sqrt( mean( ErrVec .^ 2 ) ) );
        it = it + 1;
        
        load( [ pwd '\PTemp.mat' ] )
        pCrit = 0.0075 * mean( GetFt( 'Node', 'p', location ) ); % in mmHg
        % If the pressures are excessive, crash the simulation
        if pCrit > pCritThresh + 15
            disp( [ 'Pressure threshold exceeded: ' num2str( pCrit ) ' mmHg > ' num2str( pCritThresh + 15 ) ' mmHg' ] )
            pCrit = NaN;
            load PBackup;
            break
        end
        load PTemp
        save Paux P
        delete( [ pwd '\PTemp.mat' ] ) % Prevents PTemp files hanging around and being accidentally reloaded.
        
    end
       
catch % If the simulation crashes
    disp( 'Simulation has crashed' )
    pCrit = NaN;
    load PBackup;
end

end

