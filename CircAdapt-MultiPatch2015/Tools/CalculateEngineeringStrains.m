function [ engStrain, onsetIdx ] = CalculateEngineeringStrains( onsetTime )

global P;

originalTime = P.t - P.t( 1 );

onsetIdx = find( originalTime > onsetTime, 1, 'first' );
sarcLength = P.Patch.Ls;
refSarcLength = sarcLength( onsetIdx, : );

extension = bsxfun( @rdivide, sarcLength, refSarcLength );

engStrain = 100 * ( extension - 1 );

end

