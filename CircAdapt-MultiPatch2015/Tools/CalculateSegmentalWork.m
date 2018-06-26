function sarcWork = CalculateSegmentalWork

global P;

Ef = GetFt( 'Patch', 'Ef', 'All' );
Sf = GetFt( 'Patch', 'Sf', 'All' );
VWall = GetFt( 'Patch', 'VWall', 'All' );

sarcWork = zeros( length( Ef( 1, : ) ), 1 );

for i=1:length ( Ef( 1, : ) )
    sarcWork( i, 1 ) = -trapz( Ef( :, i ), Sf( :, i ) ) * VWall( i );
end

return

