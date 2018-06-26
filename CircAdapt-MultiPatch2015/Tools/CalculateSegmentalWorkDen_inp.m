function sarcWorkDen = CalculateSegmentalWorkDen_inp(Ef,Sf)

sarcWorkDen = zeros( length( Ef( 1, : ) ), 1 );

for i=1:length ( Ef( 1, : ) )
    sarcWorkDen( i, 1 ) = -trapz( Ef( :, i ), Sf( :, i ) );
end

return

