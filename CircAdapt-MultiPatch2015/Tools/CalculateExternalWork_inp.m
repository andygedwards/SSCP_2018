function ExtWork = CalculateExternalWork_inp(V,p)

ExtWork = zeros( length( V( 1, : ) ), 1 );

for i=1:length ( V( 1, : ) )
    ExtWork( i, 1 ) = -trapz( V( :, i ), p( :, i ) );
end

return

