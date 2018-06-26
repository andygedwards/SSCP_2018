function [ EDV, ESV, volumeChange, anteAtrVentr, anteVentrArt, retroAtrVentr, retroVentrArt ] = CalculateVentrVolumes
% UNTITLED Summary of this function goes here
% Detailed explanation goes here

global P;

cavityVolumes = GetFt( 'Cavity', 'V', { 'Lv', 'Rv' } ).*1e6;

volumeChange = max( cavityVolumes ) - min( cavityVolumes );

EDV = max( cavityVolumes );
ESV = min( cavityVolumes );

% Calculate flow volumes

aorticFlow = GetFt( 'Valve', 'q', 'LvSyArt' ).*1e6;
aorticEjectionIdxs = find( aorticFlow > 0 );
aorticForwardFlow = zeros( size( aorticFlow ) );
aorticForwardFlow( aorticEjectionIdxs ) = aorticFlow( aorticEjectionIdxs );
aorticRegurgitationIdxs = find( aorticFlow < 0 );
aorticBackwardFlow = zeros( size( aorticFlow ) );
aorticBackwardFlow( aorticRegurgitationIdxs ) = aorticFlow( aorticRegurgitationIdxs );

pulmonaryFlow = GetFt( 'Valve', 'q', 'RvPuArt' ).*1e6;
pulmonaryEjectionIdxs = find( pulmonaryFlow > 0 );
pulmonaryForwardFlow = zeros( size( pulmonaryFlow ) );
pulmonaryForwardFlow( pulmonaryEjectionIdxs ) = pulmonaryFlow( pulmonaryEjectionIdxs );
pulmonaryRegurgitationIdxs = find( pulmonaryFlow < 0 );
pulmonaryBackwardFlow = zeros( size( pulmonaryFlow ) );
pulmonaryBackwardFlow( pulmonaryRegurgitationIdxs ) = pulmonaryFlow( pulmonaryRegurgitationIdxs );

anteVentrArt = [ trapz( P.t, aorticForwardFlow ) trapz( P.t, pulmonaryForwardFlow ) ];
retroVentrArt = [ trapz( P.t, aorticBackwardFlow ) trapz( P.t, pulmonaryBackwardFlow ) ];

mitralFlow = GetFt( 'Valve', 'q', 'LaLv' ).*1e6;
mitralFillingIdxs = find( mitralFlow > 0 );
mitralForwardFlow = zeros( size( mitralFlow ) );
mitralForwardFlow( mitralFillingIdxs ) = mitralFlow( mitralFillingIdxs );
mitralRegurgitationIdxs = find( mitralFlow < 0 );
mitralBackwardFlow = zeros( size( mitralFlow ) );
mitralBackwardFlow( mitralRegurgitationIdxs ) = mitralFlow( mitralRegurgitationIdxs );

tricuspidFlow = GetFt( 'Valve', 'q', 'RaRv' ).*1e6;
tricuspidFillingIdxs = find( tricuspidFlow > 0 );
tricuspidForwardFlow = zeros( size( tricuspidFlow ) );
tricuspidForwardFlow( tricuspidFillingIdxs ) = tricuspidFlow( tricuspidFillingIdxs );
tricuspidRegurgitationIdxs = find( tricuspidFlow < 0 );
tricuspidBackwardFlow = zeros( size( tricuspidFlow ) );
tricuspidBackwardFlow( tricuspidRegurgitationIdxs ) = tricuspidFlow( tricuspidRegurgitationIdxs );

anteAtrVentr = [ trapz( P.t, mitralForwardFlow ) trapz( P.t, tricuspidForwardFlow ) ];
retroAtrVentr = [ trapz( P.t, mitralBackwardFlow ) trapz( P.t, tricuspidBackwardFlow ) ];

end

