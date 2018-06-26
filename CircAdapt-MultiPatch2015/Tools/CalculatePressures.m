function [ dPdtmax, dPdtmin, meanArterialP, meanAtrialP, systVentrP ] = CalculatePressures
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global P;

ventrPressures = GetFt( 'Node', 'p', { 'Lv', 'Rv' } );
atriaPressures = GetFt( 'Node', 'p', { 'La', 'Ra' } );
arterPressures = GetFt( 'Node', 'p', { 'SyArt', 'PuArt' } );

Dt = P.General.Dt;
pFac = 133.322;

dPdtmax = max(diff(ventrPressures))/pFac/Dt;
dPdtmin = min(diff(ventrPressures))/pFac/Dt;

meanAtrialP = mean(atriaPressures)/pFac;
meanArterialP = mean(arterPressures)/pFac;

systVentrP = max(ventrPressures)/pFac;

end