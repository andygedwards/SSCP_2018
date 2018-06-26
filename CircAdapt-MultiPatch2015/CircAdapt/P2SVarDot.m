function P2SVarDot
% function P2SVarDot
% Physiological P -> scaling -> Time derivative of SVar
% stored in P.SVarDot
% time in column direction, SVar in row direction
% Theo Arts, Maastricht University, Oct 30, 2011

global P

ScaleVqY=P.General.ScaleVqY;
ScV=ScaleVqY(1); Scq=ScaleVqY(2); ScY=ScaleVqY(3); Sc0=1;

P.SVarDot= [P.tDot,P.Cavity.VDot/ScV,P.Valve.qDot/Scq,...
    P.Patch.CDot,P.Patch.LsiDot, P.TriSeg.VDot/ScV,P.TriSeg.YDot/ScY];

end

