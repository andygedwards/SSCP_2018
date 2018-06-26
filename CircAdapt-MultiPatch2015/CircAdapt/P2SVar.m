function P2SVar

global P

ScaleVqY=P.General.ScaleVqY;
ScV=ScaleVqY(1); Scq=ScaleVqY(2); ScY=ScaleVqY(3);

P.SVar=[P.t,P.Cavity.V/ScV,P.Valve.q/Scq,P.Patch.C,P.Patch.Lsi,...
    P.TriSeg.V/ScV,P.TriSeg.Y/ScY];

end

