function TriSegV2p
% function TriSegV2p
% TriSeg is a 3-wall structure (Left,Septal,Right) with 2 cavities (R,L)
% Calculates: cavity volumes V -> dimensions of the 'double bubble',
% myofiber stress Sf, wall tension T and cavity pressures p
% VS and YS repesent septal volume displacement and junction radius.
% State variables P.TriSeg V and Y represent initial estimates to calculate
% VS and YS accurately.
% Theo Arts, Maastricht University, Oct 30, 2011

global P;
TriSeg=P.TriSeg;

n  = TriSeg.n; % number of TriSeg's
ic = TriSeg.iCavity; iCavity=[ic;ic+1]; %related cavities
iw = TriSeg.iWall  ; iWall=[iw;iw+1;iw+2]; %related walls
rhob=P.General.rhob; % blood density

for i=1:n %for all TriSeg's
    iC   = iCavity(:,i); % 2 cavities
    iW   = iWall(:,i); % 3 walls
    Am0  = P.Wall.Am0(:,iW); %zero stress wall area
    DADT = P.Wall.DADT(:,iW); %wall compliance
    nt   = size(Am0,1); %number of time points
    VWall= P.Wall.VWall(iW); % 3 wall volumes
    VWL  = (VWall(1)+VWall(2))/2; % enclosed wall volume 1st cavity
    VWR  = (VWall(3)+VWall(2))/2; % enclosed wall volume 2nd cavity
    V    = max(0,P.Cavity.V(:,iC))+repmat([VWL,VWR],[nt,1]); %2 midwall volumes(t)

    VS=TriSeg.V(:,i); YS= TriSeg.Y(:,i); % 1st estimate of [V,Y]
    VRef=mean(V(:)); % reference volume
    dv=0.01*VRef; dy=0.02*VRef^(1/3); % increments for d[Txy]/dVY

    [Tx0,Ty0]=VY2Txy(VS,YS,Am0,DADT,V); % tension Txy =fu(VY)
    [TxV,TyV]=VY2Txy(VS+dv,YS,Am0,DADT,V); % partial V derivative
    [TxY,TyY]=VY2Txy(VS,YS+dy,Am0,DADT,V); % partial Y derivative
    DTxX=(TxV-Tx0)/dv; DTyX=(TyV-Ty0)/dv; % Jacobian matrix inversion
    DTxY=(TxY-Tx0)/dy; DTyY=(TyY-Ty0)/dy;
    DET= DTxX.*DTyY-DTxY.*DTyX; % determinant
    % 1st iteration to solve TriSeg geometry VY
    dV=(-DTyY.*Tx0+DTxY.*Ty0)./DET;
    dY=(+DTyX.*Tx0-DTxX.*Ty0)./DET;
    VS=VS+dV;
    YS=YS+dY;
    Tau=TriSeg.Tau;
    TriSeg.VDot(:,i)= dV/Tau; % keep track of solution
    TriSeg.YDot(:,i)= dY/Tau;

    for j=1:2 % extra iterations for solution of TriSeg geometry
        [Tx0,Ty0]=VY2Txy(VS,YS,Am0,DADT,V);
        dV=(-DTyY.*Tx0+DTxY.*Ty0)./DET;
        dY=(+DTyX.*Tx0-DTxX.*Ty0)./DET;
        VS=VS+dV;
        YS=YS+dY;
    end
    
    % writing geometric data TriSeg
    TriSeg.YS(:,i)=YS; % final solution Rv-Sv-Lv junction radius
    TriSeg.VS(:,i)=VS; % volume of septal cap
    
    % writing wall data
    [Tx0,Ty0,Am,Cm,T]=VY2Txy(VS,YS,Am0,DADT,V); % solution TriSeg
    P.Wall.Am(:,iWall)= Am; % wall area
    P.Wall.Cm(:,iWall)= Cm; % wall curvature
    P.Wall.T(:,iWall )= T ; % wall tension
    P.Wall.pTrans(:,iWall)= 2*Cm.*T; % transmural pressure
    
    % Cavity impedance properties, needed to make node connection
    Vw= repmat([VWL,VWR],[nt,1]);
    Vm= V + Vw;
    Len= 2*Vm.^(1/3);
    A  = ( V + 0.1*Vw ) ./Len;
    P.Cavity.A(:,iC) = A; % cross-sectional area for valve inflow and outflow pressure
    DADT= P.Wall.DADT(:,iWall([1,3]));
    P.Cavity.Z(:,iCavity) = 0.2*sqrt(rhob*Len./abs(DADT))./A; % Compatibitlity with Tube

end

P.TriSeg=TriSeg;
end

function [Tx,Ty,Am,Cm,Tm]=VY2Txy(VS,YS,Am0,DADT,VLR)
% 1st order approximation of TriSeg according to J. Lumens et al.
% For a wall with zero-stress area Am0 and compliance DADT
% cap volume VS, junction radius YS and cavity volumes VLR=[VL,VR]
% Result: Summed axial and radial tension components [Tx,Ty] on junction
% for 3 walls: Am= midwall area, Cm= midwall curvature, Tm= wall tension
Vm= [VLR,VS]* [...
    -1     0     0
    0     0     1
    1     1     1];
Ym=[YS,YS,YS];

% Solving 3rd order polynomial
% Mathematica: Solve[x^3 + 3y^2x  - 2V == 0, x]
% Q= (V + Sqrt(V^2 + y^6))^(1/3);  x= Q - y^2/Q;
SignVm= sign(Vm); Vm=abs(Vm);
V     = (3/pi)*Vm;
Q     = (V + sqrt(V.^2 + Ym.^6)).^(1/3);
Xm    = SignVm .* ( Q - Ym.^2 ./ Q );

%calculate midwall area Am and curvature Cm=1/rm
X2    = Xm.^2; Y2= Ym.^2;
R2    = X2+Y2;
% Am    = pi*R2; % midwall cap area
Am    = max(Am0,pi*R2); % midwall cap area, buckling with T<0
Cm    = 2*Xm./R2; % midwall cap curvature

% calculation of tension T and components Tx, Ty
Tm=(Am-Am0)./DADT;
Sin= Ym.*Cm;
Cos= (Y2-X2)./R2;
Txi = Cos.*Tm; %
Tyi = Sin.*Tm;
TRef=sqrt(sum(Tm.^2,2)); 
Tx= sum(Txi,2)./TRef; % axial tension component
Ty= sum(Tyi,2)./TRef; % radial tension component
end

