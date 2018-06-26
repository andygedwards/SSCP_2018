function SarcEf2Sf;
%Theo Arts, Maastricht University. July 30, 2006.

global P; % general time
Sarc= P.Patch;

%==== Input variables
t       = P.t;
Ef      = Sarc.Ef ;
nt      = size(Ef,1); % nt: number of times; nr: number of sarcomeres
Col1    = ones(nt,1);
tc      = Tc(t,Sarc.ActivationDelay);
Lsi     = Sarc.Lsi;
C       = Sarc.C;

LenSeriesElement= SparseDiag(Sarc.LenSeriesElement);
TR      = SparseDiag(Sarc.TR     ); %TR atrial muscle > ventricular muscle
TD      = SparseDiag(Sarc.TD     );
TimeAct = SparseDiag(Sarc.TimeAct); %time scale of contraction pulse
Ls0     = SparseDiag(Sarc.Lsi0Act); %zero active stress sarcomere length
Ls0Pas  = SparseDiag(Sarc.Ls0Pas );
dLsPas  = SparseDiag(Sarc.dLsPas );
SfPas   = SparseDiag(Sarc.SfPas  );
CRest   = SparseDiag(Sarc.CRest  ); %Resting C-value (Ca++  contractility)
LsRef   = SparseDiag(Sarc.LsRef  );
SfAct   = SparseDiag(Sarc.SfAct  );
vMax    = SparseDiag(Sarc.vMax   );

% series elasticity and sarcomere shortening velocity
Ls         = exp(Ef)*LsRef;
Sarc.Ls    = Ls;

%=== Active sarcomere
% constants related to timing are mainly based on experimental findings
L  = max(Lsi/Ls0-1,0.0001) ; % normalized sarc length for active contraction
tA = (0.65+1.0570*L)*TimeAct; % activation time lengthens with sarcomere length
tR = 0.55*TR*TimeAct        ; % rise time
tD = 0.33*TD*TimeAct        ; % decay time (default 0.22)
T  = tc/tR;
x=min(8,max(0,T)); % normalized time during rise of activation
ft1= x.^3 .* exp(-x) .* (8-x).^2 * 0.020 / tR;
%         rise of contraction, 'amount of Ca++ release'
%Integral T^n exp(-T) = Gamma(n+1) = n!
x=(tc-tA)/tD; % normalized time during decay of activation
tanhx= 0.5+0.5*sin( sign(x).*min(pi/2,abs(x)) ); %always>0
% Time confined approximation of 1/(1-e^x) function
FL= tanh(9.1204*L.^2); % regulates increase of contractility with Ls
Sarc.CDot= FL.*ft1 - C.*tanhx/tD; % 1st order rise and decay of [Ca++]
SfIso    = (C     .* L) * (1.51*SfAct) ;
SfRest   = L*(1.51*CRest*SfAct) ;

k1= 10; k2= 0.01; kk3= 2*(LsRef/dLsPas);
LfP   = exp(Ef)*(LsRef/Ls0Pas);
yEcm  = LfP.^k1;
SfEcm = 0.0349*(yEcm-1)*SfPas;% Ls=Ls0Pas, zero stress
y     = exp(log(LfP)*kk3);
SfTit = (y-1)*(k2*SfAct); % titin is softer than ecm, proportional with SfAct
SfPasT= SfTit + SfEcm;
DSfPasDEf= y*(k2*SfAct*kk3) + yEcm*(0.0349*k1*SfPas);

%=== Stress/Ls and stiffness, collected for output/problem Buckle!!+++++
Sarc.SfEcm  = SfEcm;% passive stress ECM
Sarc.SfPasT = SfPasT;% total passive stress +++++++++++ to be removed
LNormSe     = (Ls-Lsi)/LenSeriesElement;
Sarc.LsiDot = (LNormSe-1)*vMax;
Sarc.Sf     = SfPasT + (SfIso+SfRest).*LNormSe - SfRest;
Sarc.DSfDEf = DSfPasDEf+(SfIso.*Ls)/LenSeriesElement; % estimate of sarcomere stiffness

P.Patch=Sarc;
end

function tc=Tc(t,ttrig) % 'matrix proof'
% selects last trigger moment
% nt=number of samples in time
% nc=number of stored trigger moments
% nr=number of sarcomeres (or patches)
% selcts for each [nt,nr] event best trigger moment out of nc
% tc[nt,nr]= time-trigger moment
nt=size(t,1); %time length
[nc,nr]=size(ttrig); %number of triggers, number of sarc's
ttrigT=ttrig';
A=double(repmat(t,[1,nr*nc])>repmat(ttrigT(:)',[nt,1]));
DBnr=repmat((nc+1)*(0:nr-1),[nt,1]); %index shift per column nr
B=reshape( sum(reshape(A,[nt*nr,nc]),2), [nt,nr])+DBnr; %index trigger time matrix
b=[zeros(1,nr)-10;ttrig]; %trigger time matrix
tc=repmat(t,[1,nr])-b(B+1); %trigger time subtracted
end

function X=SparseDiag(x);
X=sparse(diag(x(:)));
end

