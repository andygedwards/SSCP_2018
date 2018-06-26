function TimingP

%SHOULD BE IMPROVED!!!!

global P

TauDisp=0.0; % dispersion of activation per wall, if relevant

% Reading timing parameters
G        = P.General    ;
tCycleRef= G.tCycleRest ; % Reference tCycle for body size scaling
tCycle   = G.tCycle     ; % current cycle time
TimeFac  = G.TimeFac    ; % scaling of contraction time
TauAv    = G.TauAv      ; % AV-delay, controlled in Adapt0

% TimeAct [s] duration of contraction for Ls==LsStress0Act (s)
ta= 0.15*(tCycle/0.85)*TimeFac; % atrial activation duration
tv= (0.10*tCycleRef +0.40*tCycle)*TimeFac; % ventricular " "

iWLv= P.TriSeg.iWall;
iWSv= iWLv+1; iWRv= iWLv+2;
iWLa= P.Chamber.iWall(1);
iWRa= P.Chamber.iWall(2);

iWVec= [iWLv,iWSv,iWRv,iWLa,iWRa];
tAVec= [tv,tv,tv,ta,ta];

% Delay times electrical activation
tRa2La= 0.0284 * (tCycleRef/0.85) * TimeFac;
tRa2Rv= TauAv+P.General.dTauAv; % AV-Delay
tRv2Lv= 0.0;       % pacing delay LBBB, negative-> prepacing
tRv2Sv= 0.33*tRv2Lv;  % pacing delay LBBB, negative-> prepacing
tRa2Ra= tCycle;      % cycle time
% tRa2Ra= tRa2Ra*(1+0.3*(rand-0.5)); % irregular HR


%TRIAL, CAN BE IMPROVED+++++++++++++++++++++++++++++++
% Time interval of simulation
tStart= P.t(end);
tEnd  = tStart+tCycle; % assume Ra as main trigger for activation

tRa=tStart; tLa=[]; tLv=[]; tSv=[]; tRv=[];
t=tStart;
while t<tEnd;
    tLa = [tLa; tRa(end) + tRa2La];
    tRv = [tRv; tRa(end) + tRa2Rv];
    tSv = [tSv; tRv(end) + tRv2Sv];
    tLv = [tLv; tRv(end) + tRv2Lv];
    tRa = [tRa; tRa(end) + tRa2Ra];
    t= tRa(end); % start new cycle
end
tRa=tRa(1:end-1); % tStart is begin of next cycle with Ra activation

iWVec= [iWLv,iWSv,iWRv,iWLa,iWRa]; % wall indices
tWVec= [tLv,tSv,tRv,tLa,tRa]; % trigger times of walls

%=== (multi-)patch activation, not tested yet
for i=1:5 %all walls
    iW= iWVec(i); tW=tWVec(i);
    nP= P.Wall.nPatch(iW);
    RgP= P.Wall.iPatch(iW)+(0:nP-1);
    P.Patch.TimeAct(RgP)= tAVec(i);
    tP=tW+TauDisp*((1:nP)-(nP+1)/2)/nP;
    Aux= [tP+P.Patch.dT(:,RgP);tP+P.Patch.dT(:,RgP)];% +++++++JOOST
    Aux(1,:)=Aux(1,:)-tCycle;
    P.Patch.ActivationDelay(:,RgP)= Aux;
end
end

