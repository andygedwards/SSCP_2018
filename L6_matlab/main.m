% MATLAB Implementation of the O'Hara-Rudy dynamic (ORd) model for the
% undiseased human ventricular action potential and calcium transient
%
% The ORd model is described in the article "Simulation of the Undiseased
% Human Cardiac Ventricular Action Potential: Model Formulation and
% Experimental Validation"
% by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
%
% The article and supplemental materails are freely available in the
% Open Access jounal PLoS Computational Biology
% Link to Article:
% http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
% 
% Email: tom.ohara@gmail.com / rudy@wustl.edu
% Web: http://rudylab.wustl.edu
% 
% The ORd model is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. The ORd model is distributed in the hope that
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details (www.gnu.org/licenses)

%close all
%clear all

%initial conditions for state variables
v=-87;
nai=7;
nass=nai;
ki=145;
kss=ki;
cai=1.0e-4;
cass=cai;
cansr=1.2;
cajsr=cansr;
m=0;
hf=1;
hs=1;
j=1;
hsp=1;
jp=1;
mL=0;
hL=1;
hLp=1;
a=0;
iF=1;
iS=1;
ap=0;
iFp=1;
iSp=1;
d=0;
ff=1;
fs=1;
fcaf=1;
fcas=1;
jca=1;
nca=0;
ffp=1;
fcafp=1;
xrf=0;
xrs=0;
xs1=0;
xs2=0;
xk1=1;
Jrelnp=0;
Jrelp=0;
CaMKt=0;
%X0 is the vector for initial sconditions for state variables
X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt]';

CL=500;%pacing cycle length in ms
beats=1;%number of beats in the simulation

options=[];%options for ode solver


for n=[1:beats]
    [time X]=ode15s(@model,[0 CL],X0,options,1);
    X0=X(size(X,1),:);
    n %output beat number to the screen to monitor runtime progress
end

%rename values in the state variables vector
v=X(:,1);
nai=X(:,2);
nass=X(:,3);
ki=X(:,4);
kss=X(:,5);
cai=X(:,6);
cass=X(:,7);
cansr=X(:,8);
cajsr=X(:,9);
m=X(:,10);
hf=X(:,11);
hs=X(:,12);
j=X(:,13);
hsp=X(:,14);
jp=X(:,15);
mL=X(:,16);
hL=X(:,17);
hLp=X(:,18);
a=X(:,19);
iF=X(:,20);
iS=X(:,21);
ap=X(:,22);
iFp=X(:,23);
iSp=X(:,24);
d=X(:,25);
ff=X(:,26);
fs=X(:,27);
fcaf=X(:,28);
fcas=X(:,29);
jca=X(:,30);
nca=X(:,31);
ffp=X(:,32);
fcafp=X(:,33);
xrf=X(:,34);
xrs=X(:,35);
xs1=X(:,36);
xs2=X(:,37);
xk1=X(:,38);
Jrelnp=X(:,39);
Jrelp=X(:,40);
CaMKt=X(:,41);

%calculate and name dependent variables for the final beat in the
%simulation (i.e. currents and fluxes)
for i=[1:size(X,1)];
    IsJs=model(time(i),X(i,:),0);
    INa(i)=IsJs(1);
    INaL(i)=IsJs(2);
    Ito(i)=IsJs(3);
    ICaL(i)=IsJs(4);
    IKr(i)=IsJs(5);
    IKs(i)=IsJs(6);
    IK1(i)=IsJs(7);
    INaCa_i(i)=IsJs(8);
    INaCa_ss(i)=IsJs(9);
    INaK(i)=IsJs(10);
    IKb(i)=IsJs(11);
    INab(i)=IsJs(12);
    ICab(i)=IsJs(13);
    IpCa(i)=IsJs(14);
    Jdiff(i)=IsJs(15);
    JdiffNa(i)=IsJs(16);
    JdiffK(i)=IsJs(17);
    Jup(i)=IsJs(18);
    Jleak(i)=IsJs(19);
    Jtr(i)=IsJs(20);
    Jrel(i)=IsJs(21);
    CaMKa(i)=IsJs(22);
    Istim(i)=IsJs(23);
end

%create plots showing results for the final paced beat

figure
subplot(2,3,1),plot(time,INa),title('INa')
subplot(2,3,2),plot(time,INaL),title('INaL')
subplot(2,3,3),plot(time,INaK),title('INaK')
subplot(2,3,4),plot(time,INaCa_i,time,INaCa_ss),title('INaCa_i,INaCa_ss')
subplot(2,3,5),plot(time,JdiffNa),title('JdiffNa')
subplot(2,3,6),plot(time,nai,time,nass),title('nai,nass')

figure
subplot(2,3,1),plot(time,Ito),title('Ito')
subplot(2,3,2),plot(time,IKr),title('IKr')
subplot(2,3,3),plot(time,IKs),title('IKs')
subplot(2,3,4),plot(time,IK1),title('IK1')
subplot(2,3,5),plot(time,INaK),title('INaK')
subplot(2,3,6),plot(time,ki),title('ki')

figure
subplot(2,3,1),plot(time,ICaL),title('ICaL')
subplot(2,3,2),plot(time,INaCa_i,time,INaCa_ss),title('INaCa_i,INaCa_ss')
subplot(2,3,3),plot(time,cai),title('cai')
subplot(2,3,4),plot(time,cass),title('cass')
subplot(2,3,5),plot(time,cansr,time,cajsr),title('cansr,cajsr')
subplot(2,3,6),plot(time,Jrel),title('Jrel')

figure
subplot(2,2,1),plot(time,ICaL),title('ICaL')
subplot(2,2,2),plot(time,cass),title('cass')
subplot(2,2,3),plot(time,cai),title('cai')
subplot(2,2,4),plot(time,Jrel),title('Jrel')

figure
subplot(3,3,1),plot(time,INa),title('INa')
subplot(3,3,2),plot(time,INaL),title('INaL')
subplot(3,3,3),plot(time,Ito),title('Ito')
subplot(3,3,5),plot(time,ICaL),title('ICaL')
subplot(3,3,6),plot(time,IKr),title('IKr')
subplot(3,3,7),plot(time,IKs),title('IKs')
subplot(3,3,8),plot(time,INaCa_i,time,INaCa_ss),title('INaCai,INaCass')
subplot(3,3,9),plot(time,INaK),title('INaK')

figure
subplot(2,3,1),plot(time,INab),title('INab')
subplot(2,3,2),plot(time,INaL),title('INaL')
subplot(2,3,3),plot(time,IKb),title('IKb')
subplot(2,3,4),plot(time,IpCa),title('IpCa')
subplot(2,3,5),plot(time,ICab),title('ICab')

figure
subplot(2,2,1),plot(time,Jrel),title('Jrel')
subplot(2,2,2),plot(time,Jup),title('Jup')
subplot(2,2,3),plot(time,Jleak),title('Jleak')
subplot(2,2,4),plot(time,Jtr),title('Jtr')

figure
subplot(2,2,1),plot(time,JdiffNa),title('JdiffNa')
subplot(2,2,2),plot(time,JdiffK),title('JdiffK')
subplot(2,2,3),plot(time,Jdiff),title('Jdiff')

figure
subplot(2,3,1),plot(time,nai,time,nass),title('nai,nass')
subplot(2,3,2),plot(time,cai),title('cai')
subplot(2,3,3),plot(time,cass),title('cass')
subplot(2,3,4),plot(time,cansr,time,cajsr),title('cansr,cajsr')
subplot(2,3,5),plot(time,ki),title('ki')

figure
plot(time,v,time,Istim),title('v,Istim')