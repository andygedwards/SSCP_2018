function CircDisplayP

global P

SVarDot(0,P.SVar',[]);

% Variables as f(t)

% scaling of graphics
q0          = P.General.q0;
p0          = P.General.p0;
tCycle      = P.General.tCycle;

%=== scaling reference values
qSc= Rnd(q0);
VSc= Rnd(q0*tCycle/10);
pSc= Rnd(0.1*p0);

t= P.t-P.t(1);
OFFSET= -20;

%==== Hemodynamics
figure(1);
p1=GetFt('Node','p',{'SyArt','Lv','PuArt','Rv'})/pSc;
V1=GetFt('Cavity','V',{'La','Lv','Ra','Rv'})/VSc;
q1=GetFt('Valve','q',{'LvSyArt','LaLv','PuVenLa',...
    'RvPuArt','RaRv','SyVenRa','LaRa','LvRv','SyArtPuArt'})/qSc;
p1(:,[3,4])=p1(:,[3,4])+OFFSET;% pressures
V1(:,[3,4])=V1(:,[3,4])+OFFSET;% volumes
q1(:,[4:6])=q1(:,[4:6])+OFFSET;% flows
subplot(2,4,[3,4,7,8]); plot(t,[p1,V1,q1]);
title(['Units: ',num2str(qSc*1e6),' ml/s; ',num2str(pSc/1e3),' kPa; ',...
    num2str(VSc*1e6),' ml']);

% p-V loops
Calp= 0.001; CalV= 1e6;
VT= CalV*GetFt('Cavity','V',{'La','Ra','Lv','Rv'});
pT= Calp*GetFt('Cavity','p',{'La','Ra','Lv','Rv'});
subplot(2,4,1); plot(VT,pT);
title('p(V)')

%Sf-Ef loops
EfT = GetFt('Patch','Ef','All');
SfT = GetFt('Patch','Sf','All')/1e3;
subplot(2,4,2); plot(EfT,SfT);
title('stress[kPa](strain)')

% clipped low venous, atrial and ventricular pressures
A= [GetFt('Cavity','p',{'La','Ra','Lv','Rv'}),...
    GetFt('Bag'   ,'p',{'Peri'             }),...
    GetFt('Node'  ,'p',{'PuVen','SyVen'    })]/pSc;
Maxy=1.05*max(max(A(:,[1,2])));
Miny=min(A(:));
subplot(2,4,5:6); plot(t,A);
axis([t(1),t(end),Miny,Maxy]);
title('p(t)')

figure(1);
end

%========== AUXILARY FUNCTIONS =============

function X=Rnd(x)
X1= 10^round(log(x)/log(10));
X2=  2^round(log(x/X1)/log(2));
X=X1*X2;
end

