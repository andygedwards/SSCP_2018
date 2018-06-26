%% general data
t = P.t-min(P.t);
tcycle = P.General.tCycle;
ncycle = round(t(end)/tcycle);

icycle = round(P.General.tCycle/P.General.Dt)+1;

%% pre-allocation
LVEDV   = nan(1,ncycle);
LVESV   = nan(1,ncycle);
RVEDV   = nan(1,ncycle);
RVESV   = nan(1,ncycle);
mLAP    = nan(1,ncycle);
MAP     = nan(1,ncycle);
mRAP    = nan(1,ncycle);
mPAP    = nan(1,ncycle);
SfMax_R  = nan(1,ncycle);
SfMax_S  = nan(1,ncycle);
SfMax_L  = nan(1,ncycle);
WDen_R   = nan(1,ncycle);
WDen_S   = nan(1,ncycle);
WDen_L   = nan(1,ncycle);
LVSP     = nan(1,ncycle);
RVSP     = nan(1,ncycle);
ExtW_LV  = nan(1,ncycle);
ExtW_RV  = nan(1,ncycle);

SVtot_LV   = nan(1,ncycle);
SVfwd_LV   = nan(1,ncycle);
SVreg_LV   = nan(1,ncycle);

SVtot_RV   = nan(1,ncycle);
SVfwd_RV   = nan(1,ncycle);
SVreg_RV   = nan(1,ncycle);

%% signals
VLV_all = GetFt('Cavity','V','Lv');
VRV_all = GetFt('Cavity','V','Rv');
VLA_all = GetFt('Cavity','V','La');
VRA_all = GetFt('Cavity','V','Ra');

pLV_all = GetFt('Node','p','Lv');
pRV_all = GetFt('Node','p','Rv');
pLA_all = GetFt('Node','p','La');
pRA_all = GetFt('Node','p','Ra');
pAO_all = GetFt('Node','p','SyArt');
pPA_all = GetFt('Node','p','PuArt');

qAO_all = GetFt( 'Valve', 'q', 'LvSyArt' );
qMV_all = GetFt( 'Valve', 'q', 'LaLv' );
qPV_all = GetFt( 'Valve', 'q', 'RvPuArt' );
qTV_all = GetFt( 'Valve', 'q', 'RaRv' );

% EfL_all = GetFt('Patch','Ef','Lv1');
% EfS_all = GetFt('Patch','Ef','Sv1');
% EfR_all = GetFt('Patch','Ef','Rv1');
% SfL_all = GetFt('Patch','Sf','Lv1');
% SfS_all = GetFt('Patch','Sf','Sv1');
% SfR_all = GetFt('Patch','Sf','Rv1');

EfL_all = mean(GetFt('Patch','Ef','Lv'),2);
EfS_all = mean(GetFt('Patch','Ef','Sv'),2);
EfR_all = mean(GetFt('Patch','Ef','Rv'),2);
SfL_all = mean(GetFt('Patch','Sf','Lv'),2);
SfS_all = mean(GetFt('Patch','Sf','Sv'),2);
SfR_all = mean(GetFt('Patch','Sf','Rv'),2);

%% beat separation
for j = 1:ncycle
    
    range(j,1) = 1+(j-1)*icycle;
    range(j,2) = icycle+(j-1)*icycle;
    
    VLV(:,j) = VLV_all(range(j,1):range(j,2));
    VRV(:,j) = VRV_all(range(j,1):range(j,2));
    VLA(:,j) = VLA_all(range(j,1):range(j,2));
    VRA(:,j) = VRA_all(range(j,1):range(j,2));
    
    pLV(:,j) = pLV_all(range(j,1):range(j,2));
    pRV(:,j) = pRV_all(range(j,1):range(j,2));
    pLA(:,j) = pLA_all(range(j,1):range(j,2));
    pRA(:,j) = pRA_all(range(j,1):range(j,2));
    pAO(:,j) = pAO_all(range(j,1):range(j,2));
    pPA(:,j) = pPA_all(range(j,1):range(j,2));

    qAO(:,j) = qAO_all(range(j,1):range(j,2));
    qMV(:,j) = qMV_all(range(j,1):range(j,2));
    qPV(:,j) = qPV_all(range(j,1):range(j,2));
    qTV(:,j) = qTV_all(range(j,1):range(j,2));

    EfL(:,j) = EfL_all(range(j,1):range(j,2));
    EfS(:,j) = EfS_all(range(j,1):range(j,2));
    EfR(:,j) = EfR_all(range(j,1):range(j,2));
    
    SfL(:,j) = SfL_all(range(j,1):range(j,2));
    SfS(:,j) = SfS_all(range(j,1):range(j,2));
    SfR(:,j) = SfR_all(range(j,1):range(j,2));
    
end

%% beat-to-beat analysis
for i = 1:ncycle
    
    % Hemodynamic data
    LVEDV(i)    = max(VLV(:,i))*1e6;
    LVESV(i)    = min(VLV(:,i))*1e6;
    mLAP(i)     = mean(pLA(:,i))/133.322;
    MAP(i)      = mean(pAO(:,i))/133.322;
    LVSP(i)     = max(pLV(:,i))/133.322;
    RVEDV(i)    = max(VRV(:,i))*1e6;
    RVESV(i)    = min(VRV(:,i))*1e6;
    mRAP(i)     = mean(pRA(:,i))/133.322;
    mPAP(i)     = mean(pPA(:,i))/133.322;
    RVSP(i)     = max(pRV(:,i))/133.322;
    
    % Stroke volumes
    SVtot_LV(i) = LVEDV(i) - LVESV(i);
    SVfwd_LV(i) = trapz( P.General.Dt:P.General.Dt:P.General.Dt*(length(qAO(:,i)-1)), qAO(:,i) )*1e6;
    SVreg_LV(i) = SVtot_LV(i) - SVfwd_LV(i);
    SVtot_RV(i) = RVEDV(i) - RVESV(i);
    SVfwd_RV(i) = trapz( P.General.Dt:P.General.Dt:P.General.Dt*(length(qPV(:,i)-1)), qPV(:,i) )*1e6;
    SVreg_RV(i) = SVtot_RV(i) - SVfwd_RV(i);
    
    % Myocardial mechanics
    WDen_L(i)   = CalculateSegmentalWorkDen_inp(EfL(:,i),SfL(:,i));
    WDen_S(i)   = CalculateSegmentalWorkDen_inp(EfS(:,i),SfS(:,i));
    WDen_R(i)   = CalculateSegmentalWorkDen_inp(EfR(:,i),SfR(:,i));
    SfMax_L(i)  = max(SfL(:,i));
    SfMax_S(i)  = max(SfS(:,i));
    SfMax_R(i)  = max(SfR(:,i));

    % External pump work
    ExtW_LV(i)  = CalculateExternalWork_inp(VLV(:,i),pLV(:,i));
    ExtW_RV(i)  = CalculateExternalWork_inp(VRV(:,i),pRV(:,i));
    
end

