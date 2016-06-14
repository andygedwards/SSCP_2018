function [ ] = results_ss_human_atrial(outfile, filename)
% Used to calculate some important values from human atrial myocyte model
% simulation data
%
% Code written by Topi Korhonen, Department of Physiology, University of Oulu
% Modified by Jussi Koivum??ki.
%

i_V = 1;
i_lastend = i_V;

% INa
i_start = i_lastend + 1;
i_INam = i_start; i_INah1 = i_start + 1; i_INah2 = i_start + 2;
i_lastend = i_INah2;

% ICaL 
i_start = i_lastend + 1;
i_ICaLd = i_start; i_ICaLf1 = i_start + 1; i_ICaLf2 = i_start + 2; i_ICaLfca = i_start + 3;
i_lastend = i_ICaLfca;

% It
i_start = i_lastend + 1;
i_Itr = i_start; i_Its = i_start + 1;
i_lastend = i_Its;

% Isus (Ikur)
i_start = i_lastend + 1;
i_Isusr = i_start; i_Isuss = i_start + 1;
i_lastend = i_Isuss;

% IKs
i_start = i_lastend + 1;
i_IKsn = i_start;
i_lastend = i_IKsn;

% IKr
i_start = i_lastend + 1;
i_IKrpa = i_start;
i_lastend = i_IKrpa;

% If
i_start = i_lastend + 1;
i_Ify = i_start;
i_lastend = i_Ify;

% gKCa
i_start = i_lastend + 1;
i_gKCa = i_start;
i_lastend = i_gKCa;

% RyR
i_start = i_lastend + 1;
i_RyRoss = i_start; i_RyRcss = i_start + 1; i_RyRass = i_start + 2;
i_RyRo1 = i_start + 3; i_RyRc1 = i_start + 4; i_RyRa1 =  i_start + 5;
i_RyRo2 = i_start + 6; i_RyRc2 = i_start + 7; i_RyRa2 =  i_start + 8;
i_RyRo3 = i_start + 9; i_RyRc3 = i_start + 10; i_RyRa3 =  i_start + 11;
i_lastend = i_RyRa3;

% SERCA
i_start = i_lastend + 1;
i_SERCACa1 = i_start; % first compartment from center
i_SERCACa2 = i_start + 1; % 2nd compartment from center
i_SERCACa3 = i_start + 2; % 3rd compartment from center
i_SERCACass = i_start + 3; % ss-SERCA
i_lastend = i_SERCACass;

% Nai ja Ki
i_start = i_lastend + 1;
i_Nass = i_start; i_Nai = i_start + 1;
i_Ki = i_start + 2;
i_lastend = i_Ki;

% Cai and CaSR 
i_start = i_lastend + 1;
i_Cass = i_start;
i_Cacenter = i_start + 1;

% cAF dilation
cAF_lcell = 1.0; % nSR
%cAF_lcell = 1.10; % cAF

% Physical & environmental constants
F = 96487;
R = 8314;
T = 306.15;

Nao = 130;
%Nai = y(i_Nai);

Cao = 1.8;

Ko = 5.4;
%Ki = y(i_Ki);

Cm = 0.05; %nF

% Cell dilation in cAF
Ddcell = (cAF_lcell - 1)*(20/10) + 1;
Dvcell = cAF_lcell*Ddcell^2;

% Geometry
Vss = 4.99232e-5 * Dvcell; %nL
rjunct = 6.5 * Ddcell; % mum
lcell = 122.051 * cAF_lcell; % mum

% Ca diffusion grid
dx = 1.625 * Ddcell; % mum
rstart = 0 + 0.5*dx;
rend = rjunct - 0.5*dx;
j = round(rstart/dx):1:round(rjunct/dx); % Spatial index of Cai diffusion
j = j';

Aj_nj = pi*rjunct*2*lcell*0.5; % Area between junct and nonjunct
xj_nj = 0.02/2 * Ddcell + dx/2; % diffusion distance from center to center of junct to first njunct
xj_nj_Nai = 0.02/2 * Ddcell + 2*dx; % diffusion distance from center of junct to center of njunct (between 2nd and 3rd njunct)

% Diffusion compartment volumes
Vnonjunct = zeros(length(j),1);
Vnonjunct = (pi.*(j.*dx).^2.*lcell-pi.*((j-1).*dx).^2.*lcell).*1e-6.*0.5; %nL

Vcytosol = sum(Vnonjunct) + Vss;

VSR = 0.05.*Vnonjunct./2*0.9 / Dvcell;

% SR Ca buffers
CSQN = 6.7; %6.7
KdCSQN = 0.8;

% Data import & data analysis
%style = '';

%for n = 1:(1+29+4+ 4*2 + 17) % functions, differential variables and time
%    style = [style ' %f'];
%end

%style = [style '\n'];

result = zeros(48, 1);
%dydata = zeros(29 +5 +4*2,1);

    data = load(filename);
    time = data(:,1);
%    dydata = data(:, 2:(1 + 1 + 14 + 12 + 4 + 3 + 5 + 4)); % time, V, ion channels, RyR, SERCA, Na&K, Cai and CaSR
%    dydata = data(:, 2:(1 + 1 + 14 + 12 + 4 + 3 + 5 + 4 + 1)); % IKCa
    dydata = data(:, 2:(1 + 28+4+3+1+2*4)); % time, differential variables

%    currentdata = data(1:end, (end - 24):(end)); % Ca fluxes
%    currentdata = data(1:end, (end - 25):(end)); % IKCa
    currentdata = data(1:end, (end - 26):(end)); % excluding Istim
    INa   = currentdata(:,1); 
    ICaL  = currentdata(:,2);
    It    = currentdata(:,3); 
    Isus  = currentdata(:,4);
    IKr   = currentdata(:,5);
    IKs   = currentdata(:,6);
    IK1   = currentdata(:,7);
    IKACh = currentdata(:,8);
    IKCa  = currentdata(:,9);
    If    = currentdata(:,10);
    INab  = currentdata(:,11);
    ICab  = currentdata(:,12);
    ICaP  = currentdata(:,13);
    INaK  = currentdata(:,14);
    INaCa = currentdata(:,15);

    Jrelss = currentdata(:,16); 
    Jrel1  = currentdata(:,17);
    Jrel2  = currentdata(:,18);
    Jrel3  = currentdata(:,19);
    J_bulkSERCA1 = currentdata(:,20);
    J_bulkSERCA2 = currentdata(:,21);
    J_bulkSERCA3 = currentdata(:,22);
    J_bulkSERCAss = currentdata(:,23); 
    JSRCaleak1 = currentdata(:,24);
    JSRCaleak2 = currentdata(:,25); 
    JSRCaleak3 = currentdata(:,26);
    JSRCaleakss = currentdata(:,27);  

    % Action potential.
    V = dydata(:, i_V);
    dVdt_max = max(diff(V)./diff(time));
    V_max = max(V);
    V_min = min(V);
    V_ampl = V_max - V_min;
    index_AP_huippu = find(V == V_max, 1);
    t_AP_huippu = time(find(V == max(V), 1));

    %t_shifted = time(index_AP_huippu:end);
    %V = V(index_AP_huippu:end);
    
    t_20_start = time(find(V >= (V_max - 0.20*V_ampl), 1, 'first'));
    t_20_end = time(find(V >= (V_max - 0.20*V_ampl), 1, 'last'));
    if numel(t_20_end) == 0
        APD_20 = 0;
    else
        APD_20 = t_20_end - t_20_start;
    end
    
    t_50_start = time(find(V >= (V_max - 0.50*V_ampl), 1, 'first'));
    t_50_end = time(find(V >= (V_max - 0.50*V_ampl), 1, 'last'));
    if numel(t_50_end) == 0
        APD_50 = 0;
    else
        APD_50 = t_50_end - t_50_start;
    end
    
    t_75_start = time(find(V >= (V_max - 0.75*V_ampl), 1, 'first'));
    t_75_end = time(find(V >= (V_max - 0.75*V_ampl), 1, 'last'));
    if numel(t_75_end) == 0
        APD_75 = 0;
    else
        APD_75 = t_75_end - t_75_start;
    end
        
    t_90_start = time(find(V >= (V_max - 0.90*V_ampl), 1, 'first'));
    t_90_end = time(find(V >= (V_max - 0.90*V_ampl), 1, 'last'));
    if numel(t_90_end) == 0
        APD_90 = 0;
    else
        APD_90 = t_90_end - t_90_start;
    end
        
    index_APD_20 = find(V >= (V_max - 0.20*V_ampl), 1, 'last');
    index_APD_30 = find(V >= (V_max - 0.30*V_ampl), 1, 'last');
    Plt_20 = mean(V(index_APD_20:index_APD_30));
        
    % Calcium transient.
    index_Ca_start = index_AP_huippu;
    cadata = dydata(index_Ca_start:end, i_Cacenter:i_Cacenter+3);
    cassdata = dydata(index_Ca_start:end, i_Cass);
    camean = (cadata(:, 1).*1.625 + cadata(:, 2).*1.625 + cadata(:, 3).*1.625 + cadata(:, 4).*1.625 + cassdata.*0.02)./6.52;
    
    Ca_cntr_peak = max(cadata(:, 1));
    camean_min = min(camean);
    camean_max = max(camean);
    Cai_amp = camean_max - camean_min;
    Cass_amp = max(cassdata) - min(cassdata);
    Cass_peak = max(cassdata);
    Ca_cntr_vs_ss = (max(cadata(:, 1)) - min(cadata(:, 1))) / (max(cassdata) - min(cassdata));

    % CaT time-to-peak, duration and decay
    index_Ca_huippu = find(camean == camean_max, 1);
    %t_Ca_huippu = time(camean == camean_max, 1);
    t_Ca_huippu = time(index_Ca_huippu);
    
    index_Ca_half = find(camean(index_Ca_huippu:end) <= (camean_max - 0.5*Cai_amp), 1);
    t_Ca_half = time(index_Ca_huippu + index_Ca_half);
    if numel(t_Ca_half) == 0
        t_Ca_50 = 0;
    else
        t_Ca_50 = t_Ca_half(1) - t_Ca_huippu;
    end
    
    index_Ca_eighty = find(camean(index_Ca_huippu:end) <= (camean_max - 0.8*Cai_amp), 1);
    t_Ca_eighty = time(index_Ca_huippu + index_Ca_eighty);
    if numel(t_Ca_eighty) == 0
        t_Ca_80 = 0;
    else
        t_Ca_80 = t_Ca_eighty(1) - t_Ca_huippu;
    end
    Ca_decay = 0;

    % Calcium fluxes.
    JCaL_int = -trapz(time, ICaL) / (2 * Vcytosol * F );
    Jrel_int = (trapz(time, Jrel1) + trapz(time, Jrel2) + trapz(time, Jrel3) + trapz(time, Jrelss) )/ Vcytosol;
    Jup_minus_Jleak_int = (trapz(time, (J_bulkSERCA1-JSRCaleak1)) + trapz(time, (J_bulkSERCA2-JSRCaleak2)) + trapz(time, (J_bulkSERCA3-JSRCaleak3)) + trapz(time, (J_bulkSERCAss-JSRCaleakss))) / Vcytosol;
    INaCa_abs = abs(INaCa);
    JNaCa_rev_int = trapz(time, (INaCa + INaCa_abs)/2)  / (Vcytosol * F);
    JNaCa_fwd_int = trapz(time, (INaCa - INaCa_abs)/2) / (Vcytosol * F);
    JpCa_int = trapz(time, ICaP)  / (2 * Vcytosol * F );
    JCab_int = -trapz(time, ICab)  / (2 * Vcytosol * F );
       
    Cain = JCaL_int + Jrel_int + JCab_int + JNaCa_rev_int;
    Caout = Jup_minus_Jleak_int + JpCa_int - JNaCa_fwd_int;

    %Calcium transport.
    SRATPase_fraction = Jup_minus_Jleak_int / Caout;
    NCX_fraction = -JNaCa_fwd_int / Caout;
    slow_fraction = JpCa_int / Caout;

    % SRCa
    freeCa = dydata(1:end,i_Cacenter+4:i_Cacenter+3+4);

    SRCa1 = (freeCa(1,1) + CSQN./(1+KdCSQN./freeCa(1,1))).*VSR(1)./Vcytosol;
    SRCa2 = (freeCa(1,2) + CSQN./(1+KdCSQN./freeCa(1,2))).*VSR(2)./Vcytosol;
    SRCa3 = (freeCa(1,3) + CSQN./(1+KdCSQN./freeCa(1,3))).*VSR(3)./Vcytosol;
    SRCa4 = (freeCa(1,4) + CSQN./(1+KdCSQN./freeCa(1,4))).*VSR(4)./Vcytosol;

    SRCa = SRCa1 + SRCa2 + SRCa3 + SRCa4;


    % Results.
    result(1) = V_min;
    result(2) = V_max;
    result(3) = V_ampl;
    result(4) = dVdt_max;
    result(5) = Plt_20;
    result(6) = APD_20;
    result(7) = APD_50;
    result(8) = APD_90;

    result(10) = camean_min;
    result(11) = camean_max;
    result(12) = Cai_amp;
    result(13) = Cass_amp;
    result(14) = Ca_cntr_peak;
    result(15) = Cass_peak;
    result(16) = Ca_cntr_vs_ss;
    result(17) = t_Ca_50;
    result(18) = t_Ca_80; 
    result(19) = t_Ca_huippu;
    
    result(21) = JCaL_int;
    result(22) = max(-ICaL);
    result(23) = Jrel_int;
    result(24) = Jup_minus_Jleak_int;    
    result(25) = JNaCa_rev_int;
    result(26) = max(INaCa);
    result(27) = JNaCa_fwd_int;
    result(28) = min(INaCa);
    
    result(30) = SRCa;
    result(31) = JCaL_int/(JCaL_int + Jrel_int);
    result(32) = SRATPase_fraction;
    result(33) = NCX_fraction;
    result(34) = slow_fraction;
    result(35) = Cain;
    result(36) = Caout;
    
    result(38) = max(INaK);
    result(39) = mean(dydata(:, i_Nass));
    result(40) = mean(dydata(:, i_Nai));
    result(41) = mean(dydata(:, i_Ki));

    result(43) = max(It); 
    result(44) = max(Isus); 
    result(45) = max(IKr); 
    result(46) = max(IKs); 
    result(47) = max(IK1); 
    result(48) = max(If);
    
    cadata = dydata(1:end,i_Cacenter:i_Cacenter+3);
    caSRdata = dydata(1:end,i_Cacenter+4:i_Cacenter+7);
    cassdata = dydata(1:end,i_Cass);
    x = 6.52 - [0.8125:1.625:5.6875 6.51]; % Reverse distance
    ca = [cadata cassdata]*1000;
    figure;
    [ch, ch] = contourf(time, x*cAF_lcell, ca', 100);
    colormap(jet);
    set(ch,'edgecolor','none'); 
    c = colorbar;
    caxis([0 1.2]);
    set(gca, 'FontName', 'Arial', 'FontSize', 14);
    xlabel('Time (s)');
    ylabel('Distance from sarcolemma');    

%% Save results to mat-file
outputfile_mat = [outfile '.mat'];
save(outputfile_mat, 'result');

%% Write results to a dat-file.
outputfile_dat = [outfile '.dat'];
style = '';
style = [style ' %16E'];
style = [style '\n'];

fid = fopen(outputfile_dat, 'w');

fprintf(fid, 'Variables \t');

fprintf(fid, '\nV_diast \t');
fprintf(fid, style, result(1,:));
fprintf(fid, 'V_syst \t');
fprintf(fid, style, result(2,:));
fprintf(fid, 'V_ampl \t');
fprintf(fid, style, result(3,:));
fprintf(fid, 'dVdt_max \t');
fprintf(fid, style, result(4,:));
fprintf(fid, 'Plt_20 \t');
fprintf(fid, style, result(5,:));
fprintf(fid, 'APD_20 \t');
fprintf(fid, style, result(6,:));
fprintf(fid, 'APD_50 \t');
fprintf(fid, style, result(7,:));
fprintf(fid, 'APD_90 \t');
fprintf(fid, style, result(8,:));

fprintf(fid, '-------- \t');
fprintf(fid, style, result(9,:));
fprintf(fid, 'Cai_dias \t');
fprintf(fid, style, result(10,:));
fprintf(fid, 'Cai_sys \t');
fprintf(fid, style, result(11,:));
fprintf(fid, 'Cai_amp \t');
fprintf(fid, style, result(12,:));
fprintf(fid, 'Cass_amp \t');
fprintf(fid, style, result(13,:));
fprintf(fid, 'Ca_cntr_peak \t');
fprintf(fid, style, result(14,:));
fprintf(fid, 'Cass_peak \t');
fprintf(fid, style, result(15,:));
fprintf(fid, 'Ca_cntr_vs_ss \t');
fprintf(fid, style, result(16,:));
fprintf(fid, 'CaTD_50 \t');
fprintf(fid, style, result(17,:));
fprintf(fid, 'CaTD_80 \t');
fprintf(fid, style, result(18,:));
fprintf(fid, 'Cai_TPT \t');
fprintf(fid, style, result(19,:));

fprintf(fid, '-------- \t');
fprintf(fid, style, result(20,:));
fprintf(fid, 'Jcal_int \t');
fprintf(fid, style, result(21,:));
fprintf(fid, 'Ical_peak \t');
fprintf(fid, style, result(22,:));
fprintf(fid, 'Jrel_int \t');
fprintf(fid, style, result(23,:));
fprintf(fid, 'Jup_-_leak_int \t');
fprintf(fid, style, result(24,:));
fprintf(fid, 'Jncx_rev_int \t');
fprintf(fid, style, result(25,:));
fprintf(fid, 'Jncx_rev_mx \t');
fprintf(fid, style, result(26,:));
fprintf(fid, 'Jncx_fwd_int \t');
fprintf(fid, style, result(27,:));
fprintf(fid, 'Jncx_fwd_mx \t');
fprintf(fid, style, result(28,:));

fprintf(fid, '-------- \t');
fprintf(fid, style, result(29,:));
fprintf(fid, 'SRCa_cont \t');
fprintf(fid, style, result(30,:));
fprintf(fid, 'JCaL/L_+_RyR \t');
fprintf(fid, style, result(31,:));
fprintf(fid, 'SRup_frac \t');
fprintf(fid, style, result(32,:));
fprintf(fid, 'NCX_frac \t');
fprintf(fid, style, result(33,:));
fprintf(fid, 'slow_frac \t');
fprintf(fid, style, result(34,:));
fprintf(fid, 'Ca_in \t');
fprintf(fid, style, result(35,:));
fprintf(fid, 'Ca_out \t');
fprintf(fid, style, result(36,:));

fprintf(fid, '-------- \t');
fprintf(fid, style, result(37,:));
fprintf(fid, 'Inka_peak \t');
fprintf(fid, style, result(38,:));
fprintf(fid, 'Nass \t\t');
fprintf(fid, style, result(39,:));
fprintf(fid, 'Nai \t\t');
fprintf(fid, style, result(40,:));
fprintf(fid, 'Ki \t\t');
fprintf(fid, style, result(41,:));

fprintf(fid, '-------- \t');
fprintf(fid, style, result(42,:));
fprintf(fid, 'It_max \t');
fprintf(fid, style, result(43,:));
fprintf(fid, 'Isus_max \t');
fprintf(fid, style, result(44,:));
fprintf(fid, 'IKr_mx \t');
fprintf(fid, style, result(45,:));
fprintf(fid, 'IKs_max \t');
fprintf(fid, style, result(46,:));
fprintf(fid, 'IK1_max \t');
fprintf(fid, style, result(47,:));
fprintf(fid, 'If_max \t');
fprintf(fid, style, result(48,:));

fclose(fid);

end