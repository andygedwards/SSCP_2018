
function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =44;
end
% There are a total of 15 entries in each of the rate and state variable arrays.
% There are a total of 48 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 5];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    %plot(VOI, STATES);
    % Plot only V_m
    plot(VOI, STATES(:,1));
    %xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES(1,:));
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (second)');
    LEGEND_STATES(:,1) = strpad('V in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component membrane (joule_per_kilomole_kelvin)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component membrane (coulomb_per_mole)');
    LEGEND_CONSTANTS(:,4) = strpad('C in component membrane (microF)');
    LEGEND_CONSTANTS(:,44) = strpad('RTONF in component membrane (millivolt)');
    LEGEND_ALGEBRAIC(:,26) = strpad('i_f in component hyperpolarising_activated_current (nanoA)');
    LEGEND_ALGEBRAIC(:,28) = strpad('i_K in component time_dependent_potassium_current (nanoA)');
    LEGEND_ALGEBRAIC(:,29) = strpad('i_K1 in component time_independent_potassium_current (nanoA)');
    LEGEND_ALGEBRAIC(:,30) = strpad('i_Na_b in component sodium_background_current (nanoA)');
    LEGEND_ALGEBRAIC(:,32) = strpad('i_Ca_b in component calcium_background_current (nanoA)');
    LEGEND_ALGEBRAIC(:,33) = strpad('i_p in component sodium_potassium_pump (nanoA)');
    LEGEND_ALGEBRAIC(:,34) = strpad('i_NaCa in component Na_Ca_exchanger (nanoA)');
    LEGEND_ALGEBRAIC(:,36) = strpad('i_Na in component fast_sodium_current (nanoA)');
    LEGEND_ALGEBRAIC(:,43) = strpad('i_si in component second_inward_current (nanoA)');
    LEGEND_ALGEBRAIC(:,23) = strpad('i_fNa in component hyperpolarising_activated_current (nanoA)');
    LEGEND_ALGEBRAIC(:,1) = strpad('E_Na in component hyperpolarising_activated_current (millivolt)');
    LEGEND_ALGEBRAIC(:,10) = strpad('E_K in component hyperpolarising_activated_current (millivolt)');
    LEGEND_ALGEBRAIC(:,25) = strpad('i_fK in component hyperpolarising_activated_current (nanoA)');
    LEGEND_CONSTANTS(:,5) = strpad('g_f_Na in component hyperpolarising_activated_current (microS)');
    LEGEND_CONSTANTS(:,6) = strpad('g_f_K in component hyperpolarising_activated_current (microS)');
    LEGEND_CONSTANTS(:,7) = strpad('Km_f in component hyperpolarising_activated_current (millimolar)');
    LEGEND_STATES(:,2) = strpad('Kc in component extracellular_potassium_concentration (millimolar)');
    LEGEND_STATES(:,3) = strpad('Ki in component intracellular_potassium_concentration (millimolar)');
    LEGEND_STATES(:,4) = strpad('Nai in component intracellular_sodium_concentration (millimolar)');
    LEGEND_CONSTANTS(:,8) = strpad('Nao in component extracellular_sodium_concentration (millimolar)');
    LEGEND_STATES(:,5) = strpad('y in component hyperpolarising_activated_current_y_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,11) = strpad('alpha_y in component hyperpolarising_activated_current_y_gate (per_second)');
    LEGEND_ALGEBRAIC(:,18) = strpad('beta_y in component hyperpolarising_activated_current_y_gate (per_second)');
    LEGEND_CONSTANTS(:,9) = strpad('delta_y in component hyperpolarising_activated_current_y_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,2) = strpad('E0_y in component hyperpolarising_activated_current_y_gate (millivolt)');
    LEGEND_CONSTANTS(:,10) = strpad('speed_y in component hyperpolarising_activated_current_y_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,27) = strpad('I_K in component time_dependent_potassium_current (nanoA)');
    LEGEND_CONSTANTS(:,11) = strpad('i_K_max in component time_dependent_potassium_current (nanoA)');
    LEGEND_STATES(:,6) = strpad('x in component time_dependent_potassium_current_x_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,12) = strpad('alpha_x in component time_dependent_potassium_current_x_gate (per_second)');
    LEGEND_ALGEBRAIC(:,19) = strpad('beta_x in component time_dependent_potassium_current_x_gate (per_second)');
    LEGEND_CONSTANTS(:,12) = strpad('delta_x in component time_dependent_potassium_current_x_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,3) = strpad('E0_x in component time_dependent_potassium_current_x_gate (millivolt)');
    LEGEND_CONSTANTS(:,13) = strpad('g_K1 in component time_independent_potassium_current (microS)');
    LEGEND_CONSTANTS(:,14) = strpad('Km_K1 in component time_independent_potassium_current (millimolar)');
    LEGEND_CONSTANTS(:,15) = strpad('g_Nab in component sodium_background_current (microS)');
    LEGEND_ALGEBRAIC(:,31) = strpad('E_Ca in component calcium_background_current (millivolt)');
    LEGEND_CONSTANTS(:,16) = strpad('g_Cab in component calcium_background_current (microS)');
    LEGEND_STATES(:,7) = strpad('Cai in component intracellular_calcium_concentration (millimolar)');
    LEGEND_CONSTANTS(:,17) = strpad('Cao in component extracellular_calcium_concentration (millimolar)');
    LEGEND_CONSTANTS(:,18) = strpad('I_p in component sodium_potassium_pump (nanoA)');
    LEGEND_CONSTANTS(:,19) = strpad('K_mK in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,20) = strpad('K_mNa in component sodium_potassium_pump (millimolar)');
    LEGEND_CONSTANTS(:,21) = strpad('n_NaCa in component Na_Ca_exchanger (dimensionless)');
    LEGEND_CONSTANTS(:,22) = strpad('K_NaCa in component Na_Ca_exchanger (nanoA)');
    LEGEND_CONSTANTS(:,23) = strpad('d_NaCa in component Na_Ca_exchanger (dimensionless)');
    LEGEND_CONSTANTS(:,24) = strpad('gamma in component Na_Ca_exchanger (dimensionless)');
    LEGEND_CONSTANTS(:,25) = strpad('g_Na in component fast_sodium_current (microS)');
    LEGEND_ALGEBRAIC(:,35) = strpad('E_mh in component fast_sodium_current (millivolt)');
    LEGEND_STATES(:,8) = strpad('m in component fast_sodium_current_m_gate (dimensionless)');
    LEGEND_STATES(:,9) = strpad('h in component fast_sodium_current_h_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,13) = strpad('alpha_m in component fast_sodium_current_m_gate (per_second)');
    LEGEND_ALGEBRAIC(:,20) = strpad('beta_m in component fast_sodium_current_m_gate (per_second)');
    LEGEND_CONSTANTS(:,26) = strpad('delta_m in component fast_sodium_current_m_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,4) = strpad('E0_m in component fast_sodium_current_m_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,5) = strpad('alpha_h in component fast_sodium_current_h_gate (per_second)');
    LEGEND_ALGEBRAIC(:,14) = strpad('beta_h in component fast_sodium_current_h_gate (per_second)');
    LEGEND_ALGEBRAIC(:,37) = strpad('i_siCa in component second_inward_current (nanoA)');
    LEGEND_ALGEBRAIC(:,38) = strpad('i_siK in component second_inward_current (nanoA)');
    LEGEND_ALGEBRAIC(:,40) = strpad('i_siNa in component second_inward_current (nanoA)');
    LEGEND_CONSTANTS(:,27) = strpad('P_si in component second_inward_current (nanoA_per_millimolar)');
    LEGEND_STATES(:,10) = strpad('d in component second_inward_current_d_gate (dimensionless)');
    LEGEND_STATES(:,11) = strpad('f in component second_inward_current_f_gate (dimensionless)');
    LEGEND_STATES(:,12) = strpad('f2 in component second_inward_current_f2_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,15) = strpad('alpha_d in component second_inward_current_d_gate (per_second)');
    LEGEND_ALGEBRAIC(:,21) = strpad('beta_d in component second_inward_current_d_gate (per_second)');
    LEGEND_CONSTANTS(:,28) = strpad('delta_d in component second_inward_current_d_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,6) = strpad('E0_d in component second_inward_current_d_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,16) = strpad('alpha_f in component second_inward_current_f_gate (per_second)');
    LEGEND_ALGEBRAIC(:,22) = strpad('beta_f in component second_inward_current_f_gate (per_second)');
    LEGEND_CONSTANTS(:,29) = strpad('delta_f in component second_inward_current_f_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,7) = strpad('E0_f in component second_inward_current_f_gate (millivolt)');
    LEGEND_CONSTANTS(:,30) = strpad('alpha_f2 in component second_inward_current_f2_gate (per_second)');
    LEGEND_ALGEBRAIC(:,8) = strpad('beta_f2 in component second_inward_current_f2_gate (per_second)');
    LEGEND_CONSTANTS(:,31) = strpad('K_mf2 in component second_inward_current_f2_gate (millimolar)');
    LEGEND_CONSTANTS(:,32) = strpad('radius in component intracellular_sodium_concentration (micrometre)');
    LEGEND_CONSTANTS(:,33) = strpad('length in component intracellular_sodium_concentration (micrometre)');
    LEGEND_CONSTANTS(:,34) = strpad('V_e_ratio in component intracellular_sodium_concentration (dimensionless)');
    LEGEND_CONSTANTS(:,45) = strpad('V_Cell in component intracellular_sodium_concentration (micrometre3)');
    LEGEND_CONSTANTS(:,46) = strpad('Vi in component intracellular_sodium_concentration (micrometre3)');
    LEGEND_CONSTANTS(:,47) = strpad('V_up in component intracellular_calcium_concentration (micrometre3)');
    LEGEND_CONSTANTS(:,48) = strpad('V_rel in component intracellular_calcium_concentration (micrometre3)');
    LEGEND_ALGEBRAIC(:,39) = strpad('i_up in component intracellular_calcium_concentration (nanoA)');
    LEGEND_ALGEBRAIC(:,41) = strpad('i_tr in component intracellular_calcium_concentration (nanoA)');
    LEGEND_ALGEBRAIC(:,44) = strpad('i_rel in component intracellular_calcium_concentration (nanoA)');
    LEGEND_STATES(:,13) = strpad('Ca_up in component intracellular_calcium_concentration (millimolar)');
    LEGEND_STATES(:,14) = strpad('Ca_rel in component intracellular_calcium_concentration (millimolar)');
    LEGEND_CONSTANTS(:,35) = strpad('Ca_up_max in component intracellular_calcium_concentration (millimolar)');
    LEGEND_CONSTANTS(:,36) = strpad('K_mCa in component intracellular_calcium_concentration (millimolar)');
    LEGEND_STATES(:,15) = strpad('p in component intracellular_calcium_concentration (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('alpha_p in component intracellular_calcium_concentration (per_second)');
    LEGEND_ALGEBRAIC(:,24) = strpad('beta_p in component intracellular_calcium_concentration (per_second)');
    LEGEND_ALGEBRAIC(:,9) = strpad('E0_p in component intracellular_calcium_concentration (millivolt)');
    LEGEND_CONSTANTS(:,37) = strpad('tau_up in component intracellular_calcium_concentration (second)');
    LEGEND_CONSTANTS(:,38) = strpad('tau_rep in component intracellular_calcium_concentration (second)');
    LEGEND_CONSTANTS(:,39) = strpad('tau_rel in component intracellular_calcium_concentration (second)');
    LEGEND_CONSTANTS(:,40) = strpad('rCa in component intracellular_calcium_concentration (dimensionless)');
    LEGEND_CONSTANTS(:,41) = strpad('V_e in component extracellular_potassium_concentration (micrometre3)');
    LEGEND_CONSTANTS(:,42) = strpad('Kb in component extracellular_potassium_concentration (millimolar)');
    LEGEND_ALGEBRAIC(:,42) = strpad('i_mK in component extracellular_potassium_concentration (nanoA)');
    LEGEND_CONSTANTS(:,43) = strpad('pf in component extracellular_potassium_concentration (per_second)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES(:,5) = strpad('d/dt y in component hyperpolarising_activated_current_y_gate (dimensionless)');
    LEGEND_RATES(:,6) = strpad('d/dt x in component time_dependent_potassium_current_x_gate (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt m in component fast_sodium_current_m_gate (dimensionless)');
    LEGEND_RATES(:,9) = strpad('d/dt h in component fast_sodium_current_h_gate (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt d in component second_inward_current_d_gate (dimensionless)');
    LEGEND_RATES(:,11) = strpad('d/dt f in component second_inward_current_f_gate (dimensionless)');
    LEGEND_RATES(:,12) = strpad('d/dt f2 in component second_inward_current_f2_gate (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt Nai in component intracellular_sodium_concentration (millimolar)');
    LEGEND_RATES(:,15) = strpad('d/dt p in component intracellular_calcium_concentration (dimensionless)');
    LEGEND_RATES(:,13) = strpad('d/dt Ca_up in component intracellular_calcium_concentration (millimolar)');
    LEGEND_RATES(:,14) = strpad('d/dt Ca_rel in component intracellular_calcium_concentration (millimolar)');
    LEGEND_RATES(:,7) = strpad('d/dt Cai in component intracellular_calcium_concentration (millimolar)');
    LEGEND_RATES(:,2) = strpad('d/dt Kc in component extracellular_potassium_concentration (millimolar)');
    LEGEND_RATES(:,3) = strpad('d/dt Ki in component intracellular_potassium_concentration (millimolar)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = -60;
    CONSTANTS(:,1) = 8314.472;
    CONSTANTS(:,2) = 310;
    CONSTANTS(:,3) = 96485.3415;
    CONSTANTS(:,4) = 0.006;
    CONSTANTS(:,5) = 6; %g_f_Na orig = 6
    CONSTANTS(:,6) = 6; %g_f_K orig = 6
    CONSTANTS(:,7) = 45;
    STATES(:,2) = 3;
    STATES(:,3) = 140;
    STATES(:,4) = 7.5;
    CONSTANTS(:,8) = 140;
    STATES(:,5) = 0.007;
    CONSTANTS(:,9) = 1e-5;
    CONSTANTS(:,10) = 2;
    CONSTANTS(:,11) = 20;
    STATES(:,6) = 0.54;
    CONSTANTS(:,12) = 0.0001;
    CONSTANTS(:,13) = 0.75;
    CONSTANTS(:,14) = 10;
    CONSTANTS(:,15) = 0.07;
    CONSTANTS(:,16) = 0.01;
    STATES(:,7) = 5.8e-5;
    CONSTANTS(:,17) = 2;
    CONSTANTS(:,18) = 50;
    CONSTANTS(:,19) = 1;
    CONSTANTS(:,20) = 40;
    CONSTANTS(:,21) = 3;
    CONSTANTS(:,22) = 0.002;
    CONSTANTS(:,23) = 0.0001;
    CONSTANTS(:,24) = 0.5;
    CONSTANTS(:,25) = 1.25; 
    STATES(:,8) = 0.076;
    STATES(:,9) = 0.015;
    CONSTANTS(:,26) = 1e-5;
    CONSTANTS(:,27) = 7.5;
    STATES(:,10) = 0.0011;
    STATES(:,11) = 0.785;
    STATES(:,12) = 0.785;
    CONSTANTS(:,28) = 0.0001;
    CONSTANTS(:,29) = 0.0001;
    CONSTANTS(:,30) = 10;
    CONSTANTS(:,31) = 0.0005;
    CONSTANTS(:,32) = 0.08;
    CONSTANTS(:,33) = 0.08;
    CONSTANTS(:,34) = 0.1;
    STATES(:,13) = 1.98;
    STATES(:,14) = 0.55;
    CONSTANTS(:,35) = 5;
    CONSTANTS(:,36) = 0.002;
    STATES(:,15) = 0.785;
    CONSTANTS(:,37) = 0.005;
    CONSTANTS(:,38) = 0.2;
    CONSTANTS(:,39) = 0.01;
    CONSTANTS(:,40) = 2;
    CONSTANTS(:,41) = 0.00016077;
    CONSTANTS(:,42) = 3;
    CONSTANTS(:,43) = 1;
    CONSTANTS(:,44) = ( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3);
    CONSTANTS(:,45) =  3.14159.*power(CONSTANTS(:,32), 2.00000).*CONSTANTS(:,33);
    CONSTANTS(:,46) =  CONSTANTS(:,45).*(1.00000 - CONSTANTS(:,34));
    CONSTANTS(:,47) =  CONSTANTS(:,46).*0.0500000;
    CONSTANTS(:,48) =  CONSTANTS(:,46).*0.0200000;
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    ALGEBRAIC(:,8) = ( STATES(:,7).*CONSTANTS(:,30))./CONSTANTS(:,31);
    RATES(:,12) = CONSTANTS(:,30) -  STATES(:,12).*(CONSTANTS(:,30)+ALGEBRAIC(:,8));
    ALGEBRAIC(:,5) =  20.0000.*exp(  - 0.125000.*(STATES(:,1)+75.0000));
    ALGEBRAIC(:,14) = 2000.00./( 320.000.*exp(  - 0.100000.*(STATES(:,1)+75.0000))+1.00000);
    RATES(:,9) =  ALGEBRAIC(:,5).*(1.00000 - STATES(:,9)) -  ALGEBRAIC(:,14).*STATES(:,9);
    ALGEBRAIC(:,2) = STATES(:,1)+52.0000;
    ALGEBRAIC(:,11) =  0.0500000.*exp(  - 0.0670000.*ALGEBRAIC(:,2));
    ALGEBRAIC(:,18) = piecewise({abs(ALGEBRAIC(:,2))<CONSTANTS(:,9), 2.50000 }, ALGEBRAIC(:,2)./(1.00000 -  1.00000.*exp(  - 0.200000.*ALGEBRAIC(:,2))));
    RATES(:,5) =  CONSTANTS(:,10).*( ALGEBRAIC(:,11).*(1.00000 - STATES(:,5)) -  ALGEBRAIC(:,18).*STATES(:,5));
    ALGEBRAIC(:,3) = STATES(:,1)+22.0000;
    ALGEBRAIC(:,12) = piecewise({abs(ALGEBRAIC(:,3))<CONSTANTS(:,12), 2.50000 }, ( 0.500000.*ALGEBRAIC(:,3))./(1.00000 - exp( - ALGEBRAIC(:,3)./5.00000)));
    ALGEBRAIC(:,19) = piecewise({abs(ALGEBRAIC(:,3))<CONSTANTS(:,12), 2.50000 }, ( 0.178000.*ALGEBRAIC(:,3))./(exp(ALGEBRAIC(:,3)./15.0000) - 1.00000));
    RATES(:,6) =  ALGEBRAIC(:,12).*(1.00000 - STATES(:,6)) -  ALGEBRAIC(:,19).*STATES(:,6);
    ALGEBRAIC(:,4) = STATES(:,1)+41.0000;
    ALGEBRAIC(:,13) = piecewise({abs(ALGEBRAIC(:,4))<CONSTANTS(:,26), 2000.00 }, ( 200.000.*ALGEBRAIC(:,4))./(1.00000 - exp(  - 0.100000.*ALGEBRAIC(:,4))));
    ALGEBRAIC(:,20) =  8000.00.*exp(  - 0.0560000.*(STATES(:,1)+66.0000));
    RATES(:,8) =  ALGEBRAIC(:,13).*(1.00000 - STATES(:,8)) -  ALGEBRAIC(:,20).*STATES(:,8);
    ALGEBRAIC(:,6) = (STATES(:,1)+24.0000) - 5.00000;
    ALGEBRAIC(:,15) = piecewise({abs(ALGEBRAIC(:,6))<CONSTANTS(:,28), 120.000 }, ( 30.0000.*ALGEBRAIC(:,6))./(1.00000 - exp((  - 1.00000.*ALGEBRAIC(:,6))./4.00000)));
    ALGEBRAIC(:,21) = piecewise({abs(ALGEBRAIC(:,6))<CONSTANTS(:,28), 120.000 }, ( 12.0000.*ALGEBRAIC(:,6))./(exp(ALGEBRAIC(:,6)./10.0000) - 1.00000));
    RATES(:,10) =  ALGEBRAIC(:,15).*(1.00000 - STATES(:,10)) -  ALGEBRAIC(:,21).*STATES(:,10);
    ALGEBRAIC(:,7) = STATES(:,1)+34.0000;
    ALGEBRAIC(:,16) = piecewise({abs(ALGEBRAIC(:,7))<CONSTANTS(:,29), 25.0000 }, ( 6.25000.*ALGEBRAIC(:,7))./(exp(ALGEBRAIC(:,7)./4.00000) - 1.00000));
    ALGEBRAIC(:,22) = 50.0000./(1.00000+exp((  - 1.00000.*(STATES(:,1)+34.0000))./4.00000));
    RATES(:,11) =  ALGEBRAIC(:,16).*(1.00000 - STATES(:,11)) -  ALGEBRAIC(:,22).*STATES(:,11);
    ALGEBRAIC(:,9) = (STATES(:,1)+34.0000) -  - 30.0000;
    ALGEBRAIC(:,17) = ( 0.625000.*ALGEBRAIC(:,9))./(exp(ALGEBRAIC(:,9)./4.00000) - 1.00000);
    ALGEBRAIC(:,24) = 5.00000./(1.00000+exp((  - 1.00000.*ALGEBRAIC(:,9))./4.00000));
    RATES(:,15) =  ALGEBRAIC(:,17).*(1.00000 - STATES(:,15)) -  ALGEBRAIC(:,24).*STATES(:,15);
    ALGEBRAIC(:,1) =  CONSTANTS(:,44).*log(CONSTANTS(:,8)./STATES(:,4));
    ALGEBRAIC(:,30) =  CONSTANTS(:,15).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,33) = ( (( CONSTANTS(:,18).*STATES(:,2))./(CONSTANTS(:,19)+STATES(:,2))).*STATES(:,4))./(CONSTANTS(:,20)+STATES(:,4));
    ALGEBRAIC(:,34) = ( CONSTANTS(:,22).*( exp(( CONSTANTS(:,24).*(CONSTANTS(:,21) - 2.00000).*STATES(:,1))./CONSTANTS(:,44)).*power(STATES(:,4), CONSTANTS(:,21)).*CONSTANTS(:,17) -  exp(( (CONSTANTS(:,24) - 1.00000).*(CONSTANTS(:,21) - 2.00000).*STATES(:,1))./CONSTANTS(:,44)).*power(CONSTANTS(:,8), CONSTANTS(:,21)).*STATES(:,7)))./( (1.00000+ CONSTANTS(:,23).*( STATES(:,7).*power(CONSTANTS(:,8), CONSTANTS(:,21))+ CONSTANTS(:,17).*power(STATES(:,4), CONSTANTS(:,21)))).*(1.00000+STATES(:,7)./0.00690000));
    ALGEBRAIC(:,35) =  CONSTANTS(:,44).*log((CONSTANTS(:,8)+ 0.120000.*STATES(:,2))./(STATES(:,4)+ 0.120000.*STATES(:,3)));
    ALGEBRAIC(:,36) =  CONSTANTS(:,25).*power(STATES(:,8), 3.00000).*STATES(:,9).*(STATES(:,1) - ALGEBRAIC(:,35));
    ALGEBRAIC(:,23) =  (( STATES(:,5).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,7))).*CONSTANTS(:,5).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,40) =  (( 0.0100000.*CONSTANTS(:,27).*(STATES(:,1) - 50.0000))./( CONSTANTS(:,44).*(1.00000 - exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))))).*( STATES(:,4).*exp(50.0000./CONSTANTS(:,44)) -  CONSTANTS(:,8).*exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))).*STATES(:,10).*STATES(:,11).*STATES(:,12);
    RATES(:,4) = (  - 1.00000.*(ALGEBRAIC(:,36)+ALGEBRAIC(:,30)+ALGEBRAIC(:,23)+ALGEBRAIC(:,40)+ ALGEBRAIC(:,33).*3.00000+( ALGEBRAIC(:,34).*CONSTANTS(:,21))./(CONSTANTS(:,21) - 2.00000)))./( 1.00000.*CONSTANTS(:,46).*CONSTANTS(:,3));
    ALGEBRAIC(:,39) =  (( 2.00000.*1.00000.*CONSTANTS(:,46).*CONSTANTS(:,3))./( 1.00000.*CONSTANTS(:,37).*CONSTANTS(:,35))).*STATES(:,7).*(CONSTANTS(:,35) - STATES(:,13));
    ALGEBRAIC(:,41) =  (( 2.00000.*1.00000.*CONSTANTS(:,48).*CONSTANTS(:,3))./( 1.00000.*CONSTANTS(:,38))).*STATES(:,15).*(STATES(:,13) - STATES(:,14));
    RATES(:,13) = ( 1.00000.*(ALGEBRAIC(:,39) - ALGEBRAIC(:,41)))./( 2.00000.*1.00000.*CONSTANTS(:,47).*CONSTANTS(:,3));
    ALGEBRAIC(:,27) = ( CONSTANTS(:,11).*(STATES(:,3) -  STATES(:,2).*exp( - STATES(:,1)./CONSTANTS(:,44))))./140.000;
    ALGEBRAIC(:,28) =  STATES(:,6).*ALGEBRAIC(:,27);
    ALGEBRAIC(:,10) =  CONSTANTS(:,44).*log(STATES(:,2)./STATES(:,3));
    ALGEBRAIC(:,29) = ( (( CONSTANTS(:,13).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,14))).*(STATES(:,1) - ALGEBRAIC(:,10)))./(1.00000+exp(( ((STATES(:,1)+10.0000) - ALGEBRAIC(:,10)).*2.00000)./CONSTANTS(:,44)));
    ALGEBRAIC(:,25) =  (( STATES(:,5).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,7))).*CONSTANTS(:,6).*(STATES(:,1) - ALGEBRAIC(:,10));
    ALGEBRAIC(:,38) =  (( 0.0100000.*CONSTANTS(:,27).*(STATES(:,1) - 50.0000))./( CONSTANTS(:,44).*(1.00000 - exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))))).*( STATES(:,3).*exp(50.0000./CONSTANTS(:,44)) -  STATES(:,2).*exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))).*STATES(:,10).*STATES(:,11).*STATES(:,12);
    ALGEBRAIC(:,42) = (ALGEBRAIC(:,29)+ALGEBRAIC(:,28)+ALGEBRAIC(:,25)+ALGEBRAIC(:,38)) -  2.00000.*ALGEBRAIC(:,33);
    RATES(:,2) =   - CONSTANTS(:,43).*(STATES(:,2) - CONSTANTS(:,42))+( 1.00000.*ALGEBRAIC(:,42))./( 1.00000.*CONSTANTS(:,41).*CONSTANTS(:,3));
    RATES(:,3) = (  - 1.00000.*ALGEBRAIC(:,42))./( 1.00000.*CONSTANTS(:,46).*CONSTANTS(:,3));
    ALGEBRAIC(:,26) = ALGEBRAIC(:,23)+ALGEBRAIC(:,25);
    ALGEBRAIC(:,31) =  0.500000.*CONSTANTS(:,44).*log(CONSTANTS(:,17)./STATES(:,7));
    ALGEBRAIC(:,32) =  CONSTANTS(:,16).*(STATES(:,1) - ALGEBRAIC(:,31));
    ALGEBRAIC(:,37) =  (( 4.00000.*CONSTANTS(:,27).*(STATES(:,1) - 50.0000))./( CONSTANTS(:,44).*(1.00000 - exp((  - 1.00000.*(STATES(:,1) - 50.0000).*2.00000)./CONSTANTS(:,44))))).*( STATES(:,7).*exp(100.000./CONSTANTS(:,44)) -  CONSTANTS(:,17).*exp((  - 2.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))).*STATES(:,10).*STATES(:,11).*STATES(:,12);
    ALGEBRAIC(:,43) = ALGEBRAIC(:,37)+ALGEBRAIC(:,38)+ALGEBRAIC(:,40);
    RATES(:,1) =  - (ALGEBRAIC(:,26)+ALGEBRAIC(:,28)+ALGEBRAIC(:,29)+ALGEBRAIC(:,30)+ALGEBRAIC(:,32)+ALGEBRAIC(:,33)+ALGEBRAIC(:,34)+ALGEBRAIC(:,36)+ALGEBRAIC(:,43))./CONSTANTS(:,4);
    ALGEBRAIC(:,44) = ( (( 2.00000.*1.00000.*CONSTANTS(:,48).*CONSTANTS(:,3))./( 1.00000.*CONSTANTS(:,39))).*STATES(:,14).*power(STATES(:,7), CONSTANTS(:,40)))./(power(STATES(:,7), CONSTANTS(:,40))+power(CONSTANTS(:,36), CONSTANTS(:,40)));
    RATES(:,14) = ( 1.00000.*(ALGEBRAIC(:,41) - ALGEBRAIC(:,44)))./( 2.00000.*1.00000.*CONSTANTS(:,48).*CONSTANTS(:,3));
    RATES(:,7) = (  - 1.00000.*((((ALGEBRAIC(:,37)+ALGEBRAIC(:,32)) - ( 2.00000.*ALGEBRAIC(:,34))./(CONSTANTS(:,21) - 2.00000)) - ALGEBRAIC(:,44))+ALGEBRAIC(:,39)))./( 2.00000.*1.00000.*CONSTANTS(:,46).*CONSTANTS(:,3));
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    ALGEBRAIC(:,8) = ( STATES(:,7).*CONSTANTS(:,30))./CONSTANTS(:,31);
    ALGEBRAIC(:,5) =  20.0000.*exp(  - 0.125000.*(STATES(:,1)+75.0000));
    ALGEBRAIC(:,14) = 2000.00./( 320.000.*exp(  - 0.100000.*(STATES(:,1)+75.0000))+1.00000);
    ALGEBRAIC(:,2) = STATES(:,1)+52.0000;
    ALGEBRAIC(:,11) =  0.0500000.*exp(  - 0.0670000.*ALGEBRAIC(:,2));
    ALGEBRAIC(:,18) = piecewise({abs(ALGEBRAIC(:,2))<CONSTANTS(:,9), 2.50000 }, ALGEBRAIC(:,2)./(1.00000 -  1.00000.*exp(  - 0.200000.*ALGEBRAIC(:,2))));
    ALGEBRAIC(:,3) = STATES(:,1)+22.0000;
    ALGEBRAIC(:,12) = piecewise({abs(ALGEBRAIC(:,3))<CONSTANTS(:,12), 2.50000 }, ( 0.500000.*ALGEBRAIC(:,3))./(1.00000 - exp( - ALGEBRAIC(:,3)./5.00000)));
    ALGEBRAIC(:,19) = piecewise({abs(ALGEBRAIC(:,3))<CONSTANTS(:,12), 2.50000 }, ( 0.178000.*ALGEBRAIC(:,3))./(exp(ALGEBRAIC(:,3)./15.0000) - 1.00000));
    ALGEBRAIC(:,4) = STATES(:,1)+41.0000;
    ALGEBRAIC(:,13) = piecewise({abs(ALGEBRAIC(:,4))<CONSTANTS(:,26), 2000.00 }, ( 200.000.*ALGEBRAIC(:,4))./(1.00000 - exp(  - 0.100000.*ALGEBRAIC(:,4))));
    ALGEBRAIC(:,20) =  8000.00.*exp(  - 0.0560000.*(STATES(:,1)+66.0000));
    ALGEBRAIC(:,6) = (STATES(:,1)+24.0000) - 5.00000;
    ALGEBRAIC(:,15) = piecewise({abs(ALGEBRAIC(:,6))<CONSTANTS(:,28), 120.000 }, ( 30.0000.*ALGEBRAIC(:,6))./(1.00000 - exp((  - 1.00000.*ALGEBRAIC(:,6))./4.00000)));
    ALGEBRAIC(:,21) = piecewise({abs(ALGEBRAIC(:,6))<CONSTANTS(:,28), 120.000 }, ( 12.0000.*ALGEBRAIC(:,6))./(exp(ALGEBRAIC(:,6)./10.0000) - 1.00000));
    ALGEBRAIC(:,7) = STATES(:,1)+34.0000;
    ALGEBRAIC(:,16) = piecewise({abs(ALGEBRAIC(:,7))<CONSTANTS(:,29), 25.0000 }, ( 6.25000.*ALGEBRAIC(:,7))./(exp(ALGEBRAIC(:,7)./4.00000) - 1.00000));
    ALGEBRAIC(:,22) = 50.0000./(1.00000+exp((  - 1.00000.*(STATES(:,1)+34.0000))./4.00000));
    ALGEBRAIC(:,9) = (STATES(:,1)+34.0000) -  - 30.0000;
    ALGEBRAIC(:,17) = ( 0.625000.*ALGEBRAIC(:,9))./(exp(ALGEBRAIC(:,9)./4.00000) - 1.00000);
    ALGEBRAIC(:,24) = 5.00000./(1.00000+exp((  - 1.00000.*ALGEBRAIC(:,9))./4.00000));
    ALGEBRAIC(:,1) =  CONSTANTS(:,44).*log(CONSTANTS(:,8)./STATES(:,4));
    ALGEBRAIC(:,30) =  CONSTANTS(:,15).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,33) = ( (( CONSTANTS(:,18).*STATES(:,2))./(CONSTANTS(:,19)+STATES(:,2))).*STATES(:,4))./(CONSTANTS(:,20)+STATES(:,4));
    ALGEBRAIC(:,34) = ( CONSTANTS(:,22).*( exp(( CONSTANTS(:,24).*(CONSTANTS(:,21) - 2.00000).*STATES(:,1))./CONSTANTS(:,44)).*power(STATES(:,4), CONSTANTS(:,21)).*CONSTANTS(:,17) -  exp(( (CONSTANTS(:,24) - 1.00000).*(CONSTANTS(:,21) - 2.00000).*STATES(:,1))./CONSTANTS(:,44)).*power(CONSTANTS(:,8), CONSTANTS(:,21)).*STATES(:,7)))./( (1.00000+ CONSTANTS(:,23).*( STATES(:,7).*power(CONSTANTS(:,8), CONSTANTS(:,21))+ CONSTANTS(:,17).*power(STATES(:,4), CONSTANTS(:,21)))).*(1.00000+STATES(:,7)./0.00690000));
    ALGEBRAIC(:,35) =  CONSTANTS(:,44).*log((CONSTANTS(:,8)+ 0.120000.*STATES(:,2))./(STATES(:,4)+ 0.120000.*STATES(:,3)));
    ALGEBRAIC(:,36) =  CONSTANTS(:,25).*power(STATES(:,8), 3.00000).*STATES(:,9).*(STATES(:,1) - ALGEBRAIC(:,35));
    ALGEBRAIC(:,23) =  (( STATES(:,5).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,7))).*CONSTANTS(:,5).*(STATES(:,1) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,40) =  (( 0.0100000.*CONSTANTS(:,27).*(STATES(:,1) - 50.0000))./( CONSTANTS(:,44).*(1.00000 - exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))))).*( STATES(:,4).*exp(50.0000./CONSTANTS(:,44)) -  CONSTANTS(:,8).*exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))).*STATES(:,10).*STATES(:,11).*STATES(:,12);
    ALGEBRAIC(:,39) =  (( 2.00000.*1.00000.*CONSTANTS(:,46).*CONSTANTS(:,3))./( 1.00000.*CONSTANTS(:,37).*CONSTANTS(:,35))).*STATES(:,7).*(CONSTANTS(:,35) - STATES(:,13));
    ALGEBRAIC(:,41) =  (( 2.00000.*1.00000.*CONSTANTS(:,48).*CONSTANTS(:,3))./( 1.00000.*CONSTANTS(:,38))).*STATES(:,15).*(STATES(:,13) - STATES(:,14));
    ALGEBRAIC(:,27) = ( CONSTANTS(:,11).*(STATES(:,3) -  STATES(:,2).*exp( - STATES(:,1)./CONSTANTS(:,44))))./140.000;
    ALGEBRAIC(:,28) =  STATES(:,6).*ALGEBRAIC(:,27);
    ALGEBRAIC(:,10) =  CONSTANTS(:,44).*log(STATES(:,2)./STATES(:,3));
    ALGEBRAIC(:,29) = ( (( CONSTANTS(:,13).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,14))).*(STATES(:,1) - ALGEBRAIC(:,10)))./(1.00000+exp(( ((STATES(:,1)+10.0000) - ALGEBRAIC(:,10)).*2.00000)./CONSTANTS(:,44)));
    ALGEBRAIC(:,25) =  (( STATES(:,5).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,7))).*CONSTANTS(:,6).*(STATES(:,1) - ALGEBRAIC(:,10));
    ALGEBRAIC(:,38) =  (( 0.0100000.*CONSTANTS(:,27).*(STATES(:,1) - 50.0000))./( CONSTANTS(:,44).*(1.00000 - exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))))).*( STATES(:,3).*exp(50.0000./CONSTANTS(:,44)) -  STATES(:,2).*exp((  - 1.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))).*STATES(:,10).*STATES(:,11).*STATES(:,12);
    ALGEBRAIC(:,42) = (ALGEBRAIC(:,29)+ALGEBRAIC(:,28)+ALGEBRAIC(:,25)+ALGEBRAIC(:,38)) -  2.00000.*ALGEBRAIC(:,33);
    ALGEBRAIC(:,26) = ALGEBRAIC(:,23)+ALGEBRAIC(:,25);
    ALGEBRAIC(:,31) =  0.500000.*CONSTANTS(:,44).*log(CONSTANTS(:,17)./STATES(:,7));
    ALGEBRAIC(:,32) =  CONSTANTS(:,16).*(STATES(:,1) - ALGEBRAIC(:,31));
    ALGEBRAIC(:,37) =  (( 4.00000.*CONSTANTS(:,27).*(STATES(:,1) - 50.0000))./( CONSTANTS(:,44).*(1.00000 - exp((  - 1.00000.*(STATES(:,1) - 50.0000).*2.00000)./CONSTANTS(:,44))))).*( STATES(:,7).*exp(100.000./CONSTANTS(:,44)) -  CONSTANTS(:,17).*exp((  - 2.00000.*(STATES(:,1) - 50.0000))./CONSTANTS(:,44))).*STATES(:,10).*STATES(:,11).*STATES(:,12);
    ALGEBRAIC(:,43) = ALGEBRAIC(:,37)+ALGEBRAIC(:,38)+ALGEBRAIC(:,40);
    ALGEBRAIC(:,44) = ( (( 2.00000.*1.00000.*CONSTANTS(:,48).*CONSTANTS(:,3))./( 1.00000.*CONSTANTS(:,39))).*STATES(:,14).*power(STATES(:,7), CONSTANTS(:,40)))./(power(STATES(:,7), CONSTANTS(:,40))+power(CONSTANTS(:,36), CONSTANTS(:,40)));
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end
