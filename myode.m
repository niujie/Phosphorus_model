function dCPdt = myode(~, CP)

%% 
CL   = CP(1);
CH   = CP(2);
CM   = CP(3);
% PV   = CP(4);
PL   = CP(5);
PM   = CP(6);
PH   = CP(7);
PI   = CP(8);
POCC = CP(9);
s    = CP(10);

C_over_P_L = CL / PL;
C_over_P_H = CH / PH;
C_over_P_M = CM / PM;


%% soil moisture part

sstar = 0.27;       % dimensionless, soil moisture when plants close stomata
sW    = 0.17;       % dimensionless, soil moisture at wilting point
sh    = 0.14;       % dimensionless, hygroscopic soil moisture
EW    = 0.01;       % unit: cm day^-1, evaporation rate at wilting point
EMAX  = 0.15;       % unit: cm day^-1, maximum transpiration rate
if sh < s && s <= sW
    Es = EW * (s - sh) / (sW - sh);
else
    if s <= sstar
        Es = EW + (EMAX - EW) * (s - sW) / (sstar - sW);
    else
        if s <= 1
            Es = EMAX;
        end
    end
end

TMAX = 0.13;         % unit: cm day^-1, maximum evaporation rate
if 0 <= s && s <= sW
    Ts = 0;
else
    if s <= sstar
        Ts = (s - sW) / (sstar - sW) * TMAX;
    else
        if s <= 1
            Ts = TMAX;
        end
    end
end

n   = 0.625;          % Dimensionless, porosity
Zr  = 100;            % unit: cm, soil depth
sfc = 0.45;           % dimensionless
Ks  = 452;            % unit: cm day^-1, saturated hydraulic conductivity
% b is the pore size distribution index
%                 b       c     beta
% Sand          4.05    11.1    12.1
% Loamy sand    4.38    11.7    12.7
% Sandy Load    4.90    12.8    13.8
% Loam          5.39    13.8    14.8
% Clay          11.4    25.8    26.8
b = 5.0;
beta = 2 * b + 4;
% vertical percolation with unit gradient
% Ks is the saturated hydraulic conductivity
Ls = Ks / (exp(beta*(1-sfc))-1) * (exp(beta*(s-sfc))-1);    % sfc < s <=1
% c = 2 * b + 3;
% Ls = K(s) = Ks * s^c
% soil moisture dynamics, Ist is the rate of infiltration from rainfall
dsdt = (Ist - Es - Ts - Ls) / (n * Zr);

%% C part

if s <= sfc         % sfc is soil moisture at field capacity
    fd = s / sfc;   % nondimensional factor describes the effects of soil
                    % moisture on decomposition
else
    fd = sfc / s;
end

C_over_P_THRp = 65;     % unit: kg C kg P^-1, microbial C/P ratio where nutrients are limiting
C_over_P_THRm = 56;     % unit: kg C kg P^-1, microbial C/P ratio where carbon is limiting
if C_over_P_M < C_over_P_THRm
    phi = 1;        % parameter reduces decomposition rates when nutrient poor
else
    if C_over_P_M < C_over_P_THRp
        phi = (C_over_P_THRp - C_over_P_M) / (C_over_P_THRp - C_over_P_THRm);
    else
        phi = 0;
    end
end

kd = 0.0004;        % unit: ha Mg C^-1 day^-1, litter decomposition rate
% decomposition rate of leaf litter
% CM is carbon in microbial biomass
% CL is carbon in litter
DECL = phi * fd * kd * CM * CL;


rh = 0.18;          % dimensionless, fraction of decomposing litter that 
                    % undergoes humification, isohumic coefficient
kh = 0.000003;      % unit: ha Mg C^-1 day^-1, humus decomposition rate
% decomposition rate of the humus pool
DECH = phi * fd * kh * CM * CH;
% carbon in humus
dCHdt = rh * DECL - DECH;


if s <= sfc                 % sfc is soil moisture at field capacity
    fm = (sfc - s) / sfc;   % the control that soil moisture, s, has on the 
                            % death rate
else
    fm = (s - sfc) / (1 - sfc);
end
kmax = 0.001;       % unit: day^-1, maximum microbial death rate
kopt = 0.000008;    % unit: day^-1, optimal microbial death rate
% mortality rate of microbial biomass
MD = max(kmax*fm, kopt) * CM;
rr = 0.6;           % dimensionless, fraction of decomposing C lost to respiration
% carbon in microbial biomass
dCMdt = (1-rh-rr)*DECL + (1-rr)*DECH - MD;

% the amount of leaf litter, dead roots and other plant residues that are
% present on the forest floor
% LF is the rate of litter recycled back to the forest floor
% be replaced by CLM column carbon state variable, 
% totlitc, total litter carbon, unit: gC/m2
dCLdt = LF + MD - DECL;


%% P part

kp    = 10;             % Dimensionless, role of AM in movement of P ions
                        % role of mycorrhizal hyphae in the movement of P
                        % ions to the plant
Pstar = 0.054;          % unit: kg PI ha^-1, amount of PI required before 
                        % it becomes limiting


DEMstar = kp * Ts * Pstar / (n*Zr*s);   % P demand when the amount of 
                                        % PI in the soil can satisfy
                                        % the demand

% the following logic corresponds to eqs. (11) -- (17)
if PI < Pstar           % PI, available soil, inorganic P
    % PI / (n*Zr*s) is the concentration of PI in soil water for a soil
    % with a root depth, Zr, porosity, n, and relative soil moisture
    % content, s
    UPp = kp * Ts * PI / (n*Zr*s);  % passive uptake
    
    if DEMstar - UPp <= 0       % i.e. PI >> P*
        UPa = 0;                % active uptake through enzymatic release
    else
        k_enzi = 0.014;         % unit: day^-1, maximum rate of enzymatic 
                                % dissolution from Pocc, occluded P pool
        ENZ_I_max = k_enzi * Pocc;
        k_enzo = 0.04;          % unit: day^-1, maximum rate of enzymatic 
                                % dissolution from PH, humus P pool
        ENZ_O_max = k_enzo * PH;
        
        if DEMstar - UPp <= ENZ_O_max + ENZ_I_max
            UPa = DEMstar - UPp;
        else
            if DEMstar - UPp <= 0
                ENZ_O = 0;      % enzymatic breakdown of organic P by releasing
                                % extracellular enzymes
            else
                if DEmstar - UPp <= ENZ_O_max
                    ENZ_O = DEMstar - UPp;
                else
                    ENZ_O = ENZ_O_max;
                end
            end
            if DEMstar - UPp <= ENZ_O_max
                ENZ_I = 0;      % the dissolution of P-containing minerals 
                                % through a combination of soil
                                % acidification and the release of metal
                                % complexing agents (predominantly organic 
                                % acid anions)
            else
                if DEMstar - UPp <= ENZ_O_max + ENZ_I_max
                    ENZ_I = DEMstar - UPp - ENZ_O_max;
                else
                    ENZ_I = ENZ_I_max;
                end
            end
            UPa = ENZ_O + ENZ_I;
        end
    end
    
    UP = UPp + UPa;     % rate of P uptake
else
    UPp = DEMstar;
    UP = UPp;           % rate of P uptake
end

C_over_P_LF = 1300;     % unit: kg C kg P^-1, C/P ratio of leaf litter
% organic phosphorus stored in the biomass pool (Pv)
% LF is the rate at which litter (Carbon) is recycled back to the forest 
% floor, divided by C/P ratio of leaf litter gives the rate of litter (P)
% recycled back to the forest floor
dPVdt = UP - LF / C_over_P_LF;

% P in plant litter
% MD, mortality rate of microbial biomass, 
% DECL, decomposition rate of leaf litter,
% both are calculated in C cycle
% (C/P)_M is the carbon to phosphorus ratio of the microbial biomass
% (C/P)_L is the carbon to phosphorus ratio of the decomposing leaf litter
dPLdt = LF / C_over_P_LF + MD / C_over_P_M - DECL / C_over_P_L;

kimm = 0.015;   % unit: ha Mg C^-1 day^-1, rate at which microbes immobilize PI
% imobilization
IMM = (1 - phi) * (rr* (DECL / C_over_P_L + DECH / C_over_P_H) + kimm * PI * CM);
% mineralization
MIN = phi * rr * (DECL / C_over_P_L + DECH / C_over_P_H);

% P in microbial biomass
% rr, fraction of decomposing C lost to respiration
% rh, fraction of decomposing litter that 
% undergoes humification, isohumic coefficient
% C_over_P_HC, the recalcitrant organic matter that microbes are unable to
% decompose (i.e., the fraction going into the humus pool) will also have
% a relatively constant C/P ratio
% (1-rr)*(C/P)_HC/(C/P)_LF >= rh
% just put a coefficient of 1.5 here temporarily
C_over_P_HC = 1.5 * rh * C_over_P_LF / (1 - rr);
dPMdt = ((1 - rr) / C_over_P_L - rh / C_over_P_HC) * DECL + (1 - rr) * ...
    DECH / C_over_P_H - MD / C_over_P_M + IMM;


% P in humus
% (C/P)HC must meet the condition: (1-rr)*(C/P)HC/(C/P)LF >= rh
dPHdt = rh * DECL / C_over_P_HC - DECH / C_over_P_H - ENZ_O;

ATM = 0.000016;     % unit: kg PI ha^-1 day^-1, PI deposited atmospherically
WE = kw * POCC;     % input due to weathering from mineral forms
kf = 0.0005;        % unit: day^-1, rate at which PI is fixed to POCC
FIX = kf * PI;      % available soil P, PI can be fixed to more occluded form
ke = 0.001;         % NOT SURE
ERI = ke * PI;      % losses of PI due to erosion, and soil particles that 
                    % are transported during runoff

kl = 0.035;         % unit: day^-1, leaching loss constant
% a function of deep leaching losses
LE = Ls * kl * PI / (n * Zr * s);
% inorganic P
dPIdt = MIN + ATM + WE - kimm * PI * CM - FIX - UP - ERI - LE;

ER_O = ke * POCC;
% changes in the occluded P pool
dPOCCdt = FIX - WE - ENZ_I - ER_O;

%% construct ODE

dCPdt = zeros(size(CP));
dCPdt(1) = dCLdt;
dCPdt(2) = dCHdt;
dCPdt(3) = dCMdt;
dCPdt(4) = dPVdt;
dCPdt(5) = dPLdt;
dCPdt(6) = dPMdt;
dCPdt(7) = dPHdt;
dCPdt(8) = dPIdt;
dCPdt(9) = dPOCCdt;
dCPdt(10) = dsdt;