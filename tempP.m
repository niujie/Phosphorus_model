kp    = 10;             % Dimensionless, role of AM in movement of P ions
                        % role of mycorrhizal hyphae in the movement of P
                        % ions to the plant
n     = 0.625;          % Dimensionless, porosity
Zr    = 100;            % unit: cm, soil depth
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
        K_enzo = 0.04;          % unit: day^-1, maximum rate of enzymatic 
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
    UP = DEMstar;           % rate of P uptake
end

C_over_P_LF = 1300;     % unit: kg C kg P^-1, C/P ratio of leaf litter
% organic phosphorus stored in the biomass pool (Pv)
% LF is the rate at which litter (Carbon) is recycled back to the forest 
% floor, divided by C/P ratio of leaf litter gives the rate of litter (P)
% recycled back to the forest floor
dPvdt = UP - LF / C_over_P_LF;

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
dPMdt = ((1 - rr) / C_over_P_L - rh / C_over_P_HC) * DECL + (1 - rr) * ...
    DECH / C_over_P_H - MD / C_over_P_M + IMM;


% P in humus
% (C/P)HC must meet the condition: (1-rr)*(C/P)HC/(C/P)LF >= rh
dPHdt = rh * DECL / C_over_P_HC - DECH / C_over_P_H - ENZ_O;

ATM = 0.000016;     % unit: kg PI ha^-1 day^-1, PI deposited atmospherically
WE = kw * POCC;     % input due to weathering from mineral forms
kf = 0.0005;        % unit: day^-1, rate at which PI is fixed to POCC
FIX = kf * PI;      % available soil P, PI can be fixed to more occluded form
ERI = ke * PI;      % losses of PI due to erosion, and soil particles that 
                    % are transported during runoff
kl = 0.035;         % unit: day^-1, leaching loss constant

% b is the pore size distribution index
%                 b       c     beta
% Sand          4.05    11.1    12.1
% Loamy sand    4.38    11.7    12.7
% Sandy Load    4.90    12.8    13.8
% Loam          5.39    13.8    14.8
% Clay          11.4    25.8    26.8
beta = 2 * b + 4;
c = 2 * b + 3;
% vertical percolation with unit gradient
% Ks is the saturated hydraulic conductivity
Ls = Ks / (exp(beta*(1-sfc))-1) * (exp(beta*(s-sfc))-1);    % sfc < s <=1
% or Ls = K(s) = Ks * s^c

% a function of deep leaching losses
LE = Ls * kl * PI / (n * Zr * s);
% inorganic P
dPIdt = MIN + ATM + WE - kimm * PI * CM - FIX - UP - ERI - LE;

ER_O = ke * POCC;
% changes in the occluded P pool
dPOCCdt = FIX - WE - ENZ_I - ER_O;


