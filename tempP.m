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
