
sfc = 0.45;         % dimensionless
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












