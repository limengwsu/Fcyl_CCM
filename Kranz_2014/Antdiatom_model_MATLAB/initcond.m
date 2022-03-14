function Cinit = initcond(p)
%calculate initial conditions from DIC and pH

H = 10^-p.pH;
CO2 = p.DIC ./ (1 + (p.K1./H) + ((p.K1 * p.K2)./(H^2)));        %calculate total 13CO2 from DIC and pH, see Zeebe and Wolf-Gladrow CO2 in SW Chp1 p4
B   = p.DIC - CO2;                                              %pooled HCO3- + CO32-

Cinit = [CO2; B; CO2; B; CO2; B; CO2; B; CO2; B];

return;

end

