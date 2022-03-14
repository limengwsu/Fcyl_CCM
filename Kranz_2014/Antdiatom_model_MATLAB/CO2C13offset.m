function dC13CO2 = CO2C13offset(TK)
%calculates equilibrium CO2 fractionation offset from HCO3- as a function of temp 
%(Zeebe and Wolf Gladrow 2001; from Mook 1986, Zhang 1995 is similar)
dC13CO2 = 24.12 - (9866/TK);   

end

