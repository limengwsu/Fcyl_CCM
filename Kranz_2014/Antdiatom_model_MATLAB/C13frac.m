function C13data = C13frac(del, Flux,p)
%natural abundance 13C fractionation prediction using fluxes determined
%during photosynthesis and known or assumed fractionation factors

%describe fractionation factors and convert to alpha notation
Ecdiff = 0;                  %diffusional fractionation of CO2
Ecup = 0;                    %fractionation associated with transport of CO2
Ebdiff = 0;                  %fractionation of HCO3 from diffusion
Ebup= 0;                     %fractionation of HCO3 from uptake 
Erub = -27;                  %rubisco fractionation %NOT REALLY SURE WHERE I GOT -22, THERE'S NO DATA ON DIATOMS. If you look at Tcherkez 2006 correlations between specificity and E would predict diatom RubisCO is ~-26, similar to -25 I used to use
Ecb = -13;                    %fractionation factor for production of HCO3 from CO2 (OLeary from Zeebe and Wolfe Gladrow 2001; at 25C).
Ebc = Ecb + CO2C13offset(p.T);                   %fractionation factor for production of CO2 from HCO3, not currently distinguishing between CA produced or uncatalyzed; must move toward equilibrium value 
%Ebc = -10;
%vectors associating each flux with a fractionation factor (Epsilon)
VepsDiff = [Ecdiff+del(1,1); Ecdiff; Ebdiff + del(2,1); Ebdiff;...         %diffusion from bulk solution to surface layer; include isotope composition of bulk CO2 and HCO3- because they're fixed
            Ecdiff; Ecdiff; Ebdiff; Ebdiff;...         %diffusion from surface layer to cytoplasm
            Ecdiff; Ecdiff; Ebdiff; Ebdiff;...         %diffusion from cytoplast to chloroplast stroma
            Ecdiff; Ecdiff; Ebdiff; Ebdiff];           %diffusion from chloroplast stroma to pyrenoid 
VepsHyd  = [Ecb; Ebc; Ecb; Ebc; Ecb; Ebc; Ecb; Ebc; Ecb; Ebc];
VepsAct  = [Ecup; Ebup; Ecup; Ebup; Erub];
        
%define matrix for linear system, rows are the systems dCs/dt, dBs/dt,etc
%all at steady state. Column represents constant terms

M = zeros(8,8);
%dCs/dt equation
M(1,1) = -(Flux.Diff(2,1) + Flux.Diff(5,1) + Flux.Hyd(3,1) + Flux.Active(1,1));
M(1,2) =  Flux.Hyd(4,1);
M(1,3) =  Flux.Diff(6,1);

%dBs/dt
M(2,1) = Flux.Hyd(3,1);
M(2,2) = -(Flux.Diff(4,1) + Flux.Diff(7,1) + Flux.Hyd(4,1) + Flux.Active(2,1));
M(2,4) = Flux.Diff(8,1);

%dCc/dt
M(3,1) = Flux.Diff(5,1) + Flux.Active(1,1);
M(3,3) = -(Flux.Diff(6,1) + Flux.Diff(9,1) + Flux.Hyd(5,1) + Flux.Active(3,1));
M(3,4) = Flux.Hyd(6,1);
M(3,5) = Flux.Diff(10,1);

%dBc/dt
M(4,2) = Flux.Diff(7,1) + Flux.Active(2,1);
M(4,3) = Flux.Hyd(5,1);
M(4,4) = -(Flux.Diff(8,1) + Flux.Diff(11,1) + Flux.Hyd(6,1) + Flux.Active(4,1));
M(4,6) = Flux.Diff(12,1);

%dCp/dt
M(5,3) = Flux.Diff(9,1) + Flux.Active(3,1);
M(5,5) = -(Flux.Diff(10,1) + Flux.Diff(13,1) + Flux.Hyd(7,1));
M(5,6) = Flux.Hyd(8,1);
M(5,7) = Flux.Diff(14,1);

%dBp/dt
M(6,4) = Flux.Diff(11,1) + Flux.Active(4,1);
M(6,5) = Flux.Hyd(7,1);
M(6,6) = -(Flux.Diff(12,1) + Flux.Diff(15,1) + Flux.Hyd(8,1));
M(6,8) = Flux.Diff(16,1);

%dCy/dt
M(7,5) = Flux.Diff(13,1);
M(7,7) = -(Flux.Diff(14,1) + Flux.Hyd(9,1) + Flux.Active(5,1));
M(7,8) = Flux.Hyd(10,1);

%dBy/dt
M(8,6) = Flux.Diff(15,1);
M(8,7) = Flux.Hyd(9,1);
M(8,8) = -(Flux.Diff(16,1) + Flux.Hyd(10));

%construct vector of constants
Mdiff = [1 -1  0  0 -1  1  0  0  0  0  0  0  0  0  0  0;...  %diffusion matrix (effect of diffusive fluxes on Cs, Bs, etc)
         0  0  1 -1  0  0 -1  1  0  0  0  0  0  0  0  0;...
         0  0  0  0  1 -1  0  0 -1  1  0  0  0  0  0  0;...
         0  0  0  0  0  0  1 -1  0  0 -1  1  0  0  0  0;...
         0  0  0  0  0  0  0  0  1 -1  0  0 -1  1  0  0;...
         0  0  0  0  0  0  0  0  0  0  1 -1  0  0 -1  1;...
         0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0;...
         0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1];

Mhyd = [0 0 -1  1  0  0  0  0  0  0;...          %hydration/dehydration matrix
        0 0  1 -1  0  0  0  0  0  0;...
        0 0  0  0 -1  1  0  0  0  0;...
        0 0  0  0  1 -1  0  0  0  0;...
        0 0  0  0  0  0 -1  1  0  0;...
        0 0  0  0  0  0  1 -1  0  0;...
        0 0  0  0  0  0  0  0 -1  1;...
        0 0  0  0  0  0  0  0  1 -1];
    
Mactive = [-1  0  0  0  0;...             %active flux matrix
            0 -1  0  0  0;...
            1  0 -1  0  0;...
            0  1  0 -1  0;...
            0  0  1  0  0;...
            0  0  0  1  0;...
            0  0  0  0 -1;...
            0  0  0  0  0];
            
V = -( Mdiff*(VepsDiff.*Flux.Diff) + Mhyd*(VepsHyd.*Flux.Hyd)  + Mactive*(VepsAct.*Flux.Active));  

C13data = M\V;
dC13Fix = C13data(7,1) + Erub;
fprintf(1,'13C-Cb: %f\n13C-Bb: %f\n13C-Cs: %f\n13C-Bs: %f\n13C-Cc: %f\n13C-Bc: %f\n13C-Cp: %f\n13C-Bp: %f\n13C-Cy: %f\n13C-By: %f\n13C-Fixed: %f\n',del, C13data,dC13Fix);
isofile='CCM_lowtemp_C13isotope.txt';
fiso = fopen(isofile,'w');
fprintf(fiso,'13C-Cb: %f\n13C-Bb: %f\n13C-Cs: %f\n13C-Bs: %f\n13C-Cc: %f\n13C-Bc: %f\n13C-Cp: %f\n13C-Bp: %f\n13C-Cy: %f\n13C-By: %f\n13C-Fixed: %f\n',del, C13data,dC13Fix);
fclose(fiso);
%IsotopeBalanceAcrossCM = ((C13data(3)+Ecdiff).*Flux.Diff(9,1)+(C13data(3)+Ecup).*Flux.Active(3,1) + (C13data(4)+Ebdiff).*Flux.Diff(11,1) + (C13data(4)+Ebup).*Flux.Active(4))...
%   ./((C13data(5)+Ecdiff).*Flux.Diff(10,1) + (C13data(6)+Ebdiff).*Flux.Diff(12,1))

return



