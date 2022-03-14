function D = diff_coef(TC, Sal, P)
%Diffusion coefficients in sea water system constant calculation function
%Modified from Xinping Hu's code on his website which was: Modified from FORTRAN Code by Dr. Bernie Boudreau at the Dalhousie Univ.

% Varaibles in the function
%                  D(1) = H2O
%                  D(2) = O2
%                  D(3) = CO2
%                  D(4) = NH3
%                  D(5) = H2S
%                  D(6) = H3PO4
%                  D(7) = B(OH)3
%                  D(8) = HCO3-
%                  D(9) = CO3=
%                  D(10) = NH4+
%                  D(11) = HS-
%                  D(12) = NO3-
%                  D(13) = H2PO4-
%                  D(14) = HPO4=
%                  D(15) = PO4(---)
%                  D(16) = B(0H)4-
%                  D(17) = H+
%                  D(18) = OH-
%                  D(19) = Ca++
%                  D(20) = Mg++
%                  D(21) = Fe++
%                  D(22) = Mn++
%                  D(23) = SO4=
%                  D(24) = H4SiO4
%                  D(25) = CH4
%  Note: 1) enter S in ppt, T in deg. C and P in atm.
%  2) diffusion coefficients are in units of cm^2/s
%  3) H2O viscosity is in unit of centipoise

% Water density calculation
% Pure water density (kg/m3)
    RHO_Pw = (999.842594+6.793952*0.01*TC-9.09529*0.001*TC^2+1.001685*0.0001*TC^3-1.120083*10^-6*TC^4+...
        6.536332*10^-9*TC^5);

% Seawater density (kg/L)
    RHO=(RHO_Pw+(8.24493*0.1-4.0899*0.001*TC+7.6438*10^-5*TC^2-8.2467*10^-7*TC^3+5.3875*10^-9*TC^4)*Sal...
        +(-5.72466*10^-3+1.0227*10^-4*TC-1.6546*10^-6*TC^2)*Sal^1.5+4.8314*0.0001*TC^2)/1000;

% Viscosity calculation 
% Valid for 0<T<30 and 0<S<36, Calculated viscosity is in centipoise.
    Vis = visco(Sal, P, TC); % calculate in situ viscosity
    Vis0 = visco(0, 1, TC);  % calculate pure water viscosity at in situ temperature at 1 atmospheric pressure

%  Start calculations of diffusion coefficients in pure water at sample temperature.
    TK=TC+273.15;
    A = 12.5D-09*exp(-5.22D-04*P);
    B = 925.0*exp(-2.6D-04*P);
    T0 = 95.0 + 2.61D-02*P;
    D(1) = A*sqrt(TK)*exp(-B/(TK-T0))*1.0D+04;

%  Dissolved gases : from Wilke and Chang (1955)
%    note: 1) MW = molecular weight of water
%          2) VM = molar volumes (cm^3/mol) (Sherwood et al., 1975)
%  The factor PHI is reduced to 2.26 as suggested by Hayduk and Laudie (1974).
      PHI = 2.26;
      MW = 18.0;
      A1 = sqrt(PHI*MW)*TK/Vis0;
      
%Oxygen
      VMO2 = 25.6;
      D(2) = 7.4D-08*(A1/(VMO2^0.6));
%  CO2
      VMCO2 = 34.0;
      D(3) = 7.4D-08*(A1/(VMCO2^0.6));
%  NH3
      VMNH3 = 25.8;
      D(4) = 7.4D-08*(A1/(VMNH3^0.6));
%  H2S
      VMH2S = 32.9;
      D(5) = 7.4D-08*(A1/(VMH2S^0.6));   
%  CH4
      VMCH4 = 37.7;
      D(25) = 7.4D-08*(A1/(VMCH4^0.6));

%  The coefficients in pure water for the following species are
%  calculated by the linear functions of temperature (deg C)
%  found in Boudreau (1997, Springer-Verlag).

%  i.e. NO3-,HS-,H2PO4-,CO3=,SO4=,Ca++,Mg++,Mn++,Fe++,NH4+,H+ & OH-, HCO3-, HPO4=, PO4(3-)

      FAC = 1.0D-06;
      D(8) = (5.06 + 0.275*TC)*FAC;
      D(9) = (4.33 + 0.199*TC)*FAC;
      D(10) = (9.5 + 0.413*TC)*FAC;
      D(11) = (10.4 + 0.273*TC)*FAC;
      D(12) = (9.50 + 0.388*TC)*FAC;
      D(13) = (4.02 + 0.223*TC)*FAC;
      D(14) = (3.26 + 0.177*TC)*FAC;
      D(15) = (2.62 + 0.143*TC)*FAC;
      D(17) = (54.4 + 1.555*TC)*FAC;
      D(18) = (25.9 + 1.094*TC)*FAC;
      D(19) = (3.60 + 0.179*TC)*FAC;
      D(20) = (3.43 + 0.144*TC)*FAC;
      D(21) = (3.31 + 0.150*TC)*FAC;
      D(22) = (3.18 + 0.155*TC)*FAC;
      D(23) = (4.88 + 0.232*TC)*FAC;

% H3PO4 : Least (1984) determined D(H3PO4) at 25 deg C and 0 ppt S.
%         Assume that this value can be scaled by the Stokes-Einstein
%         relationship to any other temperature.
      D(6) = 0.87D-05;
      TS1 = 25.0;
      SS1 = 0.0;

      Vis_P = visco(SS1, 1, TS1); 

      D(6) = D(6)*Vis_P/298.15*TK/Vis0;
      
% B(OH)3 : Mackin (1986) determined D(B(OH)3) at 25 deg C and
%           about 29.2 ppt S.
%           Assume that this value can be scaled by the Stokes-Einstein
%           relationship to any other temperature.

      D(7) = 1.12D-05;
      TS2 = 25.0;
      SS2 = 29.2;
      
      Vis_B = visco(SS2, 1, TS2);

      D(7) = D(7)*Vis_B/298.15*TK/Vis0;
      
% B(OH)4 : No information on this species whatsoever! Boudreau and
%           Canfield (1988) assume it is 12.5% smaller than B(OH)3.

      D(16) = 0.875*D(7);


%  H4SiO4 : Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
%           and 36.1 ppt S.
%           Assume that this value can be scaled by the Stokes-Einstein
%           relationship to any other temperature.

      D(24) = 1.0E-05;
      TS3 = 25.0;
      SS3 = 36.1;
      Vis_Si = visco(SS3, 1, TS3);

      D(24) = D(24)*Vis_Si/298.15*TK/Vis0;

%  To correct for salinity, the Stokes-Einstein relationship is used.
%  This is not quite accurate, but is at least consistent.
%  Note in this code I did not include H3PO4, H4SiO4, B(OH)3 species

      FAC1 = Vis0/Vis;
      Diff_Coef=FAC1.*D;
      
return
