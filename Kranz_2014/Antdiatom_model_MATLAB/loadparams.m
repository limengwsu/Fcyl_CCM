function p = loadparams(pFile)


fp = fopen(pFile,'r');

if fp==-1
   error(['File ' pFile ' not found or permission denied.']);
end

i=1;


while 1 
    line = fgetl(fp);
    if (strcmp('________',line))
        break
    end
    raw{i} = line;
    i = i+1;
end

p = struct('N',[],'TC',[],'pH',[],'pHc',[],'pHp',[],'kcf',[],'kpf',[],'kyf',[],...
    'mRub',[],'Vm_Bc',[],'Km_Bc',[],'Vm_Bp',[],'Km_Bp',[],'fc_sm',[],'fb_sm',[]...
    ,'fc_p',[],'fb_p',[],'fc_y',[]);

for i = 1:length(raw)
    [id] = sscanf(raw{i},'%s[\t]%*s');
    
    switch id
        case 'cells/mL'
            p.N = sscanf(raw{i},'%*s %f');        %cells/mL
        case 'temp'
            p.TC = sscanf(raw{i},'%*s %f');        %temp
        case 'DIC'
            p.DIC = sscanf(raw{i},'%*s %f');      %uM; convert to mol/cm3
            p.DIC = p.DIC.*1E-9;
        case 'pHe'
            p.pH = sscanf(raw{i},'%*s %f');        %extracellular pH
        case 'pHc'
            p.pHc = sscanf(raw{i},'%*s %f');        %cytoplasmic pH
        case 'pHp'
            p.pHp = sscanf(raw{i},'%*s %f');        %chloroplast stroma pH
        case 'ksf'
            p.ksf = sscanf(raw{i},'%*s %e');        %surface CA activity
        case 'kcf'
            p.kcf = sscanf(raw{i},'%*s %e');        %cytoplasmic CA activity
        case 'kpf'    
            p.kpf = sscanf(raw{i},'%*s %e');        %chloroplast stroma CA activity
        case 'kyf'    
            p.kyf = sscanf(raw{i},'%*s %e');        %pyrenoid CA activity
        case 'mRub'    
            p.mRub = sscanf(raw{i},'%*s %e');        %mols of RubisCO molecules
        case 'Vm_Bc'    
            p.Vm_Bc = sscanf(raw{i},'%*s %e');        %maximum uptake rate of cytoplasmic HCO3- transporter in mol/cell/s
        case 'Km_Bc'    
            p.Km_Bc = sscanf(raw{i},'%*s %e');        %Km for cytoplasmic HCO3- transporter (mol/cm3)
        case 'Vm_Bp'    
            p.Vm_Bp = sscanf(raw{i},'%*s %e');        %maximum uptake rate of chloroplast pump HCO3- transporter in mol/cell/s
        case 'Km_Bp'    
            p.Km_Bp = sscanf(raw{i},'%*s %e');        %Km for chloroplast pump HCO3- transporter (mol/cm3)
        case 'reff'
            p.reff =  sscanf(raw{i},'%*s %e');        %effective cell radius (um convert to cm)
            p.reff = p.reff.*1E-4;
        case 'fc_sm'    
            p.fc_sm = sscanf(raw{i},'%*s %e');        %MTC for CO2 across cytoplasmic membrane
        case 'fb_sm'    
            p.fb_sm = sscanf(raw{i},'%*s %e');        %MTC for HCO3- across cytoplasmic membrane. generally set to zero
        case 'fc_p'    
            p.fc_p = sscanf(raw{i},'%*s %e');        %MTC for CO2 into chloroplast, from Hopkinson et al. 2011
        case 'fb_p'    
            p.fb_p = sscanf(raw{i},'%*s %e');        %MTC for HCO3- into chloroplast, generally set to zero
        case 'fc_y'    
            p.fc_y = sscanf(raw{i},'%*s %e');        %MTC for CO2 into pyrenoid, from Hopkinson et al. 2011
    end
end

%check to make sure all the fields were filled
fn = fieldnames(p);
for i = 1:length(fn)
    if isempty(p.(fn{i}))
        error('Error in load_params, missing pameter: %s\n',fn{i})
    end
end

%define pameters for CCM mechanistic model
p.T  = p.TC + 273.15;
p.S  = 35;          %salinity
p.H  = 10^-p.pH;    %H+ of extracellular solution  
p.Hc  = 10^-p.pHc;  %H+ concentration in the cytoplasm;
p.Hp = 10^-p.pHp;   %H+ concentration in the stroma/pyrenoid

%equilibrium constants for DIC, water, and borate
p.K1 = 10^-((3633.86./p.T) - 61.2172 + 9.6777*log(p.T) - 0.011555 .* p.S + 0.0001152 .* (p.S^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
p.K2 = 10^((-471.78/p.T) - 25.929 + 3.16967 * log(p.T) + 0.01781 * p.S - 0.0001122*p.S^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
p.Kw = exp(148.96502 - (13847.26 ./ p.T) - 23.6521 .* log(p.T) + (p.S^.5).*((118.67 ./ p.T) - 5.977 + 1.0495 .* log(p.T)) - 0.01615 .* p.S); %ion product of water, CO2 methods DOE 1994
p.Kb = exp(((-8966.9 - 2890.53 .* p.S^0.5 - 77.942 .* p.S + 1.728 .* p.S^1.5 - 0.0996 .* p.S^2)./p.T) + 148.0248 + 137.1942 .* p.S^0.5 + 1.62142 .* p.S...
    -(24.4344 + 25.085 .* p.S^0.5 + 0.2474 .* p.S) .* log(p.T) + 0.053105 .* p.S^0.5 .* p.T); %Kb boric acid/borate equilibrium constant
p.bfrac_e = (1/ (1+ p.K2./p.H));       % fraction of "B" pool that is HCO3- in extracellular solution
p.bfrac_i = (1/ (1+ p.K2./p.Hc));       % fraction of "B" pool that is HCO3- in cytoplasm
p.bfrac_x = (1/ (1+ p.K2./p.Hp));       % fraction of "B" pool that is HCO3-   in chloroplast and pyrenoid

%kinetic pameters
p.kp1 = exp(1246.98 - 6.19E4 ./ p.T - 183 .* log(p.T));         % CO2 + H2O -> H+ + HCO3- Johnson 1982; as presented in Zeebe and Wolf-Gladrow
p.kp4 = 4.70E7 .* exp(-23.2 ./ (8.314E-3 .* p.T));                % CO2 + OH- -> HCO3- Johnson 1982
p.km1 = p.kp1./p.K1;                                            % H+ + HCO3- -> CO2 + H2O
p.km4 = p.kp4.*p.Kw./p.K1;                                    % HCO3- -> CO2  + OH-
p.kph5 = 5.0E10;                                                    % CO32- + H+ -> HCO3-; /(M s); Zeebe and Wolf Gladrow
p.kmh5 = p.kph5 .* p.K2;                                        % HCO3- -> CO32- + H+;
p.kpoh5 = 6.0E9;                                                    % HCO3- + OH- -> CO32- + H2O; Zeebe and Wolf-Gladrow
p.kmoh5 = p.kpoh5 .* p.Kw ./ p.K2;                            % CO32- + H2O -> HCO3- + OH-;

%diffusion coefficients
DCF = diff_coef(p.TC,p.S, 1);        %diffusion coefficients (in cm2/s) at current temp, salinity, final variable is pressure: assume 1 atm.
p.Dc  = DCF(3);
fprintf(1,'CO2 diffusivity cm2/s: %5.3E\n',p.Dc);
p.Db  = DCF(8);       %diffusivity of HCO3- in cm2/s

%CO2 hydration/HCO3- dehydration rates 
p.kuf = p.kp1 + p.kp4.*(p.Kw./p.H);     %CO2 hydration rate in bulk solution

p.kur = p.bfrac_e.* p.kuf .* (p.H./p.K1);           %HCO3- dehyration rate in bulk solution
p.ksr = p.bfrac_e.* p.ksf .* (p.H./p.K1);          %HCO3- dehyration rate in surface layer (cm3/s)
p.kcr = p.bfrac_i.* p.kcf .* (p.Hc./p.K1);           %HCO3- dehyration rate in cytoplasm (/s)
p.kpr = p.bfrac_x.* p.kpf .* (p.Hp./p.K1);           %HCO3- dehyration rate in chloroplast stroma (/s)
p.kyr = p.bfrac_x.* p.kyf .* (p.Hp./p.K1);           %HCO3- dehyration rate in pyrenoid (/s)

%volumes Large diatom
p.Ve = 1;                               %volume of solution (cm3);
p.Vc = 8.3E-9;                         %total volume of a single cell w/ 15 um radius (cm3);
p.Vp = 8.3E-10;                         %chloroplast stroma volume (estimated as 10% of cell volume)
p.r_bl = 0.1E-4;                          %1 um effective surface layer
p.Vs = (4/3)*pi*((p.reff+p.r_bl)^3 - p.reff^3); % effective surface volume in cm3
Rpyr = 2.5/1E4;                         %radius of the pyrenoid, scaled from Pt pyrenoid size estimate and RubisCO content
Npyr = 1;                             %number of pyrenoids per cell, most seem to have 1 so just assume that
p.Vy = Npyr.*(4.*pi./3).*(Rpyr^3);        %volume of the pyrenoid (cm3), as estimated in Hopkinson et al 2011 from BCA clusters;
fprintf(1,'Vy /s: %5.3E\n',p.Vy);
%enzyme kinetics
p.kcat_R = 0.12;                   %Fcyl Rubisco turnover rate at 0.5 C (/s) according to Young et al. in prep TAble 1
p.Km_R   = 16 .* 1E-9;            %Fcyl Rubisco Km at 0.5 C(uM converted to mol/cm3) according to Kranz et al. in prep

%additional mass transfer coefficients
p.fc_bl = 4*pi * p.Dc * (p.reff + p.r_bl); %parameter for diffusive CO2 flux to cell surface layer
p.fb_bl = 4*pi * p.Db * (p.reff + p.r_bl); %parameter for diffusive HCO3- flux to cell surface layer
p.fb_y = Npyr.*4.*pi.*p.Db.*Rpyr;                     %MTC for HCO3- into pyrenoid.
end


