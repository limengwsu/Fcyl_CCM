function out = CCM_lowtemp_large(infile)
%mechanistic model of the CCM to match P vs DIC data
%somethings wrong with surface layer  - can't get CO2 depleted - work on
%this
parfile = strcat(infile,'.par');
p = loadparams(parfile);          %model parameters

yinit = initcond(p);     %calculate initial inorganic carbon concentrations 

t = [];
Y = [];
%use ode solver to deterimine time course of CO2 and HCO3- behavior
    tspan = [0, 300];       %time spane in seconds
    options = odeset('RelTol', 1E-6, 'AbsTol', 1E-10,'MaxStep',5);         %set ode options to obtain smooth (non-oscillating or diverging) solution
    [t, Y] = ode15s(@Cideriv, tspan, yinit, options, p);
    Y = Y';                  %ode gives one row per time point, transpose structure in which each row is the time series of a single Cspecies
    t = t';

figure(1)    
subplot(2,5,1)
plot(t,Y(1,:),'r'), title('CO2_e');
subplot(2,5,6)
plot(t,Y(2,:),'r'), title('HCO3_e');
subplot(2,5,2)
plot(t,Y(3,:),'r'),title('CO2_s');
subplot(2,5,7)
plot(t,Y(4,:),'b'), title('HCO3_s');
subplot(2,5,3)
plot(t,Y(5,:),'b'),title('CO2_c');
subplot(2,5,8)
plot(t,Y(6,:),'b'), title('HCO3_c');
subplot(2,5,4)
plot(t,Y(7,:),'b'),title('CO2_p');
subplot(2,5,9)
plot(t,Y(8,:),'b'), title('HCO3_p');
subplot(2,5,5)
plot(t,Y(9,:),'b'),title('CO2_y');
subplot(2,5,10)
plot(t,Y(10,:),'b'), title('HCO3_y');


%compute determined P, Uptake rates
Yss(:,1) = mean(Y(:,end-5:end),2);

%calculate fluxes 
Flux = Cifluxes(Yss,p);
dC13bic = 0; %del 13 C of HCO3-
dC13CO2diff = CO2C13offset(p.T);
dC13CO2 = dC13bic + dC13CO2diff;
del = [dC13CO2; dC13bic];     %del13C of bulk solution CO2 and HCO3-
%del = [-10;0]
C13data = C13frac(del,Flux, p);

%write data out to files
dlmwrite('CCM_lowtemp_Fulloutput.txt',[t' Y'],'\t');           %full time course of Ci concentrations in all compartments

fluxfile = 'CCM_lowtemp_Fluxes.txt';
fid = fopen(fluxfile,'w');
% fprintf(fid,'CO2e\t HCO3e\t CO2c\t HCO3c\t CO2p\t HCO3p\t CO2y\t HCO3y\t P\t Bup_c\t Bup_y\t Diff_CO2_etoc\t Diff_HCO3_etoc\t Diff_CO2_ctop\t Diff_HCO3ctop\t Diff_CO2_ptoy\t Diff_HCO3_ptoy\t Hyd_e\t Dehyd_e\t Hyd_c\t Dehyd_c\t Hyd_p\t Dehyd_p\t Hyd_y\t Dehyd_y\n'); 
% fclose(fid);
% dlmwrite(ssfile,[Yss' P' Bup_c' Bup_p' Flux.cec' Flux.bec' Flux.ccp' Flux.bcp' Flux.cpy' Flux.bpy' Flux.hyd_e' Flux.dehyd_e' Flux.hyd_c' Flux.dehyd_c' Flux.hyd_p' Flux.dehyd_p' Flux.hyd_y' Flux.dehyd_y'],'-append','delimiter','\t');
%SOMETHING IS WRONG W. DIFFUSIOVE FLUXES
Label = {'F.C_e2bl','F.C_bl2e','F.B_e2bl','F.B_bl2e','F.C_bl2c','F.C_c2bl','F.B_bl2c','F.B_c2bl',...
         'F.C_c2p','F.C_p2c','F.B_c2p','F.B_p2c','F.C_p2y','F.C_y2p','F.B_p2y','F.B_y2p',...
         'F.Hyd_e','F.Dehyd_e','F.Hyd_bl','F.Dehyd_bl','F.Hyd_c','F.Dehyd_c','F.Hyd_p','F.Dehyd_p','F.Hyd_y','F.Dehyd_y',...
         'F.DiffCO2up','F.Cup_c','F.Bup_c','F.Cup_p','F.Bup_p','F.P'};
Data = [Flux.Diff; Flux.Hyd; Flux.NetCO2influx; Flux.Active];
for i =1:length(Data)
    fprintf(fid,'%s\t%e\n',Label{i},Data(i));
end
fclose(fid);
end

