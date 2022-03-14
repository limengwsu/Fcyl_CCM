function Flux = Cifluxes(Yss,p)
%calculate net diffusional fluxes and CO2/HCO3 interchange at steady state

ce = Yss(1,:);      %CO2 in extracellular solution
be = Yss(2,:);      %HCO3- in extracellular solution
cs = Yss(3,:);      %CO2 at surface
bs = Yss(4,:);      %HCO3- at surface
cc = Yss(5,:);      %CO2 in cytoplasm
bc = Yss(6,:);      %HCO3- in cytoplasm
cp = Yss(7,:);      %CO2 in chloroplast stroma
bp = Yss(8,:);      %HCO3- in chloroplast stroma
cy = Yss(9,:);      %CO2 in pyrenoid
by = Yss(10,:);      %HCO3- in pyrenoid

%Flux.NetCO2influx = (p.fc_bl + p.ksf) .*(ce - cs);
 
%mass transfer vector:
vMT = [p.fc_bl; p.fb_bl; p.fc_sm; p.fb_sm; p.fc_p; p.fb_p; p.fc_y; p.fb_y];
v1 = vMT.*Yss(1:8,1);
v2 = vMT.*Yss(3:10,1);
t = [v1'; v2'];
Flux.Diff = t(:);


%CO2 hyd/dehydration fluxes
Flux.Hyd(1,1)  = (p.kuf .*ce).* (p.Ve./p.N);
Flux.Hyd(2,1)  = (p.kur .*be).* (p.Ve./p.N);
Flux.Hyd(3,1)  = (p.ksf .*cs);
Flux.Hyd(4,1)  = (p.ksr .*bs);
Flux.Hyd(5,1)  = (p.kcf .*cc).* p.Vc;
Flux.Hyd(6,1)  = (p.kcr .*bc).* p.Vc;
Flux.Hyd(7,1)  = (p.kpf .*cp).* p.Vp;
Flux.Hyd(8,1)  = (p.kpr .*bp).* p.Vp;
Flux.Hyd(9,1)  = (p.kyf .*cy).* p.Vy;
Flux.Hyd(10,1) = (p.kyr .*by).* p.Vy;
Flux.Hyd;

%photosynthesis and active transport
Flux.NetCO2influx = p.fc_bl .*(ce - cs) + (Flux.Hyd(4,1) - Flux.Hyd(3,1));
P = p.mRub.*(p.kcat_R .* cy)./(p.Km_R + cy);        %photosynthetic rate
Bup_c = (p.Vm_Bc .* be)./(p.Km_Bc + be);    %HCO3- uptake into cell
Bup_p = (p.Vm_Bp .* bc)./(p.Km_Bp + bc);    %chloroplast pump HCO3- transport rate
Flux.Active = [0; Bup_c; 0 ; Bup_p; P];   % more or less defined across memebranes (i.e. CO2 uptake across cyt, HCO3- across cyt, ...)
Flux.Active;


end