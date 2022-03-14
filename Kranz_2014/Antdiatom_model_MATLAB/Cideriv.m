function ydot = Cideriv(t,y,p)
%derivative function for Ci species in mechanistic model of the CCM

ce = y(1,1);
be = y(2,1);
cs = y(3,1);
bs = y(4,1);
cc = y(5,1);
bc = y(6,1);
cp = y(7,1);
bp = y(8,1);
cy = y(9,1);
by = y(10,1);

%calculate uptake and photosynthetic rates
Bupc  = (p.Vm_Bc.*be)./(p.Km_Bc + be);
Bupp  = (p.Vm_Bp.*bc)./(p.Km_Bp + bc);
P     = p.mRub.*(p.kcat_R .* cy)./(p.Km_R + cy);

dce = -p.kuf.*ce  + p.kur.*be  + (p.N./p.Ve).*p.fc_bl.*(cs - ce);
dbe =  p.kuf.*ce  - p.kur.*be  + (p.N./p.Ve).*p.fb_bl.*(bs - be);
dcs = (1./p.Vs).*(-p.ksf.*cs  + p.ksr.*bs  + p.fc_bl.*(ce - cs) + p.fc_sm.*(cc - cs));
dbs = (1./p.Vs).*( p.ksf.*cs  - p.ksr.*bs  + p.fb_bl.*(be - bs) + p.fb_sm.*(bc - bs) - Bupc);
dcc = -p.kcf.*cc  + p.kcr.*bc  + (1/p.Vc).*(p.fc_sm.*(cs - cc) + p.fc_p.*(cp - cc));
dbc =  p.kcf.*cc  - p.kcr.*bc  + (1/p.Vc).*(p.fb_sm.*(bs - bc) + p.fb_p.*(bp - bc) + Bupc - Bupp);
dcp = -p.kpf.*cp  + p.kpr.*bp  + (1/p.Vp).*(p.fc_p.*(cc - cp) + p.fc_y.*(cy - cp));
dbp =  p.kpf.*cp  - p.kpr.*bp  + (1/p.Vp).*(p.fb_p.*(bc - bp) + p.fb_y.*(by - bp) + Bupp);
dcy = -p.kyf.*cy  + p.kyr.*by  + (1/p.Vy).*(p.fc_y.*(cp - cy) - P);
dby =  p.kyf.*cy  - p.kyr.*by  + (1/p.Vy).*(p.fb_y.*(bp - by));

ydot = [dce; dbe;dcs; dbs; dcc; dbc; dcp; dbp; dcy; dby];

end

