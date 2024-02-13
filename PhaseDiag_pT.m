%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  linewidth = 1; fontsize = 10;  markersize = 4;   
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);
AllEOS = {'PR','SRK','PTV','YR'};    figure(1);clf;hold on;box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define fluids to study
% Refrigerant = {'R134A','Emkarate RL32'};  MF1 = 0.9;
% Refrigerant = {'CO2','methane'};  MF1 = 0.9; 
% Refrigerant = {'methane','CO2'};  MF1 = 0.5;
% Refrigerant = {'CO2','ethane'};   MF1 = 0.1;
Refrigerant = {'Nitrogen','Emkarate RL32'};  MF1 = 0.3;
% Refrigerant = {'CO2','RENISO ACC HV'}; MF1 = 0.9;
% Refrigerant = {'RENISO ACC HV','CO2'};  MF1 = 0.1; 
% Refrigerant = {'CO2','R134A'};   MF1 = 0.5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define other parameters
Lplot = 1;                 % save the figure? 1 yes, 0 not
CubicEOS = AllEOS{4};      % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main program - Nothing needs to be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare global parameters
GL1 = GetGlobals(CubicEOS,{Refrigerant{1}});  % obtain fluid constants
GL2 = GetGlobals(CubicEOS,{Refrigerant{2}});  % obtain fluid constants
GL = GetGlobals(CubicEOS,Refrigerant);  % obtain fluid constants
Zi = [MF1,1 - MF1]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define the pressure range 
pc1 = GL1.presc;   
pc2 = GL2.presc;   
p_s = max([GL1.pres_r,GL2.pres_r])+1e4;
p_e = max([pc1,pc2]); 
Tratio = max(GL2.tempc/GL1.tempc,GL1.tempc/GL2.tempc);
pline = [linspace(p_s,p_e,15),linspace(p_e,Tratio*p_e,15)];
np = length(pline);
deleteindex = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if refpropm is avaiable and can be calculated
Tx_ref = zeros(1,np);       Ty_ref = zeros(1,np);    pline_ref = pline;  deleteindex_ref = [];
[MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,Zi); 
for ip = 1:np
    try
        Ty_ref(ip) = refpropm('T','P',pline(ip)/1000,'Q',1,Refrigerant{1},Refrigerant{2},MassFrac);
        Tx_ref(ip) = refpropm('T','P',pline(ip)/1000,'Q',0,Refrigerant{1},Refrigerant{2},MassFrac);
        if Tx_ref(ip) > Ty_ref(ip) || Tx_ref(ip) < 0  ||  Ty_ref(ip) < 0  ||  (Ty_ref(ip) - Tx_ref(ip)) < 0.3 
            error(' '); 
        end 
        if ip >= 2
            if Tx_ref(ip) < Tx_ref(ip-1)
                error(' ');
            end 
        end
    catch
        if ip ~= 1
            Ty_ref(ip) = Ty_ref(ip-1); 
            Tx_ref(ip) = Tx_ref(ip-1);
            if Tx_ref(ip) < Tx_ref(ip-1) && ip/np > 0.8
                deleteindex_ref = [deleteindex_ref,ip:np];
                break
            else
                deleteindex_ref = [deleteindex_ref,ip];
            end
        else
            deleteindex_ref = [deleteindex_ref,ip];
        end
    end
end
pline_ref(deleteindex_ref) = [];
Tx_ref(deleteindex_ref) = [];
Ty_ref(deleteindex_ref) = [];
if length(Tx_ref) > 2
    plot(Tx_ref,pline_ref/1e6,'r--',Ty_ref,pline_ref/1e6,'b--')
else
    fprintf('\n  =-- Refprop Calculation fail --=  \n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start the main loop
Tx = zeros(1,np);       Ty = zeros(1,np); 
for ip = 1:np
    if ip == 1
        ff = OilPropm('P',pline(ip),'Q',1,1,GL1,0,0);
        Tsat1 = ff.T_K;
        ff = OilPropm('P',pline(ip),'Q',1,1,GL2,0,0);
        Tsat2 = ff.T_K;
        T_y = max([Tsat1,Tsat2]) + 20; 
        T_x = min([Tsat1,Tsat2]) - 20; 
    else
        T_y = Ty(ip-1); % ffy.T_K;  
        T_x = Tx(ip-1); % ffx.T_K;
    end
    try
        ffy = OilPropm('P',pline(ip),'Q',1,Zi,GL,T_y,0);
        Ty(ip) = ffy.T_K;
        ffx = OilPropm('P',pline(ip),'Q',0,Zi,GL,T_x,0);
        Tx(ip) = ffx.T_K;
        if Tx(ip) > Ty(ip) || Tx(ip) < 0  ||  Ty(ip) < 0  ||  (Ty(ip) - Tx(ip)) < 0.3 
            error(' '); 
        end 
        if ip >= 2
            if Tx(ip) < Tx(ip-1)
                error(' ');
            end 
        end
    catch
        if ip ~= 1
            Ty(ip) = Ty(ip-1); 
            Tx(ip) = Tx(ip-1);
            if Tx(ip) < Tx(ip-1) && ip/np > 0.8
                deleteindex = [deleteindex,ip:np];
                break
            else
                deleteindex = [deleteindex,ip];
            end
        else
            deleteindex = [deleteindex,ip];
        end
    end
end
pline(deleteindex) = [];
Tx(deleteindex) = [];
Ty(deleteindex) = [];
plot(Tx,pline/1e6,'r-',Ty,pline/1e6,'b-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the results
% xlabel(['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(temp) ,' K'],'fontname','Arial','fontsize',fontsize)
ylabel('Pressure \itp\rm / MPa','fontname','Arial','fontsize',fontsize); 
xlabel('Temperature \itT\rm / K','fontname','Arial','fontsize',fontsize); 

if Lplot
    Material_mix = [Refrigerant{1},'_',Refrigerant{2},'_',num2str(MF1)];
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition',[0 0 18 8]);
    print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']); 
end