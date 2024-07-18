%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  linewidth = 1; fontsize = 10;  markersize = 4;   
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);
AllEOS = {'PR','SRK','PTV','YFR'};     figure(1);clf;hold on;box on;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define fluids to study
Refrigerant = {'R134A','Emkarate RL32'};  pres_kPa = 3.5 * 1e3;  
% Refrigerant = {'Emkarate RL32','R134A'};  pres_kPa = 3.5 * 1e3;  
% Refrigerant = {'CO2','methane'};  pres_kPa =  3 * 1e3;  
% Refrigerant = {'methane','CO2'};  pres_kPa =  5 * 1e3;  
% Refrigerant = {'Nitrogen','Emkarate RL32'};  pres_kPa = 2 * 1e3;  
% Refrigerant = {'Emkarate RL32','Nitrogen'};  pres_kPa = 3.2 * 1e3;  
% Refrigerant = {'CO2','RENISO ACC HV'}; pres_kPa = 1 * 1e3;  
% Refrigerant = {'RENISO ACC HV','CO2'}; pres_kPa = 2 * 1e3;  
% Refrigerant = {'propane','R134A'};    pres_kPa = 1 * 1e3; 
% Refrigerant = {'CO2','ethane'}; pres_kPa = 3 * 1e3;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define other parameters
Lplot = 1;                 % save the figure? 1 yes, 0 not
CubicEOS = AllEOS{4};      % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main program - Nothing needs to be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare global parameters
CubicEOS = AllEOS{3}; 
GL1 = GetGlobals(CubicEOS,{Refrigerant{1}});  % obtain fluid constants
GL2 = GetGlobals(CubicEOS,{Refrigerant{2}});  % obtain fluid constants
GL = GetGlobals(CubicEOS,Refrigerant);  % obtain fluid constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define the pressure range 
try
    ff = OilPropm('all','P',pres_kPa,'Q',1,1,GL1,0,0);  Tsat1 = ff.T_K;  lve1 = 1;
catch
    Tsat1 = GL1.tempc;    lve1 = 0;
end
try
    ff = OilPropm('all','P',pres_kPa,'Q',1,1,GL2,0,0);  Tsat2 = ff.T_K;   lve2 = 1;
catch
    Tsat2 = GL2.tempc;    lve2 = 0;
end
T_s = min([Tsat1,Tsat2]); T_e = max([Tsat1,Tsat2]); 
Tline = linspace(T_s,T_e,40);
nT = length(Tline);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialization
Tratio = max(GL2.tempc/GL1.tempc,GL1.tempc/GL2.tempc);
if  lve1 == 1  &&  lve2 == 1  
    if Tsat2 > Tsat1,  x1 = ones(1,nT);     y1 = ones(1,nT);   x1(nT) = 0; y1(nT) = 0; end
    if Tsat2 < Tsat1,  x1 = zeros(1,nT);     y1 = zeros(1,nT);   x1(nT) = 1; y1(nT) = 1; end
    deleteindex = [];
    if Tratio > 1.5 algorithm = 2; else algorithm = 1;end
elseif  lve1 == 0  && lve2 == 1  
    if Tsat2 > Tsat1,  x1 = zeros(1,nT);     y1 = zeros(1,nT);  deleteindex = 1;  end
    if Tsat2 < Tsat1,  x1 = zeros(1,nT);     y1 = zeros(1,nT);  deleteindex = nT;  end
    algorithm = 2;
elseif  lve1 == 1  && lve2 == 0  
    if Tsat2 > Tsat1,  x1 = ones(1,nT);     y1 = ones(1,nT);  deleteindex = nT; end
    if Tsat2 < Tsat1,  x1 = ones(1,nT);    y1 = ones(1,nT);  deleteindex = 1; end
    algorithm = 2;
else
    error('Input pressure two high.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if refpropm is avaiable and can be calculated
nx = 40; deleteindex_ref = [];
x1_ref = linspace(0,1,nx);       y1_ref = linspace(0,1,nx); 
Tx1_ref = zeros(1,nx);           Ty1_ref = zeros(1,nx);   
for ix = 1:nx
    MF1 = x1_ref(ix);
    Zi = [MF1,1 - MF1]';
    [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,Zi);  
    try
        Ty1_ref(ix) = refpropm('T','P',pres_kPa,'Q',1,Refrigerant{1},Refrigerant{2},MassFrac);
        Tx1_ref(ix) = refpropm('T','P',pres_kPa,'Q',0,Refrigerant{1},Refrigerant{2},MassFrac);
        if Ty1_ref(ix) < Tx1_ref(ix) || Ty1_ref(ix) < 0  ||  Tx1_ref(ix) < 0
            error(' '); 
        end 
    catch
        deleteindex_ref = [deleteindex_ref,ix];
    end
end
Ty1_ref(deleteindex_ref) = [];
Tx1_ref(deleteindex_ref) = [];
x1_ref(deleteindex_ref) = [];
y1_ref(deleteindex_ref) = [];
if length(x1_ref) > 2
    plot(x1_ref,Tx1_ref,'b--',y1_ref,Ty1_ref,'r--')
else
    fprintf('\n  =-- Refprop Calculation fail --=  \n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start the main loop
if algorithm == 1   % mainly for azeotropic mixture
    nx = 40;
    x1 = linspace(0,1,nx);       y1 = linspace(0,1,nx); 
    Tx1 = zeros(1,nx);           Ty1 = zeros(1,nx);   
    Tx1(1) = Tsat2;              Ty1(1) = Tsat2; 
    Tx1(nx) = Tsat1;             Ty1(nx) = Tsat1; 
    for ix = 2:nx-1
        MF1 = x1(ix);
        Zi = [MF1,1 - MF1]';
        [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,Zi); 
        if ix == 2
            T_y = Ty1(ix-1);   T_x = Tx1(ix-1);
        else
            T_y = ffy.T_K;    T_x = ffx.T_K;
        end
        try
            ffy = OilPropm('all','P',pres_kPa,'Q',1,MassFrac,GL,T_y,0);
            Ty1(ix) = ffy.T_K;
            ffx = OilPropm('all','P',pres_kPa,'Q',0,MassFrac,GL,T_x,0);
            Tx1(ix) = ffx.T_K;
            if Ty1(ix) < Tx1(ix) || Ty1(ix) < 0  ||  Tx1(ix) < 0
                error(' '); 
            end 
        catch
            deleteindex = [deleteindex,ix];
        end
    end
    Ty1(deleteindex) = [];
    Tx1(deleteindex) = [];
    x1(deleteindex) = [];
    y1(deleteindex) = [];
    plot(x1,Tx1,'b-',y1,Ty1,'r-')
elseif algorithm == 2
    dx = 0.01;  ddx0 = 0.003;     
    firstP = 0;
    for iT = 2:nT-1
%         if iT == 2 || firstP == 0
%             if Tsat2 > Tsat1  
%                 MF1 = 1 - dx;
%             else
%                 MF1 = dx;
%             end 
%         else
%             MF1 = (ffo.MoleF_Li(1) + ffo.MoleF_Vi(1))/2  + 2 * ((Tsat1 > Tsat2) -0.5) * dx;
            MF1 = (Tline(iT) - Tsat2) / (Tsat1 - Tsat2);
%         end
        loopcount = 0;
        while 1 
            loopcount = loopcount + 1;
            if loopcount == 1
                ddx = ddx0;   
            elseif loopcount == 5
                ddx = ddx0 * 5; 
            elseif loopcount == 10
                ddx = ddx0 * 25; 
            end
            Zi = [MF1,1 - MF1]';
            [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,Zi); 
            ff = OilPropm('all','T',Tline(iT),'P',pres_kPa,MassFrac,GL,0,0);
            if strcmpi(ff.Phase,'V')
                MF1 = MF1 + 2 * ((Tsat1 > Tsat2) -0.5) * ddx;
            elseif strcmpi(ff.Phase,'L')
                MF1 = MF1 - 2 * ((Tsat1 > Tsat2) -0.5) * ddx;
            elseif strcmpi(ff.Phase,'LV')
                firstP = 1;
                ffo = ff;
                break
            end
            if MF1 > 1 || MF1 < 0 || loopcount >= 15
                deleteindex = [deleteindex,iT];
                break;
            end
        end
        x1(iT) = ff.MoleF_Li(1);
        y1(iT) = ff.MoleF_Vi(1);
    end
    Tline(deleteindex) = [];
    x1(deleteindex) = [];
    y1(deleteindex) = [];
    plot(x1,Tline,'b-',y1,Tline,'r-')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the results
xlabel(['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(pres_kPa/1e3) ,' MPa'],'fontname','Arial','fontsize',fontsize)
ylabel('Temperature \itT\rm / K','fontname','Arial','fontsize',fontsize); 
title('Solid: OilMixProp; Dashed: REFPROP')

if Lplot
    Material_mix = [Refrigerant{1},'_',Refrigerant{2},'_',num2str(pres_kPa/1e3),' MPa'];
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition',[0 0 18 8]);
    print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']); 
end