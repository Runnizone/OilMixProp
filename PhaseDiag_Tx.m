%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  linewidth = 1; fontsize = 10;  markersize = 4;   
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);
AllEOS = {'PR','SRK','PTV','YR'};    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define fluids to study
% Refrigerant = {'R134A','Emkarate RL32'};  pres = 3.5 * 1e6;  
% Refrigerant = {'Emkarate RL32','R134A'};  pres = 3.5 * 1e6;  
% Refrigerant = {'CO2','methane'};  pres =  7 * 1e6;  
% Refrigerant = {'methane','CO2'};  pres =  7 * 1e6;  
% Refrigerant = {'Nitrogen','Emkarate RL32'};  pres = 2 * 1e6;  
% Refrigerant = {'Emkarate RL32','Nitrogen'};  pres = 3.2 * 1e6;  
% Refrigerant = {'CO2','RENISO ACC HV'}; pres = 1 * 1e6;  
% Refrigerant = {'CO2','RENISO ACC HV'}; pres = 1 * 1e6;  
% Refrigerant = {'RENISO ACC HV','CO2'}; pres = 2 * 1e6;  
% Refrigerant = {'propane','R134A'};    pres = 1 * 1e6; 
Refrigerant = {'CO2','ethane'}; pres = 3 * 1e6;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define other parameters
Lplot = 1;                 % save the figure? 1 yes, 0 not
CubicEOS = AllEOS{3};      % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};  


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
    ff = OilPropm('P',pres,'Q',1,1,GL1,0,0);  Tsat1 = ff.T_K;  lve1 = 1;
catch
    Tsat1 = GL1.tempc;    lve1 = 0;
end
try
    ff = OilPropm('P',pres,'Q',1,1,GL2,0,0);  Tsat2 = ff.T_K;   lve2 = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start the main loop
figure(1);clf;
if algorithm == 1   % mainly for azeotropic mixture
    nx = 40;
    x1 = linspace(0,1,nx);       y1 = linspace(0,1,nx); 
    Tx1 = zeros(1,nx);           Ty1 = zeros(1,nx);   
    Tx1(1) = Tsat2;              Ty1(1) = Tsat2; 
    Tx1(nx) = Tsat1;             Ty1(nx) = Tsat1; 
    for ix = 2:nx-1
        MF1 = x1(ix);
        Zi = [MF1,1 - MF1]';
        if ix == 2
            T_y = Ty1(ix-1);   T_x = Tx1(ix-1);
        else
            T_y = ffy.T_K;    T_x = ffx.T_K;
        end
        try
            ffy = OilPropm('P',pres,'Q',1,Zi,GL,T_y,0);
            Ty1(ix) = ffy.T_K;
            ffx = OilPropm('P',pres,'Q',0,Zi,GL,T_x,0);
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
    plot(x1,Tx1,'-+',y1,Ty1,'-x')
elseif algorithm == 2
    dx = 0.01;  ddx0 = 0.003;     
    firstP = 0;
    for iT = 2:nT-1
        if iT == 2 || firstP == 0
            if Tsat2 > Tsat1  
                MF1 = 1 - dx;
            else
                MF1 = dx;
            end 
        else
            MF1 = (ffo.MoleF_Li(1) + ffo.MoleF_Vi(1))/2  + 2 * ((Tsat1 > Tsat2) -0.5) * dx;
        end
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
            ff = OilPropm('T',Tline(iT),'P',pres,Zi,GL,0,0);
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
    plot(x1,Tline,'-+',y1,Tline,'-x')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the results
xlabel(['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(pres/1e6) ,' MPa'],'fontname','Arial','fontsize',fontsize)
ylabel('Temperature \itT\rm / K','fontname','Arial','fontsize',fontsize); 

if Lplot
    Material_mix = [Refrigerant{1},'_',Refrigerant{2},'_',num2str(pres/1e6),' MPa'];
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition',[0 0 18 8]);
    print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']); 
end