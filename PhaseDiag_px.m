%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  linewidth = 1; fontsize = 10;  markersize = 4;   
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);
AllEOS = {'PR','SRK','PTV','YFR'};      figure(1);clf; hold on; box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define fluids to study
% Refrigerant = {'R134A','Emkarate RL32'};  temp = 273.15 + 10;  
% Refrigerant = {'CO2','methane'};  temp = 273.15 - 50;  
% Refrigerant = {'methane','CO2'};  temp = 273.15 - 50;  
% Refrigerant = {'Nitrogen','Emkarate RL32'};  temp = 273.15 + 200;  
% Refrigerant = {'CO2','RENISO ACC HV'};  temp = 400.15;  
% Refrigerant = {'RENISO ACC HV','CO2'};  temp = 400.15;  
% Refrigerant = {'propane','R134A'};     temp = 283.15;  
Refrigerant = {'propane','C12'};     temp = 419.15; %   temp = 457.65;  
% Refrigerant = {'CO2','ethane'};  temp = 273.15 - 20;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define other parameters
Lplot = 1;                 % save the figure? 1 yes, 0 not
CubicEOS = AllEOS{4};      % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Insert your exp data here, only if they exist. 
%%%  x1    y1    P/MPa	
if strcmpi(Refrigerant{1},'propane') && strcmpi(Refrigerant{2},'C12') &&  temp == 419.15
    xyp = [
    0.0755	0.94865	0.453
    0.1048	0.95308	0.622
    0.1642	0.96513	1.011
    0.2766	0.97601	1.694
    0.3403	0.97512	2.185
    0.4651	0.97518	2.908
    0.5334	0.97598	3.397
    0.5938	0.97276	3.907
    0.7033	0.97255	4.82
    0.8322	0.95778	5.975
    0.8609	0.95338	6.255
    0.8774	0.93654	6.38
    ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main program - Nothing needs to be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare global parameters
GL1 = GetGlobals(CubicEOS,{Refrigerant{1}});  % obtain fluid constants
GL2 = GetGlobals(CubicEOS,{Refrigerant{2}});  % obtain fluid constants
GL = GetGlobals(CubicEOS,Refrigerant);  % obtain fluid constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define the pressure range 
try
    ff = OilPropm('all','T',temp,'Q',1,1,GL1,0,0);  psat1 = ff.p_Pa;  lve1 = 1;
catch
    psat1 = GL1.presc;    lve1 = 0;
end
try
    ff = OilPropm('all','T',temp,'Q',1,1,GL2,0,0);  psat2 = ff.p_Pa;   lve2 = 1;
catch
    psat2 = GL2.presc;    lve2 = 0;
end
p_s = min([psat1,psat2]); p_e = max([psat1,psat2]); 
if p_e / p_s > 20
    dp = (p_e-p_s)/1000; 
    pline = [p_s * exp(linspace(0,log((p_s+dp)/p_s),10)),  p_s+2*dp: 2*dp :p_s+9*dp,  p_s+10*dp: 20*dp :p_s+90*dp,  p_s+100*dp:200*dp:p_e];
else
    pline = linspace(p_s,p_e,20);
end
np = length(pline);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialization
Tratio = max(GL2.tempc/GL1.tempc,GL1.tempc/GL2.tempc);
if  lve1 == 1  &&  lve2 == 1  
    if psat2 > psat1,  x1 = ones(1,np);     y1 = ones(1,np);   x1(np) = 0; y1(np) = 0; end
    if psat2 < psat1,  x1 = zeros(1,np);     y1 = zeros(1,np);   x1(np) = 1; y1(np) = 1; end
    deleteindex = []; 
    if Tratio > 1.5 algorithm = 2; else algorithm = 1;end
elseif  lve1 == 0  && lve2 == 1  
    pline = [pline,p_e*linspace(1.1,Tratio*2,20)]; np = length(pline);
    if psat2 > psat1,  x1 = ones(1,np);     y1 = ones(1,np);   end
    if psat2 < psat1,  x1 = zeros(1,np);     y1 = zeros(1,np);    end
    deleteindex = np;
    algorithm = 2;
elseif  lve1 == 1  && lve2 == 0  
    pline = [pline,p_e*linspace(1.1,Tratio*2,20)];  np = length(pline); 
    if psat2 > psat1,  x1 = ones(1,np);     y1 = ones(1,np);   end
    if psat2 < psat1,  x1 = zeros(1,np);     y1 = zeros(1,np);  end
    deleteindex = np;
    algorithm = 2;
else
    error('Input temperature two high.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if refpropm is avaiable and can be calculated
nx = 60; deleteindex_ref = [];
x1_ref = linspace(0,1,nx);       y1_ref = linspace(0,1,nx); 
px1_ref = zeros(1,nx);           py1_ref = zeros(1,nx);   
for ix = 1:nx
    MF1 = x1_ref(ix);
    Zi = [MF1,1 - MF1]';
    [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,Zi);  
    try
        py1_ref(ix) = refpropm('P','T',temp,'Q',1,Refrigerant{1},Refrigerant{2},MassFrac)*1000;
        px1_ref(ix) = refpropm('P','T',temp,'Q',0,Refrigerant{1},Refrigerant{2},MassFrac)*1000;
        if py1_ref(ix) > px1_ref(ix) || py1_ref(ix) < 0  ||  px1_ref(ix) < 0  
            error(' '); 
        end 
    catch
        deleteindex_ref = [deleteindex_ref,ix];
    end
end
py1_ref(deleteindex_ref) = [];
px1_ref(deleteindex_ref) = [];
x1_ref(deleteindex_ref) = [];
y1_ref(deleteindex_ref) = [];
if length(x1_ref) > 2
    plot(x1_ref,px1_ref/1e6,'b--',y1_ref,py1_ref/1e6,'r--')
else
    fprintf('\n  =-- Refprop Calculation fail --=  \n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start the main loop
if algorithm == 1   % mainly for azeotropic mixture
    nx = 40;
    x1 = linspace(0,1,nx);       y1 = linspace(0,1,nx); 
    px1 = zeros(1,nx);           py1 = zeros(1,nx);   
    px1(1) = psat2;              py1(1) = psat2; 
    px1(nx) = psat1;             py1(nx) = psat1; 
    for ix = 2:nx-1
        MF1 = x1(ix);
        Zi = [MF1,1 - MF1]';
        if ix == 2
            p_y = py1(ix-1);   p_x = px1(ix-1);
        else
            p_y = ffy.p_Pa;    p_x = ffx.p_Pa;
        end
        try
            ffy = OilPropm('all','T',temp,'Q',1,Zi,GL,0,p_y/1000);
            py1(ix) = ffy.p_Pa;
            ffx = OilPropm('all','T',temp,'Q',0,Zi,GL,0,p_x/1000);
            px1(ix) = ffx.p_Pa;
            if py1(ix) > px1(ix) || py1(ix) < 0  ||  px1(ix) < 0  
                error(' '); 
            end 
        catch
            deleteindex = [deleteindex,ix];
        end
    end
    py1(deleteindex) = [];
    px1(deleteindex) = [];
    x1(deleteindex) = [];
    y1(deleteindex) = [];
    plot(x1,px1/1e6,'b-',y1,py1/1e6,'r-')
elseif algorithm == 2
    dx = 0.01;  ddx0 = 0.003;     
    firstPF = 0;
    MaxPFound = 0;
    for ip = 2:np-1
        if MaxPFound == 1, break; end
        if ip == 2 || firstPF == 0
            if psat2 > psat1  
                MF1 = 1 - dx;
            else
                MF1 = dx;
            end 
        else
            MF1 = (ffo.MoleF_Li(1) + ffo.MoleF_Vi(1))/2  + 2 * ((psat1 > psat2) -0.5) * dx;
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
            ff = OilPropm('all','T',temp,'P',pline(ip)/1000,Zi,GL,0,0);
            if strcmpi(ff.Phase,'V')
                MF1 = MF1 - 2 * ((psat1 > psat2) -0.5) * ddx;
            elseif strcmpi(ff.Phase,'L')
                MF1 = MF1 + 2 * ((psat1 > psat2) -0.5) * ddx;
            elseif strcmpi(ff.Phase,'LV')
                firstPF = 1;
                ffo = ff;
                break
            end
            if MF1 > 1 || MF1 < 0 || loopcount >= 20
                if pline(ip) < p_e
                    deleteindex = [deleteindex,ip];
                else
                    deleteindex = [deleteindex,ip:np]; MaxPFound = 1;
                end
                break;
            end
        end
        x1(ip) = ff.MoleF_Li(1);
        y1(ip) = ff.MoleF_Vi(1);
    end
    pline(deleteindex) = [];
    x1(deleteindex) = [];
    y1(deleteindex) = [];
    plot(x1,pline/1e6,'b-',y1,pline/1e6,'r-')
    if p_e / p_s > 20, set(gca, 'YScale', 'log'); end
end


if exist('xyp','var')
    plot(xyp(:,1),xyp(:,3),'+',"Color",'b')
    plot(xyp(:,2),xyp(:,3),'x',"Color",'r')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the results
xlabel(['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(temp) ,' K'],'fontname','Arial','fontsize',fontsize)
ylabel('Pressure \itp\rm / MPa','fontname','Arial','fontsize',fontsize); 
title('Solid: OilMixProp; Dashed: REFPROP')

if Lplot
    Material_mix = [Refrigerant{1},'_',Refrigerant{2},'_',num2str(temp),' K'];
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition',[0 0 18 8]);
    print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']); 
end