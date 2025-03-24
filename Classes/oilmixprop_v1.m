
function oilmixprop_v1()

clear;close;clc;
%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Initializ   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
path(path,[pwd,'\Classes']);format short;   warning('off');
linewidth = 1; fontsize = 10;  markersize = 4;
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);

filename = '.\Classes\Fluid_Constants.txt';
fileID = fopen(filename);
C = textscan(fileID,'%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s','headerlines',1,'delimiter',';');
fclose(fileID);
fluidList=C{2};
filename = '.\Classes\Fluid_Constants_ext.txt';
fileID = fopen(filename);
C = textscan(fileID,'%s %s %s %s %f %f %f %f %f %f %f %f %f','headerlines',1,'delimiter',';');
fclose(fileID);
fluidList=[fluidList;C{2}];
filename = '.\Classes\Fluid_Constants_Fitted.txt';
fileID = fopen(filename);
C = textscan(fileID,'%s  %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1,'delimiter',';');
fclose(fileID);
fluidList=[fluidList;C{1}];


pdlist_g2={'Pressure vs. composition','Temperature vs. composition','Pressure vs. temperature'};
AllEOS = {'Peng-Robinson (PR)','Soave-Redlich-Kwong (SRK)','Patel-Teja-Valderrama (PTV)','Yang-Frotscher-Richter (YFR)'};
shortEOS = {'PR','SRK','PTV','YFR'};
flag_of_hand=1;

%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
fraction_g1=1;
vareos_g1=4;
CubicEOS_g1 = shortEOS{vareos_g1};
dataofg1={};
Refrigerant_g1={};
percent=0;
para1=[];para2=[];
para2list1 = {'P / kPa','D / kgm3','S / kgK','H / Jkg','Q / mass fraction'};
para2list2 = {'T / K','D / kgm3','S / kgK','H / Jkg','Q / mass fraction'};
allfraction={'MASS','MOLE'};
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
vareos_g2=1;
varpd_g2=1;
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
options = optimset('Display', 'off');  %  warning('off');
CubicEOS_g3 = shortEOS{4};
Lplot = 1;             % plot and save the figure? 1 yes, 0 not
Material_g3 = [];
kij_fit=0;
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
vareos_g4=4;
R = 8.31446261815324;     kB = 1.380649e-23;  % J / K  % m2 kg s-2 K-1
PowerConst = 2/7;   epsilon_kB_K = 300; sigma_nm = 0.51;
xi0 = 1.97E-10; Gamma = 5.42E-02; qDinv = 5.98E-10; xi_mu = 1;
Densityfile=[];
CPfile = [];
VISfile = [];
TCfile = [];
VPfile = [];
Material_g4 = [];
Zc = 0.272;       % give a good guess of the compressvility factor at the critical point
MM_gmol = 200;         % give a good guess of the molar mass
res2show={};


Fig = figure('Position',[100,100,1000,600],'menu','none',...
    'NumberTitle','off','Name','OilMixProp 1.0');

Fig.WindowButtonMotionFcn = @mouseMoved;
    function mouseMoved(~,~)
        mousePos = Fig.CurrentPoint;
        if (flag_of_hand && (mousePos(1) >= Fig.Position(3)*0.84) && (mousePos(2) >= Fig.Position(4)*0.12) && (mousePos(2) <= Fig.Position(4)*0.84)) || ...
                (flag_of_hand && (mousePos(1) >= Fig.Position(3)*0.32) && (mousePos(1) <= Fig.Position(3)*0.68) && (mousePos(2) >= Fig.Position(4)*0.85) )
            Fig.Pointer = 'hand';
        else
            Fig.Pointer = 'arrow';
        end
    end

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  pnl of    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
Pnl1_g0 = uipanel(Fig,'Position',[0,0,1,0.05]);
Pnl2_g0 = uipanel(Fig,'Position',[0.07,0.15,0.75,0.7]);
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
Pnl1_g1 = uipanel(Fig,'Position',[0.02,0.06,0.3,0.92],'Visible','off');
Pnl2_g1 = uipanel(Fig,'Position',[0.33,0.06,0.65,0.92],'Visible','off');
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
Pnl1_g2 = uipanel(Fig,'Position',[0.02,0.06,0.3,0.92],'Visible','off');
Pnl2_g2 = uipanel(Fig,'Position',[0.33,0.06,0.65,0.92],'Visible','off');
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
Pnl1_g3 = uipanel(Fig,'Position',[0.02,0.06,0.3,0.92],'Visible','off');
Pnl2_g3 = uipanel(Fig,'Position',[0.33,0.06,0.65,0.92],'Visible','off');
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
Pnl1_g4 = uipanel(Fig,'Position',[0.02,0.18,0.3,0.8],'Visible','off');
Pnl2_g4 = uipanel(Fig,'Position',[0.33,0.18,0.65,0.8],'Visible','off');
Pnl3_g4 = uipanel(Fig,'Position',[0.02,0.06,0.96,0.12],'Visible','off');

fontsize1=12;
fontsize2=14;
fontsize3=14;

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-TEXT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
% A=fileread('intro.txt');
intro_gui0={'Welcome to OilMixProp 1.0, a software package for important thermophysical properties of oils, common liquids and their mixtures.',...
    ' ',...
    'Key functions of OilMixProp 1.0:',...
    '1. Users can define their own oils (fuels, lubricants or other oil variants) using the fitting tools.',...
    '2. Inputs and outputs are specifically designed for the analysis of thermodynamic cycles.',...
    '3. So far 632 pure fluids with unlimited user-defined oils are available.',...
    '4. Key properties: density, phase behavior, heat capacity, entropy, enthalpy, speed of sound, viscosity, and thermal conductivity.',...
    ' ',...
    'Model within OilMixProp 1.0:',...
    ' ',...
    '1. Cubic EoS: Yang-Frotscher-Richter, Patel-Teja Valderrama, Soave-Redlich-Kwong and Peng-Robinson.',...
    '2. State-of-the-art residual entropy scaling (RES) for viscosity and thermal conductivity.',...
    '3. A robust flash algorithm to guarantee reliable and stable liquid-vapor equilibrium calculations.',...
    ' ',...
    'Developer: Dr. Xiaoxian Yang of TUCtt',...
    'Supported by: Prof. Dr. Markus Richter and Prof. Dr.-Ing. habil. Thorsten Urbaneck',...
    ' '};
% Text1_g0 = uicontrol(Pnl2_g0,'style','text','String',A,'Fontsize',12,...
%     'Units','normalized','Position',[0.01,0.01,0.98,0.98],'Visible','on',...
%     'HorizontalAlignment','left','FontName','Arial');

Text1_g0 = uicontrol(Pnl2_g0,'style','text','String',intro_gui0,'Fontsize',12,...
    'Units','normalized','Position',[0.01,0.01,0.98,0.98],'Visible','on',...
    'HorizontalAlignment','left','FontName','Arial');
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
Text1_g1 = uicontrol(Pnl1_g1,'style','text',...
    'String','Select EoS','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.83,0.8,0.05],'HorizontalAlignment','left');
Text2_g1 = uicontrol(Pnl1_g1,'style','text',...
    'String','Guess of Temperature / K:','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.31,0.66,0.05],'HorizontalAlignment','left');
Text3_g1 = uicontrol(Pnl1_g1,'style','text',...
    'String','Guess of Pressure / kPa:','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.19,0.63,0.05],'HorizontalAlignment','left');
Text4_g1 = uicontrol(Pnl2_g1,'style','edit',...
    'String','','Fontsize',12,'max',5,...
    'Units','normalized','Position',[0,0,1,1],'HorizontalAlignment','left','FontName','FixedWidth');

Text5_g1 = uicontrol(Pnl1_g1,'style','text',...
    'String','The 1st Input:','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.66,0.8,0.05],'HorizontalAlignment','left');
Text6_g1 = uicontrol(Pnl1_g1,'style','text',...
    'String','The 2nd Input:','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.5,0.8,0.05],'HorizontalAlignment','left');
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
Text1_g2 = uicontrol(Pnl1_g2,'style','text',...
    'String','Select Fluid 1','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.93,0.8,0.05],'HorizontalAlignment','left');
Text2_g2 = uicontrol(Pnl1_g2,'style','text',...
    'String','Select Fluid 2 ','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.76,0.8,0.05],'HorizontalAlignment','left');
Text3_g2 = uicontrol(Pnl1_g2,'style','text',...
    'String','Select phase diagram','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.59,0.8,0.05],'HorizontalAlignment','left');
Text4_g2 = uicontrol(Pnl1_g2,'style','text',...
    'String','Input Temperature / K:','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.20,0.8,0.05],'HorizontalAlignment','left');
Text5_g2 = uicontrol(Pnl1_g2,'style','text',...
    'String','Select EoS','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.42,0.8,0.05],'HorizontalAlignment','left');
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
Text1_g3 = uicontrol(Pnl1_g3,'style','text',...
    'String','Import Vapor Pressure Data','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.05,0.90,0.8,0.05],'HorizontalAlignment','left');
% Text2_g3 = uicontrol(Pnl1_g3,'style','text',...
%     'String','select fluid2','Fontsize',fontsize3,...
%     'Units','normalized','Position',[0.1,0.76,0.8,0.05]);
Text3_g3 = uicontrol(Pnl1_g3,'style','text',...
    'String','Select EoS','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.05,0.65,0.8,0.05],'HorizontalAlignment','left');

intro_gui4_1={'User Guide: Steps to Use'};
intro_gui4_2={'1. Create a Fluid Folder',...
    'In the MATLAB interface, navigate to the "ExpData" folder located on the left side.',...
    ' ',...
    'Create a new folder within "ExpData" for each measured fluid. Name the folder after the fluid, e.g., "PAG68". The format for recording data in the .txt files can be referenced from the "PAG68" folder.',...
    ' ',...
    'Within this folder, add the measurement data (VaporPressure) in separate .txt files.',...
    ' ',...
    '2. Click the "Load VP data" Button',...
    ' ',...
    '3. Select the "ExpData" Folder',...
    '',...
    '4. Add to database',...
    ''};

Text4_g3 = uicontrol(Pnl2_g3,'style','text','String',intro_gui4_1,'Fontsize',12,'FontWeight','bold',...
    'Units','normalized','Position',[0.01,0.91,0.98,0.07],'Visible','on',...
    'HorizontalAlignment','left','FontName','Arial');
Text5_g3 = uicontrol(Pnl2_g3,'style','text','String',intro_gui4_2,'Fontsize',12,...
    'Units','normalized','Position',[0.01,0.01,0.98,0.9],'Visible','on',...
    'HorizontalAlignment','left','FontName','Arial');

%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
Text1_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String','Fluid name','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.68,0.4,0.05],'HorizontalAlignment','left');
Text2_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String','Zc guess','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.6,0.4,0.05],'HorizontalAlignment','left');
Text3_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String','MM guess','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.52,0.4,0.05],'HorizontalAlignment','left');
Text4_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String','g/mol','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.83,0.52,0.2,0.05],'HorizontalAlignment','left');
Text5_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String','Import data','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.90,0.4,0.05],'HorizontalAlignment','left');
Text6_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String','Select EoS:','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.4,0.4,0.05],'HorizontalAlignment','left');

intro_gui3_1={'User Guide: Steps to Use'};
intro_gu3_2={'1. Create a Fluid Folder',...
    'In the MATLAB interface, navigate to the "ExpData" folder located on the left side.',...
    ' ',...
    'Create a new folder within "ExpData" for each measured fluid. Name the folder after the fluid, e.g., "PAG68". The format for recording data in the .txt files can be referenced from the "PAG68" folder.',...
    ' ',...
    'Within this folder, add the measurement data in separate .txt files with the following naming convention:',...
    '',...
    ' •  Density (required)'...
    ' •  Fluid_Constants_xxx (xxx is the fluid name)'...
    ' •  HeatCapacity'...
    ' •  ThermalConductivity'...
    ' •  VaporPressure'...
    ' •  Viscosity'...
    '',...
    '2. Click the "Density data" Button.',...
    '',...
    '3. Select the "ExpData" Folder',...
    'In the "ExpData" folder, find the fluid file that was just created and choose to open it.',...
    '',...
    '4. Add to database',...
    ''};

Text7_g4 = uicontrol(Pnl2_g4,'style','text','String',intro_gui3_1,'Fontsize',12,'FontWeight','bold',...
    'Units','normalized','Position',[0.01,0.91,0.98,0.07],'Visible','on',...
    'HorizontalAlignment','left','FontName','Arial');
Text8_g4 = uicontrol(Pnl2_g4,'style','text','String',intro_gu3_2,'Fontsize',12,...
    'Units','normalized','Position',[0.01,0.01,0.98,0.9],'Visible','on',...
    'HorizontalAlignment','left','FontName','Arial');

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-BUTTON    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
Bt1_g0= uicontrol(Pnl1_g0,'style','pushbutton','String','Property Calculation','Fontsize',14,...
    'Units','normalized','Position',[0,0,0.25,1],'Callback',@fluidcalculate_g0);
Bt2_g0 = uicontrol(Pnl1_g0,'style','pushbutton','String','Phase Diagram of Binaries','Fontsize',14,...
    'Units','normalized','Position',[0.25,0,0.25,1],'Callback',@drawpd_g0);
Bt3_g0= uicontrol(Pnl1_g0,'style','pushbutton','String','Fitter for Mixtures','Fontsize',14,...
    'Units','normalized','Position',[0.75,0,0.25,1],'Callback',@fitckij_g0);
Bt4_g0 = uicontrol(Pnl1_g0,'style','pushbutton','String','Fitter for Pure Oils','Fontsize',14,...
    'Units','normalized','Position',[0.5,0,0.25,1],'Callback',@fitpo_g0);
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
Bt1_g1 = uicontrol(Pnl1_g1,'style','pushbutton','String','Select Fluids','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.90,0.8,0.08],'Callback',@setfluid_g1);
% Bt2_g1= uicontrol(Pnl1_g1,'style','pushbutton','String','Del','Fontsize',fontsize1,...
%     'Units','normalized','Position',[0.3,0.41,0.1,0.2],'Callback',@delfromlist_g1);
% Bt3_g1= uicontrol(Pnl1_g1,'style','pushbutton','String','Next','Fontsize',fontsize1,...
%     'Units','normalized','Position',[0.3,0.2,0.1,0.2],'Callback',@list2table_g1);
% Bt4_g1= uicontrol(Pnl1_g1,'style','pushbutton','String','Done','Fontsize',fontsize1,...
%     'Units','normalized','Position',[0.746,0.08,0.2,0.16],'Callback',@datadone_g1,Visible='off');
Bt5_g1= uicontrol(Pnl1_g1,'style','pushbutton','String','Calculate','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.03,0.8,0.1],'Callback',@fluidcalculate_g1);
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
Bt1_g2= uicontrol(Pnl1_g2,'style','pushbutton','String','Draw','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.02,0.8,0.05],'Callback',@drawline_g2);
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
Bt1_g3= uicontrol(Pnl1_g3,'style','pushbutton','String','Fitting & Calculate','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.12,0.8,0.08],'Callback',@fitandcalc_g3);
Bt2_g3= uicontrol(Pnl1_g3,'style','pushbutton','String','Add to database','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.02,0.8,0.08],'Callback',@add2database_g3);
Bt3_g3= uicontrol(Pnl1_g3,'style','pushbutton','String','Load VP data','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.82,0.6,0.06],'Callback',@loadVPdata_g3);
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
Bt1_g4= uicontrol(Pnl1_g4,'style','pushbutton','String','Density data','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.05,0.82,0.43,0.06],'Callback',@loaddendata_g4);
% Bt2_g4= uicontrol(Pnl1_g4,'style','pushbutton','String','CP data','Fontsize',fontsize2,...
%     'Units','normalized','Position',[0.55,0.52,0.43,0.06],'Callback',@loadcpdata_g4);
% Bt3_g4= uicontrol(Pnl1_g4,'style','pushbutton','String','Viscosity data','Fontsize',fontsize2,...
%     'Units','normalized','Position',[0.05,0.42,0.43,0.06],'Callback',@loadvisdata_g4);
% Bt4_g4= uicontrol(Pnl1_g4,'style','pushbutton','String','TC data','Fontsize',fontsize2,...
%     'Units','normalized','Position',[0.55,0.42,0.43,0.06],'Callback',@loadtcdata_g4);
Bt5_g4= uicontrol(Pnl1_g4,'style','pushbutton','String','Fitting & Calculate','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.12,0.8,0.08],'Callback',@fitandcalc_g4);
Bt6_g4= uicontrol(Pnl1_g4,'style','pushbutton','String','Add to database','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.02,0.8,0.08],'Callback',@add2database_g4);

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-MENU    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
% Menu1_g1 = uicontrol(Pnl1_g1,'style','popupmenu',...
%     'String',{'Mass Fraction','Mole Fraction'},'Fontsize',fontsize1,...
%     'Units','normalized','Position',[0.74,0.88,0.2,0.05],'Visible','off');
Menu2_g1 = uicontrol(Pnl1_g1,'style','popupmenu',...
    'String',AllEOS,'Fontsize',fontsize2,'Value',4,...
    'Units','normalized','Position',[0.1,0.78,0.8,0.05],'Callback',@seteos_g1);
Menu3_g1 = uicontrol(Pnl1_g1,'style','popupmenu',...
    'String',{'T / K','P / kPa'},'Fontsize',fontsize1,...
    'Units','normalized','Position',[0.51,0.67,0.3,0.05],'Callback',@setpara1_g1);
Menu4_g1 = uicontrol(Pnl1_g1,'style','popupmenu',...
    'String',para2list1,'Fontsize',fontsize1,...
    'Units','normalized','Position',[0.51,0.51,0.3,0.05]);
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
Menu1_g2 = uicontrol(Pnl1_g2,'style','popupmenu',...
    'String',fluidList,'Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.88,0.8,0.05]);
Menu2_g2 = uicontrol(Pnl1_g2,'style','popupmenu',...
    'String',fluidList,'Fontsize',fontsize2,...
    'Units','normalized','Position',[0.1,0.71,0.8,0.05]);
Menu3_g2 = uicontrol(Pnl1_g2,'style','popupmenu',...
    'String',pdlist_g2,'Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.54,0.8,0.05],'Callback',@switchdata_g2);
Menu4_g2 = uicontrol(Pnl1_g2,'style','popupmenu',...
    'String',AllEOS,'Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.37,0.8,0.05],'Value',4);
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
% Menu1_g3 = uicontrol(Pnl1_g3,'style','popupmenu',...
%     'String',fluidList,'Fontsize',fontsize3,...
%     'Units','normalized','Position',[0.1,0.84,0.8,0.05]);
% Menu2_g3 = uicontrol(Pnl1_g3,'style','popupmenu',...
%     'String',fluidList,'Fontsize',fontsize3,...
%     'Units','normalized','Position',[0.1,0.68,0.8,0.05]);
Menu3_g3 = uicontrol(Pnl1_g3,'style','popupmenu',...
    'String',AllEOS,'Fontsize',fontsize3,'Value',4,...
    'Units','normalized','Position',[0.05,0.58,0.8,0.05]);
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
Menu1_g4 = uicontrol(Pnl1_g4,'style','popupmenu',...
    'String',AllEOS,'Fontsize',fontsize2,'Value',4,...
    'Units','normalized','Position',[0.05,0.34,0.8,0.05],'Callback',@seteos_g4);

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-EDIT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
%          注 GUI1有一个edit为Text4_g1
Edit1_g1 = uicontrol(Pnl1_g1,'style','edit',...
    'String','0','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.27,0.8,0.05]);
Edit2_g1 = uicontrol(Pnl1_g1,'style','edit',...
    'String','0','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.14,0.8,0.05]);
Edit3_g1 = uicontrol(Pnl1_g1,'style','edit',...
    'String',' ','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.62,0.8,0.05]);
Edit4_g1 = uicontrol(Pnl1_g1,'style','edit',...
    'String',' ','Fontsize',fontsize1,...
    'Units','normalized','Position',[0.1,0.46,0.8,0.05]);
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
Edit1_g2 = uicontrol(Pnl1_g2,'style','edit',...
    'String','','Fontsize',fontsize3,...
    'Units','normalized','Position',[0.1,0.13,0.8,0.05]);
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
% tempfluidname=['FN',sprintf('%03d',round(100*rand(1)))];
Edit1_g4 = uicontrol(Pnl1_g4,'style','text',...
    'String',[],'Fontsize',fontsize2,...
    'Units','normalized','Position',[0.42,0.68,0.4,0.05]);%注！此控件本为edit，后经讨论更正为text
Edit2_g4 = uicontrol(Pnl1_g4,'style','edit',...
    'String','0.272','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.39,0.6,0.4,0.05]);
Edit3_g4 = uicontrol(Pnl1_g4,'style','edit',...
    'String','200','Fontsize',fontsize2,...
    'Units','normalized','Position',[0.39,0.52,0.4,0.05]);

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-Listbox1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
% Listbox1_g1 = uicontrol(Pnl1_g1,'style','listbox',...
%     'String',fluidList,'Fontsize',10,...
%     'Unit','normalized','Position',[0.02,0.02,0.25,0.96]);
% Listbox2_g1 = uicontrol(Pnl1_g1,'style','listbox',...
%     'String',[],'Fontsize',10,...
%     'Unit','normalized','Position',[0.44,0.02,0.25,0.96]);
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-AXES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
Axes1_g0 = axes(Fig,'Position',[0.1,0.8,0.8,0.2]);
axis off;
[logo1, ~, alpha] = imread('GUI\logo\logo1.png');
f = imshow(logo1,'Parent',Axes1_g0);
set(f, 'AlphaData', alpha,'buttondownFcn',{@open_webpage1});

Axes2_g0 = axes(Fig,'Position',[0.82,0.72,0.15,0.12]);
axis off;
[logo1, ~, alpha] = imread('GUI\logo\logo2.png');
f = imshow(logo1,'Parent',Axes2_g0);
set(f, 'AlphaData', alpha,'buttondownFcn',{@open_webpage2});

Axes3_g0 = axes(Fig,'Position',[0.80,0.485,0.2,0.2]);
axis off;
[logo1, ~, alpha] = imread('GUI\logo\logo3.png');
f = imshow(logo1,'Parent',Axes3_g0);
set(f, 'AlphaData', alpha,'buttondownFcn',{@open_webpage3});

Axes4_g0 = axes(Fig,'Position',[0.825,0.32,0.15,0.15]);
axis off;
[logo1, ~, alpha]  = imread('GUI\logo\logo4.png');
f = imshow(logo1,'Parent',Axes4_g0);
set(f, 'AlphaData', alpha,'buttondownFcn',{@open_webpage4});

Axes5_g0 = axes(Fig,'Position',[0.81,0.12,0.18,0.16]);
axis off;
[logo1, ~, alpha]  = imread('GUI\logo\logo5.png');
f = imshow(logo1,'Parent',Axes5_g0);
set(f, 'AlphaData', alpha,'buttondownFcn',{@open_webpage5});

%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
Axes1_g1 = axes(Pnl1_g1,'Position',[0.77,0.32,0.05,0.05]);
axis off;
[logo1, ~, alpha]  = imread('GUI\logo\mark.png');
f = imshow(logo1,'Parent',Axes1_g1);
set(f, 'AlphaData', alpha,'buttondownFcn',{@tips1});
Axes2_g1 = axes(Pnl1_g1,'Position',[0.77,0.2,0.05,0.05]);
axis off;
f2 = imshow(logo1,'Parent',Axes2_g1);
set(f2, 'AlphaData', alpha,'buttondownFcn',{@tips2});
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
% Axes1_g2 = axes(Pnl2_g2,'Position',[0.1,0.1,0.9,0.9]);
% %     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
% Axes1_g3= axes(Pnl2_g3,'Position',[0.1,0.05,0.9,0.9]);
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
% Axes1_g4= axes(Pnl2_g4,'Position',[0.1,0.15,0.9,0.9]);
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
% Axes1_g3= axes(Pnl2_g3,'Position',[0.1,0.05,0.9,0.9]);

%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  controls-TABLE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colname_g4={'Materials','MM_gmol','Power_fix','Zc','Tc_K','Dc_kgm3','pc_MPa','acentric','k0_JmolK','k1_JmolK','n1_mu','n2_mu','n3_mu','n4_mu','n1_tc','n2_tc','n3_tc','n4_tc'};
%     %%%%%%%%%%%%%%%%%%  gui0    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui1    %%%%%%%%%%%%%%%%%%
Table1_g1 =  uitable(Pnl1_g1,'Data',{},'Position',[468,180,160,300],'Fontsize',10,'ColumnName',{'Materiac','percent'},'Visible','off','ColumnEditable',[false true],'Units','normalized');
%     %%%%%%%%%%%%%%%%%%  gui2    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui3    %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  gui4    %%%%%%%%%%%%%%%%%%
% Table1_g4 =  uitable(Pnl2_g4,'Data',{},'Position',[0,0,650,60],'ColumnName',colname_g4,'Fontsize',10);









%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Function    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%  funcgui0    %%%%%%%%%%%%%%%%%%
    function open_webpage1(~,~)
        web('https://github.com/runnizone/OilMixProp');
    end
    function open_webpage2(~,~)
        web('https://github.com/runnizone/OilMixProp');
    end
    function open_webpage3(~,~)
        web('https://www.tu-chemnitz.de/');
    end
    function open_webpage4(~,~)
        web('https://www.tu-chemnitz.de/mb/TechnThDyn/');
    end
    function open_webpage5(~,~)
        web('https://ketec.online/');
    end
    function tips1(~,~)
        A=fileread('GUI/tips.txt');
        warndlg(A, 'Tips');
    end
    function tips2(~,~)
        A=fileread('GUI/tips2.txt');
        warndlg(A, 'Tips');
    end

    function fluidcalculate_g0(~,~)
        set(Pnl2_g0,'Visible','off');
        set(Pnl1_g2,'Visible','off');
        set(Pnl2_g2,'Visible','off');
        set(Pnl1_g3,'Visible','off');
        set(Pnl2_g3,'Visible','off');
        set(Pnl1_g4,'Visible','off');
        set(Pnl2_g4,'Visible','off');
        set(Pnl3_g4,'Visible','off');
        set(Pnl1_g1,'Visible','on');
        set(Pnl2_g1,'Visible','on');
        set(Axes1_g0 ,'Visible','off');
        set(Axes2_g0 ,'Visible','off');
        set(Axes3_g0 ,'Visible','off');
        set(Axes4_g0 ,'Visible','off');
        set(Axes5_g0 ,'Visible','off');
        flag_of_hand=0;
    end
    function setfluid_g1(~,~)
        % set fluid selection page
        Fig_fluidselect = figure('Position',[150,150,800,500],'menu','none',...
            'NumberTitle','off','Name','DefineMixture');
        Pnltmp_g1=uipanel(Fig_fluidselect,'Position',[0.01,0.01,0.98,0.98]);
        Listboxtmp1_g1 = uicontrol(Pnltmp_g1,'style','listbox',...
            'String',fluidList,'Fontsize',10,...
            'Unit','normalized','Position',[0.01,0.01,0.4,0.92]);
        Listboxtmp2_g1 = uicontrol(Pnltmp_g1,'style','listbox',...
            'String',[],'Fontsize',10,...
            'Unit','normalized','Position',[0.59,0.01,0.4,0.92]);
        Texttmp1_g1 = uicontrol(Pnltmp_g1,'style','text',...
            'String','Available fluids','Fontsize',10,...
            'Units','normalized','Position',[0.01,0.95,0.3,0.03],'HorizontalAlignment','left');
        Texttmp2_g1 = uicontrol(Pnltmp_g1,'style','text',...
            'String','Selected mixture components','Fontsize',10,...
            'Units','normalized','Position',[0.59,0.95,0.3,0.03],'HorizontalAlignment','left');
        Bttmp1_g1 = uicontrol(Pnltmp_g1,'style','pushbutton','String','Add→','Fontsize',fontsize1,...
            'Units','normalized','Position',[0.43,0.62,0.12,0.06],'Callback',@add2list2_g1);
        Bttmp2_g1= uicontrol(Pnltmp_g1,'style','pushbutton','String','←Remove','Fontsize',fontsize1,...
            'Units','normalized','Position',[0.43,0.41,0.12,0.06],'Callback',@delfromlist_g1);
        Bttmp3_g1= uicontrol(Pnltmp_g1,'style','pushbutton','String','OK','Fontsize',fontsize1,...
            'Units','normalized','Position',[0.43,0.2,0.12,0.06],'Callback',@list2table_g1);
        function add2list2_g1(~,~)
            selected_index=get(Listboxtmp1_g1,'value');
            name_1=get(Listboxtmp2_g1,'String');
            [m,~]=size(name_1);
            if m==0
                name_1{1,1}=strjoin(fluidList(selected_index));
            else
                name_1{m+1,1}=strjoin(fluidList(selected_index));
            end
            % name_1=cellUnique(name_1);
            name_1=unique(name_1,'stable');
            Refrigerant_g1=name_1';
            set(Listboxtmp2_g1,'String',name_1);
        end
        function delfromlist_g1(~,~)
            selected_index=get(Listboxtmp2_g1,'value');
            name_1=get(Listboxtmp2_g1,'String');
            [m,~]=size(name_1);
            if m>0
                if  selected_index==m
                    set(Listboxtmp2_g1,'value',1);
                end
                rowsToKeep = setdiff(1:m, selected_index);
                name_1= name_1(rowsToKeep);
            end
            Refrigerant_g1=name_1';
            set(Listboxtmp2_g1,'String',name_1);
        end
        function list2table_g1(~,~)
            name_1=get(Listboxtmp2_g1,'String');
            if isempty(name_1)
                warndlg('Empty Chart', 'ERROR');
            else
                m=size(name_1,1);
                mixturename=[];
                if m>0
                    for i=1:m
                        name_1{i,2}=1/m;
                        mixturename=[mixturename,' & ',name_1{i,1}];
                    end
                end
                mixturename=mixturename(4:end);
                %创建比例界面
                Fig_specifymixturecomposition = figure('Position',[300,150,500,500],'menu','none',...
                    'NumberTitle','off','Name','SpecifyMixtureComposition');
                Texttmp2_1_g1 = uicontrol(Fig_specifymixturecomposition,'style','text',...
                    'String','Name','Fontsize',10,...
                    'Units','normalized','Position',[0.01,0.94,0.1,0.05],'HorizontalAlignment','left');
                Edittmp2_1_g1 = uicontrol(Fig_specifymixturecomposition,'style','edit',...
                    'String',mixturename,'Fontsize',10,...
                    'Units','normalized','Position',[0.11,0.94,0.78,0.05],'HorizontalAlignment','left');
                Pnltmp2_1_g1=uipanel(Fig_specifymixturecomposition,'Position',[0.01,0.2,0.98,0.72],'Title','Components');
                Menutmp2_1_g1 = uicontrol(Pnltmp2_1_g1,'style','popupmenu',...
                    'String',{'Mass Fraction','Mole Fraction'},'Fontsize',10,...
                    'Units','normalized','Position',[0.02,0.9,0.3,0.05],'Callback',@setfraction_g1);
                Tabletmp2_1_g1 =  uitable(Pnltmp2_1_g1,'Data',name_1,'Position',[20,20,460,240],...
                    'Fontsize',10,'ColumnName',{'Fluids','percent'},'ColumnEditable',[false true],'Units','normalized','ColumnWidth',{260,172},'CellEditCallback',@tablecalculate_g1);
                Bttmp2_1_g1= uicontrol(Fig_specifymixturecomposition,'style','pushbutton','String','OK','Fontsize',fontsize1,...
                    'Units','normalized','Position',[0.746,0.08,0.2,0.06],'Callback',@datadone_g1);
            end
            function setfraction_g1(~,~)
                fraction_g1=get(Menutmp2_1_g1,'Value');
            end

            function datadone_g1(~,~)
                dataofg1=get(Tabletmp2_1_g1,'Data');
                data=dataofg1(:,2);
                ptmp=1;
                for j = 1:length(data)-1
                    percent(j,1)=data{j};
                    ptmp=ptmp-data{j};
                end
                percent(length(data),1)=ptmp;
                if sum(percent)<0.98
                    warndlg('illegal Percent input', 'ERROR');
                else
                    %关闭两个fig
                    close SpecifyMixtureComposition;
                    close DefineMixture;
                end
                % 将流体名字及比例展示在Text4_g1
                if length(Refrigerant_g1)== 2
                    set(Text4_g1,'String',[num2str(percent(1)),' ',Refrigerant_g1{1},' + ',num2str(percent(2)),' ',Refrigerant_g1{2},' in ',allfraction{fraction_g1},' Fraction']);%
                elseif length(Refrigerant_g1)== 3
                    set(Text4_g1,'String',[num2str(percent(1)),' ',Refrigerant_g1{1},' + ',num2str(percent(2)),' ',Refrigerant_g1{2},' + ',num2str(percent(3)),' ',Refrigerant_g1{3},' in ',allfraction{fraction_g1},' Fraction']);%
                else
                    set(Text4_g1,'String',[Refrigerant_g1{1}]);%fprintf(Text4_g1,['01']);%Refrigerant{1},': ' ,CubicEOS_g1,' \n']);
                end

            end

            set(Table1_g1,'Data',name_1);


            function tablecalculate_g1(~,~)
                % 1 确保所有数据均在0-1之间
                m=size(name_1,1);
                nametmp=get(Tabletmp2_1_g1,'Data');
                legal=1;
                dataend=1;
                for i=1:m-1
                    dataend=dataend-nametmp{i,2};
                    legal=legal && nametmp{i,2}>0 && nametmp{i,2}<1;
                end
                if dataend>0 && dataend<1 && legal
                    nametmp{m,2}=dataend;
                    set(Tabletmp2_1_g1,'Data',nametmp);
                else
                    warndlg('illegal Percent input', 'ERROR');
                    set(Tabletmp2_1_g1,'Data',name_1);
                end
            end
        end
    end



    function drawpd_g0(~,~)
        set(Pnl2_g0,'Visible','off');
        set(Pnl1_g1,'Visible','off');
        set(Pnl2_g1,'Visible','off');
        set(Pnl1_g3,'Visible','off');
        set(Pnl2_g3,'Visible','off');
        set(Pnl1_g4,'Visible','off');
        set(Pnl2_g4,'Visible','off');
        set(Pnl3_g4,'Visible','off');
        set(Pnl1_g2,'Visible','on');
        set(Pnl2_g2,'Visible','on');
        set(Axes1_g0 ,'Visible','off');
        set(Axes2_g0 ,'Visible','off');
        set(Axes3_g0 ,'Visible','off');
        set(Axes4_g0 ,'Visible','off');
        set(Axes5_g0 ,'Visible','off');
        flag_of_hand=0;
        try
            close DefineMixture;
        catch
        end
    end

    function fitckij_g0(~,~)
        set(Pnl2_g0,'Visible','off');
        set(Pnl1_g2,'Visible','off');
        set(Pnl2_g2,'Visible','off');
        set(Pnl1_g1,'Visible','off');
        set(Pnl2_g1,'Visible','off');
        set(Pnl1_g4,'Visible','off');
        set(Pnl2_g4,'Visible','off');
        set(Pnl3_g4,'Visible','off');
        set(Pnl1_g3,'Visible','on');
        set(Pnl2_g3,'Visible','on');
        set(Axes1_g0 ,'Visible','off');
        set(Axes2_g0 ,'Visible','off');
        set(Axes3_g0 ,'Visible','off');
        set(Axes4_g0 ,'Visible','off');
        set(Axes5_g0 ,'Visible','off');
        flag_of_hand=0;
        try
            close DefineMixture;
        catch
        end
    end
    function fitpo_g0(~,~)
        set(Pnl2_g0,'Visible','off');
        set(Pnl1_g2,'Visible','off');
        set(Pnl2_g2,'Visible','off');
        set(Pnl1_g3,'Visible','off');
        set(Pnl2_g3,'Visible','off');
        set(Pnl1_g1,'Visible','off');
        set(Pnl2_g1,'Visible','off');
        set(Pnl1_g4,'Visible','on');
        set(Pnl2_g4,'Visible','on');
        set(Pnl3_g4,'Visible','on');
        set(Axes1_g0 ,'Visible','off');
        set(Axes2_g0 ,'Visible','off');
        set(Axes3_g0 ,'Visible','off');
        set(Axes4_g0 ,'Visible','off');
        set(Axes5_g0 ,'Visible','off');
        flag_of_hand=0;
        try
            close DefineMixture;
        catch
        end
    end
%     %%%%%%%%%%%%%%%%%%  funcgui1    %%%%%%%%%%%%%%%%%%
    function seteos_g1(~,~)
        vareos_g1=get(Menu2_g1,'Value');
        CubicEOS_g1 = shortEOS{vareos_g1};
    end
    function setpara1_g1(~,~)
        if get(Menu3_g1,'Value')==1
            set(Menu4_g1,'String',para2list1);
        else
            set(Menu4_g1,'String',para2list2);
        end
    end
% function add2list2_g1(~,~)
%     selected_index=get(Listboxtmp1_g1,'value');
%     name_1=get(Listboxtmp2_g1,'String');
%     [m,~]=size(name_1);
%     if m==0
%         name_1{1,1}=strjoin(fluidList(selected_index));
%     else
%         name_1{m+1,1}=strjoin(fluidList(selected_index));
%     end
%     name_1=cellUnique(name_1);
%     set(Listboxtmp2_g1,'String',name_1);
% end

% function setA_g1(~,~)
%     tmp=get(Menu3_g1,'Value');
%     if tmp>1
%         set(Text3_g1,'String','PSK');
%     else
%         set(Text3_g1,'String','K');
%     end
% end


% function delfromlist_g1(~,~)
%     selected_index=get(Listbox2_g1,'value');
%     name_1=get(Listbox2_g1,'String');
%     [m,~]=size(name_1);
%     if m>0
%         if  selected_index==m
%             set(Listbox2_g1,'value',1);
%         end
%         rowsToKeep = setdiff(1:m, selected_index);
%         name_1= name_1(rowsToKeep);
%     end
%     set(Listbox2_g1,'String',name_1);
% end

% function list2table_g1(~,~)
%     set(Table1_g1,'Visible','on');
%     set(Bt4_g1,'Visible','on');
%
%     set(Menu1_g1,'Visible','on');
%     name_1=get(Listbox2_g1,'String');
%     [m,~]=size(name_1);
%     for i=1:m
%         name_1{i,2}=0;
%     end
%     set(Table1_g1,'Data',name_1);
% end

% function datadone_g1(~,~)
%     data=get(Table1_g1,'Data');
%     data=data(:,2);
%     for i = 1:length(data)
%         percent(1,i)=data{i};
%     end
%     if sum(percent)~=1
%         warndlg('illegal Percent input', 'ERROR');
%     else
%         set(Listbox1_g1,'Visible','off');
%         set(Listbox2_g1,'Visible','off');
%         set(Bt1_g1,'Visible','off');
%         set(Bt2_g1,'Visible','off');
%         set(Bt3_g1,'Visible','off');
%
%         % 让1t中组件显示出来
%         set(Menu2_g1,'Visible','on');
%         set(Menu3_g1,'Visible','on');
%         set(Text1_g1,'Visible','on');
%         set(Text2_g1,'Visible','on');
%         set(Text3_g1,'Visible','on');
%         set(Edit1_g1,'Visible','on');
%         set(Bt5_g1,'Visible','on');
%     end
% end


    function fluidcalculate_g1(~,~)
        para1=str2num(get(Edit3_g1,'String'));
        para2=str2num(get(Edit4_g1,'String'));
        if length(para1) && length(para2) && length(Refrigerant_g1)
            T_K_guess=str2num(get(Edit1_g1,'String'));
            p_kPa_guess=str2num(get(Edit2_g1,'String'));
            para_of_g1 = get(Menu3_g1,'Value') * 10 + get(Menu4_g1,'Value');
            text2put=fluidcalculate_g1_out(CubicEOS_g1,Refrigerant_g1,fraction_g1,percent,para_of_g1,para1,para2,T_K_guess,p_kPa_guess,'OtherParameter');
            set (Text4_g1,'String',text2put);
        else
            warndlg('Finish Setting First', 'ERROR');
        end

    end
%     %%%%%%%%%%%%%%%%%%  funcgui2    %%%%%%%%%%%%%%%%%%
    function switchdata_g2(~,~)
        varpd_g2=get(Menu3_g2,'value');
        switch varpd_g2
            case 1
                set(Text4_g2,'String','Input Temperature / K:','fontsize',fontsize3);
            case 2
                set(Text4_g2,'String','Input Pressure / kPa:','fontsize',fontsize3);
            case 3
                set(Text4_g2,'String','Input Mass Frac. of Comp. 1:','fontsize',12);
        end

    end

    function drawline_g2(~,~)
        try
            delete(allchild(Pnl2_g2));
        end
        Axes1_g2 = axes(Pnl2_g2,'Position',[0.09,0.09,0.9,0.9]);

        varpd_g2=get(Menu3_g2,'value');
        var1=get(Menu1_g2,'value');
        fluid1=fluidList{var1};
        var2=get(Menu2_g2,'value');
        fluid2=fluidList{var2};
        Refrigerant={fluid1,fluid2};
        hold on; box on;
        algorithm = 2;




        %% 分三种情况进行绘图，暂未修改成function，待后续优化
        switch varpd_g2
            %% px
            case 1
                temp = str2num(get(Edit1_g2,'String'));
                %%% define other parameters
                Lplot = 1;                 % save the figure? 1 yes, 0 not
                vareos_g2=get(Menu4_g2,'Value');
                CubicEOS = shortEOS{vareos_g2};     % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};

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
                    pline = [p_s: 2*dp :p_s+9*dp,  p_s+10*dp: 20*dp :p_s+90*dp,  p_s+100*dp:200*dp:p_e];
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
                    warndlg('Input temperature too high.', 'ERROR');
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
                refcal=1;
                if length(x1_ref) > 2
                    plot(Axes1_g2,x1_ref,px1_ref/1e6,'b--',y1_ref,py1_ref/1e6,'r--');hold on;
                    % legend('RefProp','RefProp','Location','best');
                else
                    refcal=0;
                    warndlg(' =-- Refprop Calculation fail --= ', 'ERROR');
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
                    plot(Axes1_g2,x1,px1/1e6,'b-',y1,py1/1e6,'r-');hold on;

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
                    plot(Axes1_g2,x1,pline/1e6,'b-',y1,pline/1e6,'r-');hold on;

                    if p_e / p_s > 20, set(gca, 'YScale', 'log'); end
                end


                if exist('xyp','var')
                    plot(Axes1_g2,xyp(:,1),xyp(:,3),'+',"Color",'b');hold on;
                    plot(Axes1_g2,xyp(:,2),xyp(:,3),'x',"Color",'r');hold on;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% plot the results
                if refcal
                    legend(Axes1_g2,'RefProp','RefProp','OilMixProp','OilMixProp','Location','best');
                else
                    legend(Axes1_g2,'OilMixProp','OilMixProp','Location','best');
                end
                xlabel(Axes1_g2,['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(temp) ,' K'],'fontname','Arial','fontsize',fontsize)
                ylabel(Axes1_g2,'Pressure \itp\rm / MPa','fontname','Arial','fontsize',fontsize);
                hold off;


                %% tx
            case 2
                pres_kPa =str2num(get(Edit1_g2,'String'));
                Lplot = 1;                 % save the figure? 1 yes, 0 not
                vareos_g2=get(Menu4_g2,'Value');
                CubicEOS = shortEOS{vareos_g2};     % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Main program - Nothing needs to be changed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% prepare global parameters
                % CubicEOS = AllEOS{3};
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
                    warndlg('Input pressure too high.', 'ERROR');
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
                refcal=1;
                if length(x1_ref) > 2
                    plot(Axes1_g2,x1_ref,Tx1_ref,'b--',y1_ref,Ty1_ref,'r--');hold on;
                else
                    refcal=0;
                    warndlg('=-- Refprop Calculation fail --=', 'ERROR');
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
                    plot(Axes1_g2,x1,Tx1,'b-',y1,Ty1,'r-');hold on;

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
                    plot(Axes1_g2,x1,Tline,'b-',y1,Tline,'r-');hold on;
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% plot the results
                if refcal
                    legend(Axes1_g2,'RefProp','RefProp','OilMixProp','OilMixProp','Location','best');
                else
                    legend(Axes1_g2,'OilMixProp','OilMixProp','Location','best');
                end
                xlabel(Axes1_g2,['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(temp) ,' K'],'fontname','Arial','fontsize',fontsize)
                ylabel(Axes1_g2,'Pressure \itp\rm / MPa','fontname','Arial','fontsize',fontsize);
                hold off;

                % if Lplot
                %     Material_mix = [Refrigerant{1},'_',Refrigerant{2},'_',num2str(pres_kPa/1e3),' MPa'];
                %     set(gcf,'paperunits','centimeters');
                %     set(gcf,'paperposition',[0 0 18 8]);
                %     print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']);
                % end


                %% pt
            case 3
                MF1 =str2num(get(Edit1_g2,'String'));
                if MF1 <=0 || MF1 >=1
                    warndlg('illegal Percent input', 'ERROR');
                else
                    %%% define other parameters
                    Lplot = 1;                 % save the figure? 1 yes, 0 not
                    vareos_g2=get(Menu4_g2,'Value');
                    CubicEOS = shortEOS{vareos_g2};     % choose the cubic eos  AllEOS = {'PR','SRK','PTV','YR'};

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
                    refcal=1;
                    if length(Tx_ref) > 2
                        plot(Axes1_g2,Tx_ref,pline_ref/1e6,'r--',Ty_ref,pline_ref/1e6,'b--');hold on;
                    else
                        refcal=0;
                        warndlg(' =-- Refprop Calculation fail --= ', 'ERROR');
                    end


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% start the main loop
                    Tx = zeros(1,np);       Ty = zeros(1,np);
                    for ip = 1:np
                        if ip == 1
                            ff = OilPropm('all','P',pline(ip)/1000,'Q',1,1,GL1,0,0);
                            Tsat1 = ff.T_K;
                            ff = OilPropm('all','P',pline(ip)/1000,'Q',1,1,GL2,0,0);
                            Tsat2 = ff.T_K;
                            T_y = max([Tsat1,Tsat2]) + 20;
                            T_x = min([Tsat1,Tsat2]) - 20;
                        else
                            T_y = Ty(ip-1); % ffy.T_K;
                            T_x = Tx(ip-1); % ffx.T_K;
                        end
                        try
                            ffy = OilPropm('all','P',pline(ip)/1000,'Q',1,Zi,GL,T_y,0);
                            Ty(ip) = ffy.T_K;
                            ffx = OilPropm('all','P',pline(ip)/1000,'Q',0,Zi,GL,T_x,0);
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
                    plot(Axes1_g2,Tx,pline/1e6,'r-',Ty,pline/1e6,'b-');hold on;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% plot the results
                    % xlabel(['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(temp) ,' K'],'fontname','Arial','fontsize',fontsize)
                    if refcal
                        legend(Axes1_g2,'RefProp','RefProp','OilMixProp','OilMixProp','Location','best');
                    else
                        legend(Axes1_g2,'OilMixProp','OilMixProp','Location','best');
                    end
                    xlabel(Axes1_g2,['Mole Frac. of ',Refrigerant{1}, ' mixed with ',Refrigerant{2},' at ' ,num2str(temp) ,' K'],'fontname','Arial','fontsize',fontsize)
                    ylabel(Axes1_g2,'Pressure \itp\rm / MPa','fontname','Arial','fontsize',fontsize);
                    hold off;

                    % if Lplot
                    %     Material_mix = [Refrigerant{1},'_',Refrigerant{2},'_',num2str(MF1)];
                    %     set(gcf,'paperunits','centimeters');
                    %     set(gcf,'paperposition',[0 0 18 8]);
                    %     print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']);
                    % end
                end
        end
    end
%     %%%%%%%%%%%%%%%%%%  funcgui3    %%%%%%%%%%%%%%%%%%
    function add2database_g3(~,~)
        %% save the parameters to database

        copyfile .\Classes\Bin_kij_fit.txt '.\Classes\Bin_kij_fit_2.txt'
        %     try
        [Material_1,Material_2,kij_all] = textread('.\Classes\Bin_kij_fit_2.txt','%s%s%f','headerlines',1,'delimiter',';');
        nM = length(Material_1);
        fid_out = fopen(['.\Classes\Bin_kij_fit.txt'],'w');
        fprintf(fid_out,['          Material_1;          Material_2; Cubic_kij\n']);
        Lfound = 0;
        if nM == 0
            fprintf(fid_out,'%20s;%20s;%10.4f\n',Material_g3{1},Material_g3{2},kij_fit);
        end
        for im = 1:nM
            if (strcmpi(Material_g3{1},Material_1{im}) && strcmpi(Material_g3{2},Material_2{im})) ...
                    || (strcmpi(Material_g3{1},Material_2{im}) && strcmpi(Material_g3{2},Material_1{im}))
                Lfound = 1;
                fprintf(fid_out,'%20s;%20s;%10.4f\n',Material_g3{1},Material_g3{2},kij_fit);
            else
                fprintf(fid_out,'%20s;%20s;%10.4f\n',Material_1{im},Material_2{im},kij_all(im));
            end
            if im == nM && Lfound == 0
                fprintf(fid_out,'%20s;%20s;%10.4f\n',Material_g3{1},Material_g3{2},kij_fit);
            end
        end
        fclose(fid_out);
        delete '.\Classes\Bin_kij_fit_2.txt'
        warndlg('Done', 'Tips');%内容，标题
        %     catch
        %         copyfile 'Classes/Bin_kij_fit_2.txt' Classes/Bin_kij_fit.txt
        %         delete 'Classes/Bin_kij_fit_2.txt'
        %         error(['Error occurs'])
        %     end


    end


    function fitandcalc_g3(~,~)

        OilNum = 1;            % oil to be study, see the fluid definition below
        Lplot = 1;             % plot and save the figure? 1 yes, 0 not
        % Lsave2database = 1;    % save parameters to database (Classes/Bin_kij_fit.txt)? 1 yes, 0 not
        % var1=get(Menu1_g3,'value');
        % fluid1=fluidList{var1};
        % var2=get(Menu2_g3,'value');
        % fluid2=fluidList{var2};
        % Material_g3={fluid1,fluid2};
        var3=get(Menu3_g3,'value');
        CubicEOS = shortEOS{var3};



        % Oil mixture details
        if OilNum == 1
            Material_g3 = {'R1233zde','Emkarate RL32'};     % make sure oil is the second component
            FitIndex_VP = [1,9];                        % Index of data used for fit, only one point needed
        elseif OilNum == 2

        elseif OilNum == 3

        else
            error(['not defined yet'])
        end

        %%% parameter preperation
        GL = GetGlobals(CubicEOS,Material_g3);  % obtain fluid constants
        Material_mix = [Material_g3{1},'_',Material_g3{2}];
        VPfile = ['.\ExpData\',Material_mix,'\VaporPressure.txt'];
        LvpExist = 0;
        if exist(VPfile,'file')
            LvpExist = 1;
            [w100_oil_all,T_K_all,vp_Pa_all] = textread(VPfile,'%f%f%f','headerlines',1);
            massf_oil_all = w100_oil_all/100;
            molef_oil_all = massf_oil_all;
            for it = 1:length(massf_oil_all)
                massfrac = [1-massf_oil_all(it),massf_oil_all(it)]';
                [~,Zi] = EOSmodel.MassF_2_MoleF(GL.MM_gmol,massfrac);
                molef_oil_all(it) = Zi(2);
            end
            %% fit kij
            T_K = T_K_all(FitIndex_VP);
            vp_Pa = vp_Pa_all(FitIndex_VP);
            molef_oil = molef_oil_all(FitIndex_VP);
            xdata = [molef_oil, T_K, vp_Pa];   ydata = vp_Pa;
            n0 = [0];    lb = [-1];    ub = [1];
            n_fit = lsqcurvefit(@(n,xdata) EOSmodel.fit_CKij(n,xdata,GL),n0,xdata,ydata,lb,ub,options);
            vp_fit = EOSmodel.fit_CKij(n_fit,xdata,GL);
            kij_fit = n_fit(1);
            ReDev = (vp_Pa - vp_fit)./vp_fit;
            xdata = [molef_oil_all, T_K_all, vp_Pa_all];
            VP_fit_all = EOSmodel.fit_CKij(n_fit,xdata,GL);
            ReDev_all = (vp_Pa_all - VP_fit_all)./VP_fit_all;
        else
            error(['No VaporPressure.txt found in ExpData/',Material_mix,'/']);
        end

        %% plot
        if Lplot && LvpExist == 1
            try
                delete(allchild(Pnl2_g3));
            end
            Axes1_g3 = axes(Pnl2_g3,'Position',[0.09,0.09,0.9,0.9]);
            hold on; box on;
            plot(Axes1_g3,T_K_all,ReDev_all*100,'o','color','r','markersize',markersize,'linewidth',linewidth-0.4);
            plot(Axes1_g3,T_K,ReDev*100,'o','markerfacecolor','r','markeredgecolor','r','markersize',markersize,'linewidth',linewidth-0.4);
            plot(Axes1_g3,[min(T_K_all)-10,max(T_K_all)+10],[0,0],'k-');
            text(min(T_K_all)-5,3,['\itk\rm_i_j = ',num2str(kij_fit)])
            % axis([min(T_K_all),max(T_K_all)+10,-2.0,2.0])
            ylabel(['100\cdot(\itp\rm_b_u_b_,_e_x_p\it',char(hex2dec('2212')),'\itp\rm_b_u_b_,_f_i_t)/\itp\rm_b_u_b_,_f_i_t'],'rotation',90,'fontsize',fontsize,'FontName','Arial','linewidth',linewidth);
            % set(gcf,'paperunits','centimeters');
            % set(gcf,'paperposition',[0 0 18 8]);
            % print(gcf,'-dtiff','-r200',['.\Figures\',thisfilename,'_',Material_mix,'.tiff']);
            hold on; box on;
        end
    end

    function loadVPdata_g3(~,~)
        [file,path2] = uigetfile('VaporPressure.txt','Please select the VP file');
        VPfile = fullfile(path2,file);
    end
%     %%%%%%%%%%%%%%%%%%  funcgui4    %%%%%%%%%%%%%%%%%%
    function loaddendata_g4(~,~)
        [file,path2] = uigetfile('Density.txt','Please select the density file');
        Densityfile = fullfile(path2,file);
        CPfile = fullfile(path2,'HeatCapacity.txt');
        VISfile = fullfile(path2,'Viscosity.txt');
        TCfile = fullfile(path2,'ThermalConductivity.txt');

        % 给material设置名字
        Material=[];
        if ~isempty(strfind(Densityfile,'PAG68'))
            Material = 'PAG68';
        elseif ~isempty(strfind(Densityfile,'Emkarate RL32'))
            Material = 'Emkarate RL32';
        elseif ~isempty(strfind(Densityfile,'DIDP'))
            Material = 'DIDP';
        elseif ~isempty(strfind(Densityfile,'ISO VG32'))
            Material = 'ISO VG32';
        elseif ~isempty(strfind(Densityfile,'LAB ISO5'))
            Material = 'LAB ISO5';
        elseif ~isempty(strfind(Densityfile,'PEC5'))
            Material = 'PEC5';
        elseif ~isempty(strfind(Densityfile,'PEC7'))
            Material = 'PEC7';
        elseif ~isempty(strfind(Densityfile,'PEC9'))
            Material = 'PEC9';
        elseif ~isempty(strfind(Densityfile,'POE ISO68'))
            Material = 'POE ISO68';
        elseif ~isempty(strfind(Densityfile,'PEB8'))
            Material = 'PEB8';
        elseif ~isempty(strfind(Densityfile,'RENISO ACC HV'))
            Material = 'RENISO ACC HV';
        else
            warndlg('not defined yet', 'ERROR');
        end
        set(Edit1_g4,'String',Material)

    end
% function loadcpdata_g4(~,~)
%     [file,path3] = uigetfile('HeatCapacity.txt','Please select the CP file');
%     CPfile = fullfile(path3,file);
% end
% function loadvisdata_g4(~,~)
%     [file,path4] = uigetfile('Viscosity.txt','Please select the Viscosity file');
%     VISfile = fullfile(path4,file);
% end
% function loadtcdata_g4(~,~)
%     [file,path5] = uigetfile('ThermalConductivity.txt','Please select the TC file');
%     TCfile = fullfile(path5,file);
% end
    function seteos_g4(~,~)
        vareos_g4=get(Menu1_g4,'Value');
    end
    function add2database_g4(~,~)
        %% save the parameters to database
        Material_g4=get(Edit1_g4,'String');
        copyfile .\Classes\Fluid_Constants_Fitted.txt '.\Classes\Fluid_Constants_Fitted_2.txt'
        try
            [TheMaterial,MM_gmol,Power_fix,Zc,Tc_fit,Dc_fit,pc_fit,af_fit,k0_fit,k1_fit,...
                n1_mu_fit,n2_mu_fit,n3_mu_fit,n4_mu_fit,n1_tc_fit,n2_tc_fit,n3_tc_fit,n4_tc_fit] = ...
                textread(['.\ExpData\',Material_g4,'\Fluid_Constants_',Material_g4,'.txt'],'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'delimiter',';');
            [Material_all,MM_gmol_all,Power_fix_all,Zc_all,Tc_fit_all,Dc_fit_all,pc_fit_all,af_fit_all,k0_fit_all,k1_fit_all,...
                n1_mu_fit_all,n2_mu_fit_all,n3_mu_fit_all,n4_mu_fit_all,n1_tc_fit_all,n2_tc_fit_all,n3_tc_fit_all,n4_tc_fit_all] = ...
                textread('.\Classes\Fluid_Constants_Fitted_2.txt','%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'delimiter',';');
            nM = length(Material_all);
            fid_out = fopen(['.\Classes\Fluid_Constants_Fitted.txt'],'w');
            fprintf(fid_out,['     Materials;    MM_gmol;  Power_fix;      Zc;      Tc_K;   Dc_kgm3;    pc_MPa;   acentric;  k0_JmolK;  k1_JmolK;    n1_mu;     n2_mu;     n3_mu;    n4_mu;        n1_tc;       n2_tc;       n3_tc;       n4_tc\n']);
            Lfound = 0;
            if nM == 0
                fprintf(fid_out,'%14s;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.6f;%10.6f;%10.6f;%10.6f;%12.6f;%12.6f;%12.6f;%12.6f\n',...
                    TheMaterial{:},MM_gmol,Power_fix,Zc,Tc_fit,Dc_fit,pc_fit,af_fit,k0_fit,k1_fit,n1_mu_fit,n2_mu_fit,n3_mu_fit,n4_mu_fit,n1_tc_fit,n2_tc_fit,n3_tc_fit,n4_tc_fit);
            end
            for im = 1:nM
                if strcmpi(TheMaterial,Material_all{im})
                    Lfound = 1;
                    fprintf(fid_out,'%14s;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.6f;%10.6f;%10.6f;%10.6f;%12.6f;%12.6f;%12.6f;%12.6f\n',...
                        TheMaterial{:},MM_gmol,Power_fix,Zc,Tc_fit,Dc_fit,pc_fit,af_fit,k0_fit,k1_fit,n1_mu_fit,n2_mu_fit,n3_mu_fit,n4_mu_fit,n1_tc_fit,n2_tc_fit,n3_tc_fit,n4_tc_fit);
                else
                    fprintf(fid_out,'%14s;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.6f;%10.6f;%10.6f;%10.6f;%12.6f;%12.6f;%12.6f;%12.6f\n',...
                        Material_all{im},MM_gmol_all(im),Power_fix_all(im),Zc_all(im),Tc_fit_all(im),Dc_fit_all(im),pc_fit_all(im),af_fit_all(im),...
                        k0_fit_all(im),k1_fit_all(im),n1_mu_fit_all(im),n2_mu_fit_all(im),n3_mu_fit_all(im),n4_mu_fit_all(im),n1_tc_fit_all(im),n2_tc_fit_all(im),n3_tc_fit_all(im),n4_tc_fit_all(im));
                end
                if im == nM && Lfound == 0
                    fprintf(fid_out,'%14s;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.6f;%10.6f;%10.6f;%10.6f;%12.6f;%12.6f;%12.6f;%12.6f\n',...
                        TheMaterial{:},MM_gmol,Power_fix,Zc,Tc_fit,Dc_fit,pc_fit,af_fit,k0_fit,k1_fit,n1_mu_fit,n2_mu_fit,n3_mu_fit,n4_mu_fit,n1_tc_fit,n2_tc_fit,n3_tc_fit,n4_tc_fit);
                end
            end
            fclose(fid_out);
            delete '.\Classes\Fluid_Constants_Fitted_2.txt'
            warndlg('Done', 'Tips');%内容，标题
        catch
            copyfile '.\Classes\Fluid_Constants_Fitted_2.txt' .\Classes\Fluid_Constants_Fitted.txt
            delete '.\Classes\Fluid_Constants_Fitted_2.txt'
            error(['Error occurs, check  .\Classes\Fluid_Constants_',Material_g4,'.txt'])
        end
    end


    function fitandcalc_g4(~,~)
        try
            delete(allchild(Axes1_g4));
        end
        if exist(Densityfile,'file')==0
            warndlg('Import Density data first.', 'ERROR');
        else
            %% 此处暂时直接使用m文件

            % 20240507版本 根据0412版本的Fitter_PureOil.m文件调整
            %% Define the fluid to be fitted
            CubicEOS = shortEOS{vareos_g4};
            % OilNum = 1;            % oil to be study, see the fluid definition below
            Zc=str2num(get(Edit2_g4,'String'));
            MM_gmol = str2num(get(Edit3_g4,'String'));
            % Lplot = 1;             % plot and save the figure? 1 yes, 0 not
            % Lsave2database = 1;    % save parameters to database (Classes/Fluid_Constants_Fitted.txt)? 1 yes, 0 not

            % define the fitting parameters of the fluid
            if ~isempty(strfind(Densityfile,'PAG68'))
                Material = 'PAG68';
                FitIndex_D = [1,6];            % Index of density data used for fit, at least two
                FitIndex_cp = [1,26];          % Index of isobaric heat capacity data used for fit, at least two
                FitIndex_V = [1,5,8,10];       % Index of viscosity data used for fit, at least four
                FitIndex_TC = [1,2,3,4,5,6];   % Index of thermal conductivity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'Emkarate RL32'))
                Material = 'Emkarate RL32';
                FitIndex_D = [1,3,4];            % Index of density data used for fit, at least two
                FitIndex_V = [1,2,3,4];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'DIDP'))
                Material = 'DIDP';
                MM_gmol = 446.672;          % give a good guess of the molar mass
                FitIndex_D = [1,7,15];            % Index of density data used for fit, at least two
                FitIndex_cp = [1,3,7];          % Index of isobaric heat capacity data used for fit, at least two
                FitIndex_V = [1,2,3,4];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'ISO VG32'))
                Material = 'ISO VG32';
                MM_gmol = 560.26;          % give a good guess of the molar mass
                FitIndex_D = [1,3,5];            % Index of density data used for fit, at least two
                FitIndex_V = [1,3,5];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'LAB ISO5'))
                Material = 'LAB ISO5';
                MM_gmol = 240;          % give a good guess of the molar mass
                FitIndex_D = [1,5,9];            % Index of density data used for fit, at least two
                FitIndex_V = [1,4,7,9];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'PEC5'))
                Material = 'PEC5';
                MM_gmol = 472.6120;          % give a good guess of the molar mass
                Zc = 0.2598;
                FitIndex_D = [30,75,120,165];            % Index of density data used for fit, at least two
                FitIndex_V = [12, 42,103,118,132,147,169,199,229,254];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'PEC7'))
                Material = 'PEC7';
                MM_gmol = 584.8350;          % give a good guess of the molar mass
                Zc = 0.2590;
                FitIndex_D = [26,71,116,161];            % Index of density data used for fit, at least two
                FitIndex_V = [9,39,56,120,136,162,178,187 194,210,225,247];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'PEC9'))
                Material = 'PEC9';
                MM_gmol = 697.0510;          % give a good guess of the molar mass
                Zc = 0.2560;
                FitIndex_D = [10,55,100,145];         % Index of density data used for fit, at least two
                FitIndex_V = [1,31,44,60,92];       % Index of viscosity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'POE ISO68'))
                Material = 'POE ISO68';
                MM_gmol = 500;          % give a good guess of the molar mass
                FitIndex_D = [1,9,17];         % Index of density data used for fit, at least two
                FitIndex_cp = [1,9,17];          % Index of isobaric heat capacity data used for fit, at least two
                FitIndex_TC = [1,8,9,17];   % Index of thermal conductivity data used for fit, at least four
            elseif ~isempty(strfind(Densityfile,'PEB8'))
                Material = 'PEB8';
                Zc = 0.2600;
                MM_gmol = 640.9;          % give a good guess of the molar mass
                FitIndex_D = [1,34,67,89];         % Index of density data used for fit, at least two
                FitIndex_V = [1, 15, 29,43, 57,71];          % Index of isobaric heat capacity data used for fit, at least two
            elseif ~isempty(strfind(Densityfile,'RENISO ACC HV'))
                Material = 'RENISO ACC HV';
                MM_gmol = 200;          % give a good guess of the molar mass
                FitIndex_D = [1,9,19];         % Index of density data used for fit, at least two
            else
                warndlg('not defined yet', 'ERROR');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~exist('Zc','var')
                if strcmpi(CubicEOS,'PTV'), Zc = 0.2563; elseif strcmpi(CubicEOS,'YR'), Zc = 0.2640; else, end
            end
            if ~exist('MM_gmol','var')
                MM_gmol = 200;
            end

            %%  constants and prepare output file
            fid_out = fopen(['ExpData/',Material,'/Fluid_Constants_',Material,'.txt'],'w');
            fprintf(fid_out,['     Materials;    MM_gmol;  Power_fix;      Zc;      Tc_K;   Dc_kgm3;    pc_MPa;   acentric;  k0_JmolK;  k1_JmolK;    n1_mu;     n2_mu;     n3_mu;    n4_mu;        n1_tc;       n2_tc;       n3_tc;       n4_tc\n']);
            LdenExist = 0; LcpExist = 0; LvisExist = 0; LtcExist = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% get 'experimetnal data'   --- density
            if exist(Densityfile,'file')
                LdenExist = 1;
                [T_D_K_all,p_D_Pa_all,D_kgm3_all] = textread(Densityfile,'%f%f%f','headerlines',1);
                T_D_K = T_D_K_all(FitIndex_D);
                p_D_Pa = p_D_Pa_all(FitIndex_D);
                den_kg = D_kgm3_all(FitIndex_D);
                den_mol = den_kg/MM_gmol*1000;   %  mol/m3
                %% fit Tc and Dc
                const = [Zc,PowerConst];
                xdata = T_D_K';   ydata = den_kg';
                n0 = [500, 600];       lb = [0,0];    ub = [1500,1000];
                n_fit = lsqcurvefit(@(n,xdata) EOSmodel.fit_Reckett_TcDc(n,xdata,const),n0,xdata,ydata,lb,ub,options);
                Tc_fit = n_fit(1);
                Dc_fit = n_fit(2);
                %% fit acentric factor
                pc_fit = Zc / (MM_gmol / R / Tc_fit/ Dc_fit / 1000);
                const = [Tc_fit, pc_fit, Zc,R];
                xdata = [T_D_K, p_D_Pa];   ydata = den_mol;
                n0 = [0.5];    lb = [-0.3];    ub = [1.5];
                n_fit = lsqcurvefit(@(n,xdata) EOSmodel.fit_CEOS_rho(n,xdata,const,CubicEOS),n0,xdata,ydata,lb,ub,options);
                den_mol_fit = EOSmodel.fit_CEOS_rho(n_fit,xdata,const,CubicEOS);
                af_fit = n_fit(1);
                den_kg_fit = den_mol_fit*MM_gmol/1000;
                ReDev_D = (den_kg - den_kg_fit)./den_kg_fit;
                xdata = [T_D_K_all, p_D_Pa_all];
                den_kg_fit_all = EOSmodel.fit_CEOS_rho(n_fit,xdata,const,CubicEOS)*MM_gmol/1000;
                ReDev_D_all = (D_kgm3_all - den_kg_fit_all)./den_kg_fit_all;
            else
                error(['No Density.txt found in ExpData/',Material,'/']);
            end

            N_A = R / kB;
            epsilon_kB_K = Tc_fit/1.2593;
            sigma_nm = (0.3189 / Dc_fit*MM_gmol*1e24/N_A)^(1/3);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% get 'experimetnal data'  --- cp
            if exist(CPfile,'file')
                LcpExist = 1;
                [T_cp_K_all,p_cp_Pa_all,cp_JgK_all] = textread(CPfile,'%f%f%f','headerlines',1);
                T_cp_K = T_cp_K_all(FitIndex_cp);
                p_cp_Pa = p_cp_Pa_all(FitIndex_cp);
                cp_JkgK = cp_JgK_all(FitIndex_cp)*1000;
                cp_JmolK = cp_JkgK*MM_gmol/1000;   %  [J/mol K)]
                %% fit cp
                const = [Tc_fit, pc_fit, Zc, af_fit, R];
                xdata = [T_cp_K, p_cp_Pa];   ydata = cp_JmolK;
                n0 = [200,200];
                lb = [10, 10];       ub = [1000, 3000];
                n_fit = lsqcurvefit(@(n,xdata) EOSmodel.fit_CEOS_cp(n,xdata,const,CubicEOS),n0,xdata,ydata,lb,ub,options);
                cp_mol_fit = EOSmodel.fit_CEOS_cp(n_fit,xdata,const,CubicEOS);
                cp_kg_fit = cp_mol_fit/MM_gmol*1000;
                ReDev_cp = (cp_JkgK - cp_kg_fit)./cp_kg_fit;
                k0_JmolK = n_fit(1);
                k1_JmolK = n_fit(2);
                xdata = [T_cp_K_all, p_cp_Pa_all];
                cp_kg_fit_all = EOSmodel.fit_CEOS_cp(n_fit,xdata,const,CubicEOS)/MM_gmol*1000;
                ReDev_cp_all = (cp_JgK_all*1000 - cp_kg_fit_all)./cp_kg_fit_all;
            else
                disp(['Warning: No HeatCapacity.txt found in ExpData/',Material,'/']);
                disp(['         Set k0_JmolK = 300, and k1_JmolK = 400']);
                k0_JmolK = 300;        k1_JmolK = 400;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% get 'experimetnal data'   --- viscosity
            if exist(VISfile,'file')
                LvisExist = 1;
                [T_V_K_all,p_V_Pa_all,V_mPas_all] = textread(VISfile,'%f%f%f','headerlines',1);
                T_V_K = T_V_K_all(FitIndex_V);
                p_V_Pa = p_V_Pa_all(FitIndex_V);
                V_mPas = V_mPas_all(FitIndex_V);
                %% fit viscosity n
                xdata = [T_V_K, p_V_Pa];   ydata = V_mPas;
                n0 = [0.220249  -0.070232    0.011963    0.000000];
                lb = [-10 -10 -10 0];       ub = [10 10 10 0];
                const = [Tc_fit, pc_fit, Zc, af_fit, MM_gmol,epsilon_kB_K,sigma_nm, kB,R];
                n_fit = lsqcurvefit(@(n,xdata) EOSmodel.fit_CEOS_vis_mPas_n(n,xdata,const,CubicEOS),n0,xdata,ydata,lb,ub,options);
                v_mPas_fit = EOSmodel.fit_CEOS_vis_mPas_n(n_fit,xdata,const,CubicEOS);
                ReDev_v = (V_mPas - v_mPas_fit)./v_mPas_fit;
                n_mu = n_fit;
                xdata = [T_V_K_all, p_V_Pa_all];
                v_mPas_fit_all = EOSmodel.fit_CEOS_vis_mPas_n(n_fit,xdata,const,CubicEOS);
                ReDev_v_all = (V_mPas_all - v_mPas_fit_all)./v_mPas_fit_all;
            else
                disp(['Warning: No Viscosity.txt found in ExpData/',Material,'/']);
                disp([' Set n_mu as defaulted values']);
                n_mu = [0.220249  -0.070232    0.011963    0.000000];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% get 'experimetnal data'   --- ThermalConductivity
            if exist(TCfile,'file')
                LtcExist = 1;
                [T_L_K_all,p_L_Pa_all,L_WmK_all] = textread(TCfile,'%f%f%f','headerlines',1);
                T_L_K = T_L_K_all(FitIndex_TC);
                p_L_Pa = p_L_Pa_all(FitIndex_TC);
                L_WmK = L_WmK_all(FitIndex_TC);
                Dc_mol = Dc_fit /MM_gmol * 1000;
                %% fit thermal conductivity  n
                xdata = [T_L_K, p_L_Pa];   ydata = L_WmK;
                n0 = [0.001000   0.895835    2.305079   -0.322008];
                lb = [0 -99 -99 -99];       ub = [99 99 99 99];
                %         n0 = [0.00100   4.835594   -0.365323    0.000000];
                %         lb = [0 -99 -99 0];       ub = [99 99 99 0];
                const = [Tc_fit, pc_fit, Dc_mol, Zc, af_fit, MM_gmol, k0_JmolK, k1_JmolK, xi_mu,epsilon_kB_K,sigma_nm, xi0, Gamma, qDinv, kB,R];
                n_fit = lsqcurvefit(@(n,xdata) EOSmodel.fit_CEOS_tc_n(n,xdata,const,n_mu,CubicEOS),n0,xdata,ydata,lb,ub,options);
                L_WmK_fit = EOSmodel.fit_CEOS_tc_n(n_fit,xdata,const,n_mu,CubicEOS);
                ReDev_L = (L_WmK - L_WmK_fit)./L_WmK_fit;
                n_tc = n_fit;
                xdata = [T_L_K_all, p_L_Pa_all];
                L_WmK_fit_all = EOSmodel.fit_CEOS_tc_n(n_fit,xdata,const,n_mu,CubicEOS);
                ReDev_tc_all = (L_WmK_all - L_WmK_fit_all)./L_WmK_fit_all;
            else
                disp(['Warning: No ThermalConductivity.txt found in ExpData/',Material,'/']);
                disp([' Set n_tc as defaulted values']);
                n_tc = [0.001000   0.895835    2.305079   -0.322008];
            end

            %% save parameters
            fprintf(fid_out,'%14s;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.6f;%10.6f;%10.6f;%10.6f;%12.6f;%12.6f;%12.6f;%12.6f\n',...
                Material,MM_gmol,PowerConst,Zc,Tc_fit,Dc_fit,pc_fit/1e6,af_fit,k0_JmolK,k1_JmolK,n_mu,n_tc);
            fclose(fid_out);
            % res2show={'Materials','MM_gmol','Power_fix','Zc','Tc_K','Dc_kgm3','pc_MPa','acentric','k0_JmolK','k1_JmolK','n1_mu','n2_mu','n3_mu','n4_mu','n1_tc','n2_tc','n3_tc','n4_tc'};
            res2show={Material,MM_gmol,PowerConst,Zc,Tc_fit,Dc_fit,pc_fit/1e6,af_fit,k0_JmolK,k1_JmolK,n_mu(1),n_mu(2),n_mu(3),n_mu(4),n_tc(1),n_tc(2),n_tc(3),n_tc(4)};


            %% plot the results
            if Lplot
                try
                    delete(allchild(Pnl2_g4));
                    delete(allchild(Pnl3_g4));
                end
                Axes1_g4 = axes(Pnl2_g4,'Position',[0.05,0.05,0.9,0.9]);
                Table1_g4 =  uitable(Pnl3_g4,'Data',{},'Position',[0,0,958,69],'ColumnName',colname_g4,'Fontsize',10,'Units','normalized');
                set(Table1_g4,'Data',res2show);
                subplot(221); hold on; box on;
                if LdenExist == 1
                    plot(T_D_K_all,ReDev_D_all*100,'o','color','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(T_D_K,ReDev_D*100,'o','markerfacecolor','r','markeredgecolor','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(get(gca,'xlim'),[0,0],'k-');
                    % axis([min(T_D_K_all),max(T_D_K_all)+10,-2.0,2.0])
                    ylabel(['100\cdot(\it\rho\rm_e_x_p\it',char(hex2dec('2212')),'\it\rho\rm_f_i_t)/\it\rho\rm_f_i_t'],'rotation',90,'fontsize',fontsize,'FontName','Arial','linewidth',linewidth);
                end

                if LcpExist == 1
                    subplot(222); hold on; box on;
                    plot(T_cp_K_all,ReDev_cp_all*100,'o','color','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(T_cp_K,ReDev_cp*100,'o','markerfacecolor','r','markeredgecolor','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(get(gca,'xlim'),[0,0],'k-');
                    % axis([min(T_cp_K)-10,max(T_cp_K)+10,-0.6,0.6])
                    ylabel(['100\cdot(\itc_p\rm_,_e_x_p\it',char(hex2dec('2212')),'\itc_p\rm_,_f_i_t)/\itc_p\rm_,_f_i_t'],'rotation',90,'fontsize',fontsize,'FontName','Arial','linewidth',linewidth);
                end

                if LvisExist == 1
                    subplot(223); hold on; box on;
                    plot(T_V_K_all,ReDev_v_all*100,'o','color','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(T_V_K,ReDev_v*100,'o','markerfacecolor','r','markeredgecolor','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(get(gca,'xlim'),[0,0],'k-');
                    % axis([min(T_v_C(DataRange_v_all))+273.15-10,max(T_v_C(DataRange_v_all))+273.15+10,-5,5])
                    ylabel(['100\cdot(\it\eta\rm_e_x_p\it',char(hex2dec('2212')),'\it\eta\rm_f_i_t)/\it\eta\rm_f_i_t'],'rotation',90,'fontsize',fontsize,'FontName','Arial','linewidth',linewidth);
                    xlabel('\itT\rm / K','FontSize',fontsize,FontName='Arial')
                end

                if LtcExist == 1
                    subplot(224); hold on; box on;
                    plot(T_L_K_all,ReDev_tc_all*100,'o','color','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(T_L_K,ReDev_L*100,'o','markerfacecolor','r','markeredgecolor','r','markersize',markersize,'linewidth',linewidth-0.4);
                    plot(get(gca,'xlim'),[0,0],'k-');
                    % axis([min(T_L)-10,max(T_L)+10,-1,1])
                    ylabel(['100\cdot(\it\lambda\rm_e_x_p\it',char(hex2dec('2212')),'\it\lambda\rm_f_i_t)/\it\lambda\rm_f_i_t'],'rotation',90,'fontsize',fontsize,'FontName','Arial','linewidth',linewidth);
                    xlabel('\itT\rm / K','FontSize',fontsize,FontName='Arial')
                end

                % set(gcf,'paperunits','centimeters');
                % set(gcf,'paperposition',[0 0 18 8]);
                % print(gcf,'-dtiff','-r200',[pwd,'\Figures\',thisfilename,'_',Material_g4,'.tiff']);
                hold off;
            end


        end

    end
end



