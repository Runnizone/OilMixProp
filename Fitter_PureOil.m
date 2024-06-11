%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Preparation
clear; clc; format short; path(path,'C:\Xiaoxian\Github\CubicEOS\Classes'); 
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);
AllEOS = {'PR','SRK','PTV','YFR'};    
linewidth = 1; fontsize = 10;  markersize = 4;
options = optimset('Display',  'off');   R = 8.31446261815324;     kB = 1.380649e-23;  % J / K  % m2 kg s-2 K-1
PowerConst = 2/7;  xi0 = 1.97E-10; Gamma = 5.42E-02; qDinv = 5.98E-10; xi_mu = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the fluid to be fitted
iOilNum = 2;            % oil to be study, see the fluid definition below, if 0, run all oils.
CubicEOS = AllEOS{4}; 

N_Oil = 11;            % the number of avaiable oil 
Lplot = 1;             % plot and save the figure? 1 yes, 0 not
Lsave2database = 1;    % save parameters to database (Classes/Fluid_Constants_Fitted.txt)? 1 yes, 0 not

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start the fitting
if iOilNum == 0, OilRun = 1:N_Oil; else, OilRun = iOilNum;  end
for OilNum = OilRun
    % define the fitting parameters of the fluid
    if OilNum == 1
        Material = 'PAG68'; 
        FitIndex_D = [1,6];            % Index of density data used for fit, at least two
        FitIndex_cp = [1,26];          % Index of isobaric heat capacity data used for fit, at least two
        FitIndex_V = [1,5,8,10];       % Index of viscosity data used for fit, at least four
        FitIndex_TC = [1,2,3,4,5,6];   % Index of thermal conductivity data used for fit, at least four
    elseif OilNum == 2
        Material = 'Emkarate RL32'; 
        FitIndex_D = [1,5,9];            % Index of density data used for fit, at least two
        FitIndex_V = [1,3,5,7];       % Index of viscosity data used for fit, at least four
        MM_gmol = 180;  
    elseif OilNum == 3
        Material = 'DIDP'; 
        MM_gmol = 446.672;          % give a good guess of the molar mass
        FitIndex_D = [1,7,15];            % Index of density data used for fit, at least two
        FitIndex_cp = [1,3,7];          % Index of isobaric heat capacity data used for fit, at least two
        FitIndex_V = [1,2,3,4];       % Index of viscosity data used for fit, at least four
    elseif OilNum == 4
        Material = 'ISO VG32'; 
        MM_gmol = 560.26;          % give a good guess of the molar mass
        FitIndex_D = [1,3,5];            % Index of density data used for fit, at least two
        FitIndex_V = [1,3,5];       % Index of viscosity data used for fit, at least four
    elseif OilNum == 5
        Material = 'LAB ISO5'; 
        MM_gmol = 240;          % give a good guess of the molar mass
        FitIndex_D = [1,5,9];            % Index of density data used for fit, at least two
        FitIndex_V = [1,4,7,9];       % Index of viscosity data used for fit, at least four
    elseif OilNum == 6
        Material = 'PEC5'; 
        MM_gmol = 472.6120;          % give a good guess of the molar mass
        Zc = 0.2598;
        FitIndex_D = [30,75,120,165];            % Index of density data used for fit, at least two
        FitIndex_V = [12, 42,103,118,132,147,169,199,229,254];       % Index of viscosity data used for fit, at least four
    elseif OilNum == 7
        Material = 'PEC7'; 
        MM_gmol = 584.8350;          % give a good guess of the molar mass
        Zc = 0.2590;
        FitIndex_D = [26,71,116,161];            % Index of density data used for fit, at least two
        FitIndex_V = [9,39,56,120,136,162,178,187 194,210,225,247];       % Index of viscosity data used for fit, at least four
    elseif OilNum == 8
        Material = 'PEC9'; 
        MM_gmol = 697.0510;          % give a good guess of the molar mass
        Zc = 0.2560;
        FitIndex_D = [10,55,100,145];         % Index of density data used for fit, at least two
        FitIndex_V = [1,31,44,60,92];       % Index of viscosity data used for fit, at least four
    elseif OilNum == 9
        Material = 'HATCOL 4467'; 
        MM_gmol = 500;          % give a good guess of the molar mass
        FitIndex_D = [1,9,17];         % Index of density data used for fit, at least two
        FitIndex_cp = [1,9,17];          % Index of isobaric heat capacity data used for fit, at least two
        FitIndex_TC = [1,8,9,17];   % Index of thermal conductivity data used for fit, at least four
    elseif OilNum == 10
        Material = 'PEB8'; 
        Zc = 0.2600;
        MM_gmol = 640.9;          % give a good guess of the molar mass
        FitIndex_D = [1,34,67,89];         % Index of density data used for fit, at least two
        FitIndex_V = [1, 15, 29,43, 57,71];          % Index of isobaric heat capacity data used for fit, at least two
    elseif OilNum == 11
        Material = 'RENISO ACC HV';
        MM_gmol = 200;          % give a good guess of the molar mass
        FitIndex_D = [1,9,19];         % Index of density data used for fit, at least two
    else
        error(['not defined yet'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('Zc','var')
        if strcmpi(CubicEOS,'PTV'), Zc = 0.2563; elseif strcmpi(CubicEOS,'YFR'), Zc = 0.2640; else, end
    end
    if ~exist('MM_gmol','var')
        MM_gmol = 200;  
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  constants and prepare output file
    fid_out = fopen(['ExpData/',Material,'/Fluid_Constants_',Material,'.txt'],'w');
    fprintf(fid_out,['     Materials;    MM_gmol;  Power_fix;      Zc;      Tc_K;   Dc_kgm3;    pc_MPa;   acentric;  k0_JmolK;  k1_JmolK;    n1_mu;     n2_mu;     n3_mu;    n4_mu;        n1_tc;       n2_tc;       n3_tc;       n4_tc\n']);
    LdenExist = 0; LcpExist = 0; LvisExist = 0; LtcExist = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get 'experimetnal data'   --- density
    Densityfile = ['ExpData/',Material,'/Density.txt'];
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
    CPfile = ['ExpData/',Material,'/HeatCapacity.txt'];
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
    VISfile = ['ExpData/',Material,'/Viscosity.txt'];
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
    TCfile = ['ExpData/',Material,'/ThermalConductivity.txt'];
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

% 
%     const = [Tc_fit, pc_fit, Dc_mol, Zc, af_fit, MM_gmol, k0_JmolK, k1_JmolK, xi_mu,epsilon_kB_K,sigma_nm, xi0, Gamma, qDinv, kB,R];
%     L_WmK_fit = EOSmodel.fit_CEOS_tc_n(n_fit,xdata,const,n_mu,CubicEOS);
%     xdata = [273, 100000]; 
%     L_WmK_fit_all = EOSmodel.fit_CEOS_tc_n(n_fit,xdata,const,n_mu,CubicEOS)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save parameters
    fprintf(fid_out,'%14s;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.4f;%10.4f;%10.3f;%10.3f;%10.6f;%10.6f;%10.6f;%10.6f;%12.6f;%12.6f;%12.6f;%12.6f\n',...
        Material,MM_gmol,PowerConst,Zc,Tc_fit,Dc_fit,pc_fit/1e6,af_fit,k0_JmolK,k1_JmolK,n_mu,n_tc);
    fclose(fid_out);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot the results
    if Lplot
        figure(1);clf;
    
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
    
        set(gcf,'paperunits','centimeters');
        set(gcf,'paperposition',[0 0 18 8]);
        print(gcf,'-dtiff','-r200',['figures/',thisfilename,'_',Material,'.tiff']); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save the parameters to database
    if Lsave2database
        copyfile Classes/Fluid_Constants_Fitted.txt 'Classes/Fluid_Constants_Fitted_2.txt'
        try
            [TheMaterial,MM_gmol,Power_fix,Zc,Tc_fit,Dc_fit,pc_fit,af_fit,k0_fit,k1_fit,...
                n1_mu_fit,n2_mu_fit,n3_mu_fit,n4_mu_fit,n1_tc_fit,n2_tc_fit,n3_tc_fit,n4_tc_fit] = ...
                textread(['ExpData/',Material,'/Fluid_Constants_',Material,'.txt'],'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'delimiter',';');
            [Material_all,MM_gmol_all,Power_fix_all,Zc_all,Tc_fit_all,Dc_fit_all,pc_fit_all,af_fit_all,k0_fit_all,k1_fit_all,...
                n1_mu_fit_all,n2_mu_fit_all,n3_mu_fit_all,n4_mu_fit_all,n1_tc_fit_all,n2_tc_fit_all,n3_tc_fit_all,n4_tc_fit_all] = ...
                textread('Classes/Fluid_Constants_Fitted_2.txt','%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'delimiter',';');
            nM = length(Material_all);
            fid_out = fopen(['./Classes/Fluid_Constants_Fitted.txt'],'w');
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
            delete 'Classes/Fluid_Constants_Fitted_2.txt'
        catch
            copyfile 'Classes/Fluid_Constants_Fitted_2.txt' Classes/Fluid_Constants_Fitted.txt 
            delete 'Classes/Fluid_Constants_Fitted_2.txt'
            error(['Error occurs, check  ./Classes/Fluid_Constants_',Material,'.txt'])
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

