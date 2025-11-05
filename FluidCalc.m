%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  AllEOS = {'PR','SRK','PTV','YFR'};      warning('off'); 
% rmpath('C:\Xiaoxian\Github\OilMixProp\Classes')
rmpath('C:\Xiaoxian\Github\HEoS\Classes')

%%%  Define cubic EoS, PTV and YFR are recommended %%%
CubicEOS = AllEOS{4}; 

%%% define fluids to study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see Classes/Fluid_Constants(_xxx).txt files for all available fluids
% Both mole fraction and mass fraction are available
% Following are a few examples. 

% TheFluid = {'R1233zde','Emkarate RL32'};   MassFrac1 = 0.5;  pres_kPa = 1.2e3;  temp_K = 273.15 + 150;  
% TheFluid = {'CO2','RENISO ACC HV'};   MassFrac1 = 0.6;  pres_kPa = 1e3;  temp_K = 273.15 + 10;  
% TheFluid = {'propane','R32'};   MassFrac1 = 0.4588;  pres_kPa = 0.84e3;  temp_K = 283.15;  
% TheFluid = {'CO2'};   MassFrac1 = 1;  pres_kPa = 1e3;  temp_K = 273.15 + 10 ;  
% TheFluid = {'CO2','TUD MO10'};   MassFrac1 = 0.99;  pres_kPa = 1e3;  temp_K = 273.15 + 10 ;  
% TheFluid = {'hydrogen'};   MassFrac1 = 1;  pres_kPa = 5e3;  temp_K = 19;  
% TheFluid = {'CO2','POEiso68'};   MassFrac1 = 0.6;  pres_kPa = 1e4;  temp_K = 273.15 + 10;  
% TheFluid = {'1-Methylnaphthalene','R32'};   MassFrac1 = 0.9;  pres_kPa = 1.2e3;  temp_K = 273.15 + 150;  
% TheFluid = {'water','nitrogen'};   MassFrac1 = 0.5;  pres_kPa = 1.0e3;  temp_K = 273.15 + 30;  
% TheFluid = {'EGLYCOL'};   MassFrac1 = 1;  pres_kPa = 5e5;  temp_K = 273.15 + 1000;  
% TheFluid = {'PAG68','propane'};   MassFrac1 = 0.8068;  pres_kPa = 80;  temp_K = 232.11; 
% TheFluid = {'PAG68'};   MassFrac1 = 1;  pres_kPa = 100;  temp_K = 273; 
% TheFluid = {'propane','R32'};   MoleFrac1 = 0.1499;  pres_kPa = 3.4e1;  temp_K = 290; 
% TheFluid = {'water','nitrogen','argon'};   MassFrac1 = 0.5; MassFrac2 = 0.2;  pres_kPa = 1.0e3;  temp_K = 273.15 + 30;  
% TheFluid = {'benzene','ethanol'};   MassFrac1 = 0.5;  pres_kPa = 1.0e3;  temp_K = 273.15 + 30;  
% TheFluid = {'Shrieve POE68','R32'}; MassFrac1 = 0.6;  pres_kPa = 15e3;  temp_K = 273.15;  
% TheFluid = {'HATCOL 4467','propane'}; MassFrac1 = 0.55;  pres_kPa = 0.1e3;  temp_K = 298.15;  
% TheFluid = {'TUD POE68','CO2'}; MassFrac1 = 0.5;  pres_kPa = 11e3;  temp_K = 350.15;  
% TheFluid = {'13BUTADIENE','ARGON'}; MassFrac1 = 0.5;  pres_kPa = 300;  temp_K = 400;  
% TheFluid = {'13BUTADIENE','N-butyric acid'}; MassFrac1 = 0.5;  pres_kPa = 300;  temp_K = 400;  
% TheFluid = {'TUD POE68','CO2'}; MassFrac1 = 1 - 0.6172;  pres_kPa = 13000;  temp_K = 70.32 + 273.15;  
% TheFluid = {'Pentane','propane'}; MassFrac1 = 0.04;  pres_kPa = 800;  temp_K = 293;  
% TheFluid = {'decane','water'}; MoleFrac1 = 0.5;  pres_kPa = 2900;  temp_K = 548;  
% TheFluid = {'water','1-butanol'}; MoleFrac1 = 0.7;  pres_kPa = 101;  temp_K = 320;  
% TheFluid = {'CO2','PEC8'};   MoleFrac1 = 0.95;  pres_kPa = 8e3;  temp_K = 303.15;  
% TheFluid = {'PEC8'};   MassFrac1 = 0.1;  pres_kPa = 100;  temp_K = 283.15;  
% TheFluid = {'CO2','methanol'};   MoleFrac1 = 0.95;  pres_kPa = 8e3;  temp_K = 298.15;  
% TheFluid = {'R134a','nitrogen'}; MoleFrac1 = 0.8;    temp_K = 300;  pres_kPa  = 1000; 
% TheFluid = {'CO2','ethane'}; MoleFrac1 = 0.3;    temp_K = 320;  pres_kPa  = 100;  D_kgm3 = 37.0394;
% TheFluid = {'CO2'};   MassFrac1 = 1;  pres_kPa = 3558.3162;  temp_K = 3.341282000000000e+02 ;    D_kgm3 = 66.0147;
% TheFluid = {'R134a','nitrogen','CO2'}; MoleFrac1 = 0.8;  MoleFrac2 = 0.1;    temp_K = 300;  pres_kPa  = 1000; 


% TheFluid = {...
%     'METHANE',      'ETHANE',       'PROPANE',      'BUTANE',      ...
%     'ISOBUTAN',    'PENTANE',     'IPENTANE',     'CO2',   ...
%     'HEXANE',      'HEPTANE',     'OCTANE',      'DECANE',      ...
%     'ETHYLENE',     'PROPYLEN',    '1BUTENE',      'C2BUTENE',     ...
%     'T2BUTENE',     'IBUTENE',    'ACETYLENE',    'CYCLOHEX'      ...
% };   MoleFrac0 = ones(20,1) * 0.05;   MoleFrac1 = 0; pres_kPa = 3.4e1;  temp_K = 290; 

TheFluid = { 'METHANE',   'ETHANE',   'ETHYLENE', 'PROPYLEN'}; MoleFrac0 = ones(4,1) * 0.25;   MoleFrac1 = 0;    temp_K = 1400;  pres_kPa  = 7e3; 


%%% parameter preperation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% You do not need to do anything here
GL = GetGlobals(CubicEOS,TheFluid);  % obtain fluid constants
ncomp = length(TheFluid);            % get over all mole fraction Zi
if ~exist('MassFrac1','var') && exist('MoleFrac1','var')
    if ncomp == 3
        MoleFrac = [MoleFrac1,MoleFrac2,1 - MoleFrac1-MoleFrac2]'; 
    elseif ncomp == 2 
        MoleFrac = [MoleFrac1,1 - MoleFrac1]'; 
    elseif  ncomp == 1
        MoleFrac = 1; 
    else
        MoleFrac = MoleFrac0;
    end
    [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleFrac); 
else
    if ncomp == 3
        MassFrac = [MassFrac1,MassFrac2,1 - MassFrac1 - MassFrac2]'; 
    elseif ncomp == 2 
        MassFrac = [MassFrac1,1 - MassFrac1]'; 
    elseif  ncomp == 1
        MassFrac = 1; 
    else
        MassFrac = MassFrac0;
    end
end
T_K_guess = 0;   % if T to be solved and a good guess is known, otherwise set 0
p_kPa_guess = 0;  % if p to be solved and a good guess is known, otherwise set 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% A list of example inputs %%%%%%%%%%%%%%%%%%%%%%%%%
ff = OilPropm('All','T',temp_K,'P',pres_kPa,MassFrac,GL,0,0);   

% ff = OilPropm('All','T',temp_K,'D',D_kgm3,MassFrac,GL,0,0); 

% ff = OilPropm('All','T',temp_K,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);   
% ff = OilPropm('All','T',temp_K,'H',ff.hh_Jkg_all,MassFrac,GL,0,0); 
% ff = OilPropm('All','P',pres_kPa,'D',ff.rho_kgm3_all,MassFrac,GL,0,0);   
% ff = OilPropm('All','P',pres_kPa,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);
% ff = OilPropm('All','P',pres_kPa,'H',ff.hh_Jkg_all,MassFrac,GL,0,0);  
% if ff.Phase == 2
%     ff = OilPropm('All','T',temp_K,'Q',ff.FracV_mass,MassFrac,GL,0,0);
%     ff = OilPropm('All','P',pres_kPa,'Q',ff.Frac_mass(4),MassFrac,GL,0,0);                   
% end
%  ff = OilPropm('All','P',pres_kPa,'Q',1,MassFrac,GL,0,0);
 % ff = OilPropm('All','P',pres_kPa,'Q',0,MassFrac,GL,0,0);  
 % ff = OilPropm('All','T',temp_K,'Q',1,MassFrac,GL,0,0);    
 % ff = OilPropm('All','T',temp_K,'Q',0,MassFrac,GL,0,0); 

% ti
%  ff = OilPropm('All','T',temp_K,'Q',1,MassFrac,GL,0,0);    
%  ff = OilPropm('All','P',ff.p_Pa/1000,'Q',1,MassFrac,GL,0,0);
% toc
% ffi = OilPropm('A+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['SoS: ',num2str(ffi),' m/s']);
% ffi = OilPropm('C+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Cp: ',num2str(ffi),' xxx']);
% ffi = OilPropm('D+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Den: ',num2str(ffi),' xxx']);
% ffi = OilPropm('H+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['HH: ',num2str(ffi),' xxx']);
% ffi = OilPropm('K+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Kappa: ',num2str(ffi),' xxx']);
% ffi = OilPropm('L+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Vis: ',num2str(ffi),' xxx']);
% ffi = OilPropm('O+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Cv: ',num2str(ffi),' xxx']);
% ffi = OilPropm('P+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Pres: ',num2str(ffi),' xxx']);
% ffi = OilPropm('Q+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['FracV: ',num2str(ffi),' xxx']);
% ffi = OilPropm('S+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['SS: ',num2str(ffi),' xxx']);
% ffi = OilPropm('T+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Temp: ',num2str(ffi),' xxx']);
% ffi = OilPropm('U+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['UU: ',num2str(ffi),' xxx']);
% ffi = OilPropm('V+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Vol: ',num2str(ffi),' xxx']);
% ffi = OilPropm('X+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); % disp(['XX: ',num2str(ffi),' xxx']);
% ffi = OilPropm('Z+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['Z: ',num2str(ffi),' xxx']);
% ffi = OilPropm('#+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['xxxx: ',num2str(ffi),' xxx']);
% ffi = OilPropm('R+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['xxxx: ',num2str(ffi),' xxx']);
% ffi = OilPropm('W+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); disp(['xxxx: ',num2str(ffi),' xxx']);
% % 

% dT = 0.001;
% ff0 = OilPropm('All','T',temp_K,'P',pres_kPa,MassFrac,GL,0,0);
% ff1 = OilPropm('All','T',temp_K+dT,'D',ff0.rho_kgm3_all,MassFrac,GL,0,0);  
% dp_dT = (ff1.p_Pa - ff0.p_Pa)/dT / 1000; 
% 
% tic
% ff = OilPropm('All','P',pres_kPa,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);
% toc
% tic
% ff = OilPropm('All+','P',pres_kPa,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% show the results  %%%%%%%%%%%%%%%%%%%%%%%%%
% try 
    % disp([num2str(MoleFrac1,'%0.2f'),' ',TheFluid{1},' + ',num2str(1-MoleFrac1,'%0.2f'),' ' ,TheFluid{2},' at ',num2str(temp_K,'%0.3f'),' K and ',num2str(pres_kPa/1e3,'%0.3f'),' MPa',]);
    if ncomp == 2
        disp([num2str(MassFrac(1),'%0.2f'),' ',TheFluid{1},' + ',num2str(MassFrac(2),'%0.2f'),' ' ,TheFluid{2},' in mass frac. with ' ,CubicEOS,' EoS ']);
    elseif ncomp == 1
        disp([TheFluid{1},': ' ,CubicEOS]);
    else
        disp(['Mixture : ' ,CubicEOS]);
    end
    
    if length(ff.Phase) == 2
        fprintf(' T:             %10.6f   K  \n ',ff.T_K);
        fprintf('p:             %10.6f   kPa \n ',ff.p_Pa/1e3);
        fprintf('MM:            %10.3f   g/mol \n ',ff.MM_gmol);
        fprintf('Phase:               %4s \n\n',ff.Phase);
        % fprintf('Phase behavior in mass fraction. \n');
        % fprintf('    Solid       Liquid2     Liquid1     Vapor\n');
        % fprintf('%12.6f%12.6f%12.6f%12.6f\n', ff.Frac_mass);
        % fprintf('          %10s      %10s  \n' ,TheFluid{:}  )
        % fprintf('Solid:         %10.6f  %10.6f  \n',ff.MassF(:,1));
        % fprintf('Liquid2:       %10.6f  %10.6f  \n',ff.MassF(:,2));
        % fprintf('Liquid1:       %10.6f  %10.6f  \n',ff.MassF(:,3));
        % fprintf('Vapor:         %10.6f  %10.6f  \n\n',ff.MassF(:,4));

        fprintf('Phase behavior in mole fraction. \n');
        fprintf('                Solid       Liquid2     Liquid1     Vapor\n');
        fprintf('%10s%12.6f%12.6f%12.6f%12.6f\n\n','all:', ff.Frac_mole);
        % fprintf('          %10s      %10s  \n' ,TheFluid{:}  )

        for icomp = 1:ncomp
            fprintf('%10s%12.6f%12.6f%12.6f%12.6f\n',TheFluid{icomp}, ff.MoleF(icomp,:));

        end


        % fprintf('Solid:         %10.6f  %10.6f  \n',ff.MoleF(:,1));
        % fprintf('Liquid2:       %10.6f  %10.6f  \n',ff.MoleF(:,2));
        % fprintf('Liquid1:       %10.6f  %10.6f  \n',ff.MoleF(:,3));
        % fprintf('Vapor:         %10.6f  %10.6f  \n\n',ff.MoleF(:,4));

        fprintf('       ---- Properties in each phase  ----       \n ');
        if strcmpi(ff.Phase,'LV')
            fprintf('                Liquid1       Vapor  \n '   )
        elseif strcmpi(ff.Phase,'LL')
            fprintf('                  Liquid2     Liquid1 \n '   )
        else
            fprintf(' to be completed \n '   )
        end
        fprintf('Z:             %10.6f  %10.6f\n ',ff.Z);
        fprintf('MM:            %10.3f  %10.3f   g/mol \n ',ff.MM_phase_gmol);
        fprintf('rho:           %10.3f  %10.3f   kg/m3 \n ',ff.rho_kgm3);
        fprintf('Entropy:       %10.3f  %10.3f   kJ/K/kg\n ',ff.ss_JkgK/1000);
        fprintf('Enthalpy:      %10.3f  %10.3f   kJ/kg\n ',ff.hh_Jkg/1000);
        fprintf('cp:            %10.3f  %10.3f   kJ/K/kg\n ',ff.cp_JkgK/1000);
        fprintf('cv:            %10.3f  %10.3f   kJ/K/kg\n ',ff.cv_JkgK/1000);
        fprintf('SoS:           %10.3f  %10.3f   m/s\n ',ff.sos_ms);
        fprintf('vis:           %10.6f  %10.6f   mPa s\n ',ff.vis_Pas*1000);
        fprintf('TC:            %10.6f  %10.6f   W/m/K\n ',ff.lambda_WmK);
        fprintf('Kappa          %10.6f  %10.6f   \n ',ff.kappa);
        fprintf('d(rho)/dP:     %10.6f  %10.6f   kg/m3/kPa \n ',ff.drhokg_dpkPa);
        fprintf('d(rho)/dT:     %10.6f  %10.6f   kg/m3/K \n ',ff.drhokg_dT);
        fprintf('dp/dT_rho:     %10.6f  %10.6f   kPa/K \n ',ff.dp_dT_kPaK);
        fprintf('\n   -- Properties with all phases combined --\n ');
        fprintf('rho_all:             %10.3f         kg/m3 \n ', ff.rho_kgm3_all  );
        fprintf('Entropy_all:         %10.3f         kJ/K/kg \n ', ff.ss_JkgK_all/1000   );
        fprintf('Enthalpy_all:        %10.3f         kJ/kg \n ', ff.hh_Jkg_all/1000   );
        fprintf('cp_all:              %10.3f         kJ/K/kg \n ', ff.cp_JkgK_all/1000   );
        fprintf('cv_all:              %10.3f         kJ/K/kg \n ', ff.cv_JkgK_all/1000   );
        fprintf('Kappa_all            %10.6f     \n ',ff.kappa_all);
    else
    %     fprintf('Fluid at %1s phase \n ',ff.Phase);
        fprintf('Fluid highly likely in the %1s phase\n ',ff.Phase);
        fprintf('T:             %10.6f   K  \n ',ff.T_K);
        fprintf('p:             %10.6f   kPa \n ',ff.p_Pa/1e3);
        fprintf('Z:             %10.6f\n ',ff.Z);
        fprintf('MM:            %10.3f   g/mol \n ',ff.MM_gmol);
        fprintf('rho_kg:        %10.3f   kg/m3 \n ',ff.rho_kgm3);
    %     fprintf('rho_mol:       %10.3f   kg/mol \n ',ff.rho_molm3);
        fprintf('d(rho)/dP:     %10.6f   kg/m3/kPa \n ',ff.drhokg_dpkPa);
        fprintf('d(rho)/dT:     %10.6f   kg/m3/K \n ',ff.drhokg_dT);
        fprintf('dp/dT_rho:     %10.6f   kPa/K \n ',ff.dp_dT_kPaK);
        fprintf('Entropy:       %10.6f   kJ/K/kg\n ',ff.ss_JkgK/1000);
    %     fprintf('Entropy:       %10.6f   J/K/mol\n ',ff.ss_JmolK);
        fprintf('Enthalpy:      %10.3f   kJ/kg\n ',ff.hh_Jkg/1000);
    %     fprintf('Enthalpy:      %10.3f   J/mol\n ',ff.hh_Jmol);
        fprintf('cp:            %10.6f   kJ/K/kg\n ',ff.cp_JkgK/1000);
        fprintf('cv:            %10.6f   kJ/K/kg\n ',ff.cv_JkgK/1000);
        fprintf('SoS:           %10.3f   m/s\n ',ff.sos_ms);
        fprintf('vis:           %10.6f   mPa s\n ',ff.vis_Pas*1000);
        fprintf('TC:            %10.6f   W/m/K\n ',ff.lambda_WmK);
        fprintf('Kappa          %10.6f   \n ',ff.kappa);
    end
% catch
%     disp(ffi);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if refpropm are in your PC, you can do the following calculation
% Contact Xiaoxian Yang if you need help on this
try
    if ncomp == 1 
        D_ref = refpropm('D','T',temp_K,'P',pres_kPa,TheFluid{:});
        Cp_ref = refpropm('C','T',temp_K,'P',pres_kPa,TheFluid{:});
        Cv_ref = refpropm('O','T',temp_K,'P',pres_kPa,TheFluid{:});
        sos_ref = refpropm('A','T',temp_K,'P',pres_kPa,TheFluid{:});
        vis_ref = refpropm('V','T',temp_K,'P',pres_kPa,TheFluid{:});
        tc_ref = refpropm('L','T',temp_K,'P',pres_kPa,TheFluid{:});
        drhodP = refpropm('R','T',temp_K,'P',pres_kPa,TheFluid{:});
        drhodT = refpropm('W','T',temp_K,'P',pres_kPa,TheFluid{:});
        dpdT = refpropm('#','T',temp_K,'P',pres_kPa,TheFluid{:});
    elseif ncomp == 2 
        D_ref = refpropm('D','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        Cp_ref = refpropm('C','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        Cv_ref = refpropm('O','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        sos_ref = refpropm('A','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        vis_ref = refpropm('V','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        tc_ref = refpropm('L','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        drhodP = refpropm('R','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        drhodT = refpropm('W','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
        dpdT = refpropm('#','T',temp_K,'P',pres_kPa,TheFluid{1},TheFluid{2},MassFrac);
    end
        fprintf('\n\n----- Refprop Calculations --- \n ');
        fprintf('rho_kg:        %10.3f   kg/m3 \n ',D_ref);
        fprintf('cp:            %10.6f   kJ/K/kg\n ',Cp_ref/1000);
        fprintf('cv:            %10.6f   kJ/K/kg\n ',Cv_ref/1000);
        fprintf('SoS:           %10.3f   m/s\n ',sos_ref);
%         fprintf('sr:            %10.3f   J/(kg K)\n ',sr);
        fprintf('vis:           %10.6f   mPa s\n ',vis_ref*1000);
        fprintf('TC:            %10.6f   W/m/K\n ',tc_ref);
        fprintf('d(rho)/dP:     %10.6f   kg/m3/kPa \n ',drhodP);
        fprintf('d(rho)/dT:     %10.6f   kg/m3/K \n ',drhodT);
        fprintf('dp/dT:         %10.6f   kPa/K \n ',dpdT);
        fprintf('------------------------------- \n ');
catch
    fprintf('\n\n--- Refprop Calculation fail --=  \n ');
end
