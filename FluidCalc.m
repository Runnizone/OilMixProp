%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  AllEOS = {'PR','SRK','PTV','YFR'};      warning('off'); 
% rmpath('C:\Xiaoxian\Github\OilMixProp\Classes')

%%%  Define cubic EoS, PTV and YFR are recommended %%%
CubicEOS = AllEOS{4}; 


%%% define fluids to study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see Classes/Fluid_Constants(_xxx).txt files for all available fluids
% Both mole fraction and mass fraction are available
% Following are a few examples. 

% Refrigerant = {'R1233zde','Emkarate RL32'};   MassFrac1 = 0.5;  pres_kPa = 1.2e3;  temp_K = 273.15 + 150;  
% Refrigerant = {'CO2','RENISO ACC HV'};   MassFrac1 = 0.6;  pres_kPa = 1e3;  temp_K = 273.15 + 10;  
 Refrigerant = {'propane','R32'};   MassFrac1 = 0.120553544760139;  pres_kPa = 3.4e3;  temp_K = 290;  
% Refrigerant = {'CO2'};   MassFrac1 = 1;  pres_kPa = 1e3;  temp_K = 273.15  + 10;  
% Refrigerant = {'hydrogen'};   MassFrac1 = 1;  pres_kPa = 5e3;  temp_K = 19;  
% Refrigerant = {'CO2','POEiso68'};   MassFrac1 = 0.6;  pres_kPa = 1e4;  temp_K = 273.15 + 10;  
% Refrigerant = {'1-Methylnaphthalene','R32'};   MassFrac1 = 0.9;  pres_kPa = 1.2e3;  temp_K = 273.15 + 150;  
% Refrigerant = {'water','nitrogen'};   MassFrac1 = 0.5;  pres_kPa = 1.0e3;  temp_K = 273.15 + 30;  
% Refrigerant = {'EGLYCOL'};   MassFrac1 = 1;  pres_kPa = 5e5;  temp_K = 273.15 + 1000;  
% Refrigerant = {'PAG68','propane'};   MassFrac1 = 0.8068;  pres_kPa = 80;  temp_K = 232.11; 
% Refrigerant = {'PAG68'};   MassFrac1 = 1;  pres_kPa = 100;  temp_K = 273; 
% Refrigerant = {'propane','R32'};   MoleFrac1 = 0.1499;  pres_kPa = 3.4e3;  temp_K = 290; 
% Refrigerant = {'water','nitrogen','argon'};   MassFrac1 = 0.5; MassFrac2 = 0.2;  pres_kPa = 1.0e3;  temp_K = 273.15 + 30;  
% Refrigerant = {'benzene','ethanol'};   MassFrac1 = 0.5;  pres_kPa = 1.0e3;  temp_K = 273.15 + 30;  

%%% parameter preperation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% You do not need to do anything here
GL = GetGlobals(CubicEOS,Refrigerant);  % obtain fluid constants
ncomp = length(Refrigerant);            % get over all mole fraction Zi
if ~exist('MassFrac1','var') && exist('MoleFrac1','var')
    if ncomp == 3
        MoleFrac = [MoleFrac1,MoleFrac2,1 - MoleFrac1-MoleFrac2]'; 
    elseif ncomp == 2 
        MoleFrac = [MoleFrac1,1 - MoleFrac1]'; 
    else 
        MoleFrac = 1; 
    end
    [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleFrac); 
else
    if ncomp == 3
        MassFrac = [MassFrac1,MassFrac2,1 - MassFrac1 - MassFrac2]'; 
    elseif ncomp == 2 
        MassFrac = [MassFrac1,1 - MassFrac1]'; 
    else 
        MassFrac = 1; 
    end
end
T_K_guess = 0;   % if T to be solved and a good guess is known, otherwise set 0
p_kPa_guess = 0;  % if p to be solved and a good guess is known, otherwise set 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% A list of example inputs %%%%%%%%%%%%%%%%%%%%%%%%%
ff = OilPropm('H','T',temp_K,'P',pres_kPa,MassFrac,GL,0,0);   
% ff = OilPropm('All','T',temp_K,'D',ff.rho_kgm3_all,MassFrac,GL,0,0); 
% ff = OilPropm('All','T',temp_K,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);   
% ff = OilPropm('All','T',temp_K,'H',ff.hh_Jkg_all,MassFrac,GL,0,0); 
% ff = OilPropm('All','P',pres_kPa,'D',ff.rho_kgm3_all,MassFrac,GL,0,0);   
% ff = OilPropm('All','P',pres_kPa,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);
% ff = OilPropm('All','P',pres_kPa,'H',ff.hh_Jkg_all,MassFrac,GL,0,0);  
% if ff.Phase == 2
%     ff = OilPropm('All','T',temp_K,'Q',ff.FracV_mass,MassFrac,GL,0,0);
%     ff = OilPropm('All','P',pres_kPa,'Q',ff.FracV_mass,MassFrac,GL,0,0);                   
% end
%  ff = OilPropm('All','P',pres_kPa,'Q',1,MassFrac,GL,0,0);
%  ff = OilPropm('All','P',pres_kPa,'Q',0,MassFrac,GL,0,0);  
%  ff = OilPropm('All','T',temp_K,'Q',1,MassFrac,GL,0,0);      
%  ff = OilPropm('All','T',temp_K,'Q',0,MassFrac,GL,0,0); 

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
try 
    % disp([num2str(MoleFrac1,'%0.2f'),' ',Refrigerant{1},' + ',num2str(1-MoleFrac1,'%0.2f'),' ' ,Refrigerant{2},' at ',num2str(temp_K,'%0.3f'),' K and ',num2str(pres_kPa/1e3,'%0.3f'),' MPa',]);
    if ncomp == 2
        disp([num2str(MassFrac(1),'%0.2f'),' ',Refrigerant{1},' + ',num2str(MassFrac(2),'%0.2f'),' ' ,Refrigerant{2},' in mass frac. with ' ,CubicEOS,' EoS ']);
    else
        disp([Refrigerant{1},': ' ,CubicEOS]);
    end
    
    if length(ff.Phase) == 2
        fprintf(' T:             %10.6f   K  \n ',ff.T_K);
        fprintf('p:             %10.6f   kPa \n ',ff.p_Pa/1e3);
        fprintf('MM:            %10.3f   g/mol \n ',ff.MM_gmol);
        fprintf('Phase:               %4s \n\n',ff.Phase);
        fprintf('Phase behavior in mass fraction.  Vapor Frac: %8.6f\n', ff.FracV_mass);
        fprintf('          %10s      %10s  \n' ,Refrigerant{:}  )
        fprintf('Liquid:        %10.6f  %10.6f  \n',ff.MassF_Li);
        fprintf('Vapor:         %10.6f  %10.6f  \n\n',ff.MassF_Vi);
    %     fprintf('Phase behavior in mole fraction.  Vapor Frac: %8.6f\n', ff.FracV_mole);
    %     fprintf('          %10s      %10s  \n' ,Refrigerant{:}  )
    %     fprintf('Liquid:        %10.6f  %10.6f  \n',ff.MoleF_Li);
    %     fprintf('Vapor:         %10.6f  %10.6f  \n\n',ff.MoleF_Vi);
        fprintf('       ---- Properties in each phase  ----       \n ');
        fprintf('                    Liquid       Vapor\n '   )
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
catch
    disp(ffi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if refpropm are in your PC, you can do the following calculation
% Contact Xiaoxian Yang if you need help on this
try
    if ncomp == 1 
        D_ref = refpropm('D','T',temp_K,'P',pres_kPa,Refrigerant{:});
        Cp_ref = refpropm('C','T',temp_K,'P',pres_kPa,Refrigerant{:});
        Cv_ref = refpropm('O','T',temp_K,'P',pres_kPa,Refrigerant{:});
        sos_ref = refpropm('A','T',temp_K,'P',pres_kPa,Refrigerant{:});
        vis_ref = refpropm('V','T',temp_K,'P',pres_kPa,Refrigerant{:});
        tc_ref = refpropm('L','T',temp_K,'P',pres_kPa,Refrigerant{:});
        drhodP = refpropm('R','T',temp_K,'P',pres_kPa,Refrigerant{:});
        drhodT = refpropm('W','T',temp_K,'P',pres_kPa,Refrigerant{:});
        dpdT = refpropm('#','T',temp_K,'P',pres_kPa,Refrigerant{:});
    elseif ncomp == 2 
        D_ref = refpropm('D','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        Cp_ref = refpropm('C','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        Cv_ref = refpropm('C','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        sos_ref = refpropm('A','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        vis_ref = refpropm('V','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        tc_ref = refpropm('L','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        drhodP = refpropm('R','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        drhodT = refpropm('W','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
        dpdT = refpropm('#','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},MassFrac);
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
