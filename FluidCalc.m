%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  AllEOS = {'PR','SRK','PTV','YR'};    %  warning('off'); 
CubicEOS = AllEOS{4}; 
%%% define fluids to study
Refrigerant = {'R1233zde','Emkarate RL32'};   MoleFrac1 = 0.900;  pres = 1.2e6;  temp = 273.15 + 150;  
% Refrigerant = {'CO2','RENISO ACC HV'};   MoleFrac1 = 0.55;  pres = 1.2e6;  temp = 373.15;  
% Refrigerant = {'propane','R32'};   MoleFrac1 = 0.5;  pres = 0.85e6;  temp = 283.15;  
% Refrigerant = {'CO2'};   MoleFrac1 = 1;  pres = 9.2e6;  temp = 273.15 + 90;  

%%% parameter preperation
GL = GetGlobals(CubicEOS,Refrigerant);  % obtain fluid constants
ncomp = length(Refrigerant);            % get over all mole fraction Zi
if ncomp == 3
    Zi = [MoleFrac1,MoleFrac2,1 - MoleFrac1-MoleFrac2]'; 
elseif ncomp == 2 
    Zi = [MoleFrac1,1 - MoleFrac1]'; 
else 
    Zi = 1; 
end
T_K_0 = 0;   % if T to be solved and a good guess is known, otherwise set 0
p_Pa_0 = 0;  % if p to be solved and a good guess is known, otherwise set 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% call the function %%%%%%%%%%%%%%%%%%%%%%%%%
ff = OilPropm('T',temp,'P',pres,Zi,GL,0,0);
if ff.Phase == 2
    ff = OilPropm('T',temp,'Q',ff.FracV,Zi,GL,0,0);
    ff = OilPropm('P',pres,'Q',ff.FracV,Zi,GL,temp,0);                    % sometimes rely on a good guess.
end
ff = OilPropm('T',temp,'D',ff.rho_molar_all,Zi,GL,300,1e6);     
ff = OilPropm('P',pres,'D',ff.rho_molar_all,Zi,GL,300,1e6);   
ff = OilPropm('P',pres,'S',ff.ss_JmolK_all,Zi,GL,300,1e6);
ff = OilPropm('P',pres,'H',ff.hh_Jmol_all,Zi,GL,300,1e6);  

% ff = OilPropm('T',temp,'Q',1,Zi,GL,0,0);   % V
% ff = OilPropm('T',temp,'Q',0,Zi,GL,0,0);   % L
% ff = OilPropm('P',pres,'Q',1,Zi,GL,0,0);   % V                  % sometimes rely on a good guess.  
% ff = OilPropm('P',pres,'Q',0,Zi,GL,0,0);   % L                   % sometimes rely on a good guess. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% show the results  %%%%%%%%%%%%%%%%%%%%%%%%%
% disp([num2str(MoleFrac1,'%0.2f'),' ',Refrigerant{1},' + ',num2str(1-MoleFrac1,'%0.2f'),' ' ,Refrigerant{2},' at ',num2str(temp,'%0.3f'),' K and ',num2str(pres/1e6,'%0.3f'),' MPa',]);
if ncomp == 2
    disp([num2str(MoleFrac1,'%0.2f'),' ',Refrigerant{1},' + ',num2str(1-MoleFrac1,'%0.2f'),' ' ,Refrigerant{2},': ' ,CubicEOS]);
else
    disp([Refrigerant{1},': ' ,CubicEOS]);
end

if length(ff.Phase) == 2
    fprintf(' T:             %10.6f   K  \n ',ff.T_K);
    fprintf('p:             %10.6f   MPa \n ',ff.p_Pa/1e6);
    fprintf('MM:            %10.3f   g/mol \n ',ff.mm_gmol);
    fprintf('\nPhase:      %4s,    Vapor Frac: %8.6f\n\n',ff.Phase, ff.FracV);
    fprintf('          %10s      %10s  \n' ,Refrigerant{:}  )
    fprintf('Liquid:        %10.6f  %10.6f  \n',ff.MoleF_Li);
    fprintf('Vapor:         %10.6f  %10.6f  \n\n',ff.MoleF_Vi);
    fprintf('                    Liquid       Vapor\n '   )
    fprintf('Z:             %10.6f  %10.6f\n ',ff.Z);
    fprintf('MM_phase:      %10.3f  %10.3f   g/mol \n ',ff.mm_phase_gmol);
    fprintf('rho_kg:        %10.3f  %10.3f   kg/m3 \n ',ff.rho_kgm3);
%     fprintf('rho_mol:       %10.3f  %10.3f   mol/m3 \n ',ff1.rho_molarm3);
    fprintf('rho_all:             %10.3f         mol/m3 \n ', ff.rho_molar_all  );
    fprintf('d(rho)/dP:     %10.6f  %10.6f   kg/m3/kPa \n ',ff.drho_dp.*ff.mm_phase_gmol);
    fprintf('d(rho)/dT:     %10.6f  %10.6f   kg/m3/K \n ',ff.drho_dT.*ff.mm_phase_gmol/1000);
    fprintf('Entropy:       %10.3f  %10.3f   kJ/K/kg\n ',ff.ss_JkgK/1000);
    fprintf('Entropy_all:         %10.3f         J/K/mol \n ', ff.ss_JmolK_all   );
    fprintf('Enthalpy:      %10.3f  %10.3f   kJ/kg\n ',ff.hh_Jkg/1000);
    fprintf('Enthalpy_all:        %10.3f         J/mol \n ', ff.hh_Jmol_all   );
    fprintf('cp:            %10.3f  %10.3f   kJ/K/kg\n ',ff.cp_JkgK/1000);
    fprintf('cv:            %10.3f  %10.3f   kJ/K/kg\n ',ff.cv_JkgK/1000);
    fprintf('SoS:           %10.3f  %10.3f   m/s\n ',ff.sos_ms);
    fprintf('vis:           %10.6f  %10.6f   mPa s\n ',ff.vis_Pas*1000);
    fprintf('TC:            %10.6f  %10.6f   W/m/K\n ',ff.lambda_WmK);
    fprintf('Kappa          %10.6f  %10.6f   \n ',ff.kappa);
else
%     fprintf('Fluid at %1s phase \n ',ff.Phase);
    fprintf('Fluid at single phase, probably %1s \n ',ff.Phase);
    fprintf('T:             %10.6f   K  \n ',ff.T_K);
    fprintf('p:             %10.6f   MPa \n ',ff.p_Pa/1e6);
    fprintf('Z:             %10.6f\n ',ff.Z);
    fprintf('MM:            %10.3f   g/mol \n ',ff.mm_gmol);
    fprintf('rho_kg:        %10.3f   kg/m3 \n ',ff.rho_kgm3);
    fprintf('rho_mol:       %10.3f   kg/mol \n ',ff.rho_molar);
    fprintf('d(rho)/dP:     %10.6f   kg/m3/kPa \n ',ff.drho_dp.*ff.mm_gmol);
    fprintf('d(rho)/dT:     %10.6f   kg/m3/K \n ',ff.drho_dT.*ff.mm_gmol/1000);
    fprintf('Entropy:       %10.6f   kJ/K/kg\n ',ff.ss_JkgK/1000);
    fprintf('Entropy:       %10.6f   J/K/mol\n ',ff.ss_JmolK);
    fprintf('Enthalpy:      %10.3f   kJ/kg\n ',ff.hh_Jkg/1000);
    fprintf('Enthalpy:      %10.3f   J/mol\n ',ff.hh_Jmol);
    fprintf('cp:            %10.6f   kJ/K/kg\n ',ff.cp_JkgK/1000);
    fprintf('cv:            %10.6f   kJ/K/kg\n ',ff.cv_JkgK/1000);
    fprintf('SoS:           %10.3f   m/s\n ',ff.sos_ms);
    fprintf('vis:           %10.6f   mPa s\n ',ff.vis_Pas*1000);
    fprintf('TC:            %10.6f   W/m/K\n ',ff.lambda_WmK);
    fprintf('Kappa          %10.6f   \n ',ff.kappa);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    if ncomp == 1 
        D_ref = refpropm('D','T',temp,'P',pres/1000,Refrigerant{:});
        Cp_ref = refpropm('C','T',temp,'P',pres/1000,Refrigerant{:});
        sos_ref = refpropm('A','T',temp,'P',pres/1000,Refrigerant{:});
        vis_ref = refpropm('V','T',temp,'P',pres/1000,Refrigerant{:});
        tc_ref = refpropm('L','T',temp,'P',pres/1000,Refrigerant{:});
        drhodP = refpropm('R','T',temp,'P',pres/1000,Refrigerant{:});
        drhodT = refpropm('W','T',temp,'P',pres/1000,Refrigerant{:});
    elseif ncomp == 2 
        [MM_mix_gmol,MassFrac] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,Zi);  
        D_ref = refpropm('D','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
        Cp_ref = refpropm('C','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
        sos_ref = refpropm('A','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
        vis_ref = refpropm('V','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
        tc_ref = refpropm('L','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
        drhodP = refpropm('R','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
        drhodT = refpropm('W','T',temp,'P',pres/1000,Refrigerant{1},Refrigerant{2},MassFrac);
    end
        fprintf('\n\n----- Refprop Calculations --- \n ');
        fprintf('rho_kg:        %10.3f   kg/m3 \n ',D_ref);
        fprintf('cp:            %10.6f   kJ/K/kg\n ',Cp_ref/1000);
        fprintf('SoS:           %10.3f   m/s\n ',sos_ref);
        fprintf('vis:           %10.6f   mPa s\n ',vis_ref*1000);
        fprintf('TC:            %10.6f   W/m/K\n ',tc_ref);
        fprintf('d(rho)/dP:     %10.6f   kg/m3/kPa \n ',drhodP);
        fprintf('d(rho)/dT:     %10.6f   kg/m3/K \n ',drhodT);
        fprintf('------------------------------- \n ');
catch
    fprintf('\n\n--- Refprop Calculation fail --=  \n ');
end
