% 0507版本
function text2put=fluidcalculate_g1_out(CubicEOS_g1,Refrigerant,fraction,percent,para_of_g1,para1,para2,T_K_guess,p_kPa_guess,~);
%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('tmp.txt', 'w');
if fid == -1
    warndlg('无法创建临时txt文件.', 'ERROR');
end

fraction_type = {'mass','mole'};
MassFrac1 = percent(1);
if length(Refrigerant)>2
    MassFrac2 = percent(2);
end
% pres_kPa = 1e3;  temp_K = 273.15 + 10;
format short;warning('off');

%%% parameter preperation
% CubicEOS_g1='YR';Refrigerant={'13BUTADIENE'};
GL = GetGlobals(CubicEOS_g1,Refrigerant);  % obtain fluid constants
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

% T_K_guess = 0;   % if T to be solved and a good guess is known, otherwise set 0
% p_kPa_guess = 0;  % if p to be solved and a good guess is known, otherwise set 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch para_of_g1
    case 11 %TP
        input1='T';input2='P';
    case 12 %TD
        input1='T';input2='D';
    case 13 %TS
        input1='T';input2='S';
    case 14 %TH
        input1='T';input2='H';
    case 15 %TQ
        input1='T';input2='Q';
    case 21 %PT
        input1='P';input2='T';
    case 22 %PD
        input1='P';input2='D';
    case 23 %PS
        input1='P';input2='S';
    case 24 %PH
        input1='P';input2='H';
    case 25 %PQ
        input1='P';input2='Q';
end

% input1='T';input2='D';para1=222;para2=0.029;T_K_guess = 0;p_kPa_guess = 0;MassFrac = 1;


ff = OilPropm('ALL',input1,para1,input2,para2,MassFrac,GL,T_K_guess,p_kPa_guess);
% ff = OilPropm('All','T',temp_K,'P',pres_kPa,MassFrac,GL,0,0);
%
% tic
% ff = OilPropm('All','P',pres_kPa,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);
% toc
% tic
% ff = OilPropm('All+','P',pres_kPa,'S',ff.ss_JkgK_all,MassFrac,GL,0,0);
% toc


%%%%%%%%%%%%%%%%%%%%%%%%% call the function %%%%%%%%%%%%%%%%%%%%%%%%%
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

%
%  ff = OilPropm('All','P',pres_kPa,'Q',1,MassFrac,GL,0,0);
%  ff = OilPropm('All','P',pres_kPa,'Q',0,MassFrac,GL,0,0);
%  ff = OilPropm('All','T',temp_K,'Q',1,MassFrac,GL,0,0);
%  ff = OilPropm('All','T',temp_K,'Q',0,MassFrac,GL,0,0);
%


% temp_K=para1;pres_kPa=para2;
% ffi = OilPropm('A+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['SoS: ',num2str(ffi),' m/s \n']);
% ffi = OilPropm('C+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Cp: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('D+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Den: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('H+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['HH: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('K+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Kappa: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('L+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Vis: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('O+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Cv: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('P+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Pres: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('Q+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['FracV: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('S+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['SS: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('T+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Temp: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('U+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['UU: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('V+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Vol: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('X+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); % % fprintf(fid,['XX: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('Z+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['Z: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('#+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['xxxx: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('R+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['xxxx: ',num2str(ffi),' xxx \n']);
% ffi = OilPropm('W+','T',temp_K,'P',pres_kPa,MassFrac,GL,T_K_guess,p_kPa_guess); %fprintf(fid,['xxxx: ',num2str(ffi),' xxx \n']);
% %

% dT = 0.001;
% ff0 = OilPropm('All','T',temp_K,'P',pres_kPa,MassFrac,GL,0,0);
% ff1 = OilPropm('All','T',temp_K+dT,'D',ff0.rho_kgm3_all,MassFrac,GL,0,0);
% dp_dT = (ff1.p_Pa - ff0.p_Pa)/dT / 1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CubicEOS_g1,'YR') CubicEOS_g1='YFR';end;
%%%%%%%%%%%%%%%%%%%%%%%%% show the results  %%%%%%%%%%%%%%%%%%%%%%%%%
try
    % disp([num2str(MoleFrac1,'%0.2f'),' ',Refrigerant{1},' + ',num2str(1-MoleFrac1,'%0.2f'),' ' ,Refrigerant{2},' at ',num2str(temp_K,'%0.3f'),' K and ',num2str(pres_kPa/1e3,'%0.3f'),' MPa',]);
    if ncomp == 2
        fprintf(fid,[num2str(MassFrac(1),'%0.2f'),' ',Refrigerant{1},' + ',num2str(MassFrac(2),'%0.2f'),' ' ,Refrigerant{2},' in mass frac. with ' ,CubicEOS_g1,' EoS  \n']);
    elseif ncomp == 3
        fprintf(fid,[num2str(MassFrac(1),'%0.2f'),' ',Refrigerant{1},' + ',num2str(MassFrac(2),'%0.2f'),' ' ,Refrigerant{2},' + ',num2str(MassFrac(3),'%0.2f'),' ' ,Refrigerant{3},' in mass frac. with ' ,CubicEOS_g1,' EoS  \n']);
    else
        fprintf(fid,[Refrigerant{1},': ' ,CubicEOS_g1,' \n']);
    end

    if length(ff.Phase) == 2
        fprintf(fid,' T:             %10.6f   K  \n ',ff.T_K);
        fprintf(fid,'P:             %10.6f   kPa \n ',ff.p_Pa/1e3);
        fprintf(fid,'MM:            %10.3f   g/mol \n ',ff.MM_gmol);
        fprintf(fid,'Phase:               %4s \n\n',ff.Phase);
        fprintf(fid,'Phase behavior in mass fraction.  Vapor Frac: %8.6f\n', ff.FracV_mass);
        if ncomp == 3
            fprintf(fid,'            %10s     %10s   %10s     \n' ,Refrigerant{:}  );
            fprintf(fid,'Liquid:     %10.6f    %10.6f  %10.6f \n',ff.MassF_Li);
            fprintf(fid,'Vapor:      %10.6f    %10.6f  %10.6f \n\n',ff.MassF_Vi);
        else
            fprintf(fid,'            %10s      %10s  \n' ,Refrigerant{:}  );
            fprintf(fid,'Liquid:     %10.6f    %10.6f \n',ff.MassF_Li);
            fprintf(fid,'Vapor:      %10.6f    %10.6f \n\n',ff.MassF_Vi);
        end
        %     fprintf(fid,'Phase behavior in mole fraction.  Vapor Frac: %8.6f\n', ff.FracV_mole);
        %     fprintf(fid,'          %10s      %10s  \n' ,Refrigerant{:}  )
        %     fprintf(fid,'Liquid:        %10.6f  %10.6f  \n',ff.MoleF_Li);
        %     fprintf(fid,'Vapor:         %10.6f  %10.6f  \n\n',ff.MoleF_Vi);
        fprintf(fid,'       ---- Properties in each phase  ----       \n ');
        fprintf(fid,'                    Liquid       Vapor\n '   );
        fprintf(fid,'Z:             %10.6f  %10.6f\n ',ff.Z);
        fprintf(fid,'MM:            %10.3f  %10.3f   g/mol \n ',ff.MM_phase_gmol);
        fprintf(fid,'rho:           %10.3f  %10.3f   kg/m3 \n ',ff.rho_kgm3);
        fprintf(fid,'Entropy:       %10.3f  %10.3f   kJ/K/kg\n ',ff.ss_JkgK/1000);
        fprintf(fid,'Enthalpy:      %10.3f  %10.3f   kJ/kg\n ',ff.hh_Jkg/1000);
        fprintf(fid,'cp:            %10.3f  %10.3f   kJ/K/kg\n ',ff.cp_JkgK/1000);
        fprintf(fid,'cv:            %10.3f  %10.3f   kJ/K/kg\n ',ff.cv_JkgK/1000);
        fprintf(fid,'SoS:           %10.3f  %10.3f   m/s\n ',ff.sos_ms);
        fprintf(fid,'vis:           %10.6f  %10.6f   mPa s\n ',ff.vis_Pas*1000);
        fprintf(fid,'TC:            %10.6f  %10.6f   W/m/K\n ',ff.lambda_WmK);
        fprintf(fid,'Kappa          %10.6f  %10.6f   \n ',ff.kappa);
        fprintf(fid,'d(rho)/dP:     %10.6f  %10.6f   kg/m3/kPa \n ',ff.drhokg_dpkPa);
        fprintf(fid,'d(rho)/dT:     %10.6f  %10.6f   kg/m3/K \n ',ff.drhokg_dT);
        fprintf(fid,'dp/dT_rho:     %10.6f  %10.6f   kPa/K \n ',ff.dp_dT_kPaK);
        fprintf(fid,'\n   -- Properties with all phases combined --\n ');
        fprintf(fid,'rho_all:             %10.3f         kg/m3 \n ', ff.rho_kgm3_all  );
        fprintf(fid,'Entropy_all:         %10.3f         kJ/K/kg \n ', ff.ss_JkgK_all/1000   );
        fprintf(fid,'Enthalpy_all:        %10.3f         kJ/kg \n ', ff.hh_Jkg_all/1000   );
        fprintf(fid,'cp_all:              %10.3f         kJ/K/kg \n ', ff.cp_JkgK_all/1000   );
        fprintf(fid,'cv_all:              %10.3f         kJ/K/kg \n ', ff.cv_JkgK_all/1000   );
        fprintf(fid,'Kappa_all            %10.6f     \n ',ff.kappa_all);
    else
        %     fprintf(fid,'Fluid at %1s phase \n ',ff.Phase);
        fprintf(fid,'Fluid highly likely in the %1s phase\n ',ff.Phase);
        fprintf(fid,'T:             %10.6f   K  \n ',ff.T_K);
        fprintf(fid,'p:             %10.6f   kPa \n ',ff.p_Pa/1e3);
        fprintf(fid,'Z:             %10.6f\n ',ff.Z);
        fprintf(fid,'MM:            %10.3f   g/mol \n ',ff.MM_gmol);
        fprintf(fid,'rho_kg:        %10.3f   kg/m3 \n ',ff.rho_kgm3);
        %     fprintf(fid,'rho_mol:       %10.3f   kg/mol \n ',ff.rho_molm3);
        fprintf(fid,'d(rho)/dP:     %10.6f   kg/m3/kPa \n ',ff.drhokg_dpkPa);
        fprintf(fid,'d(rho)/dT:     %10.6f   kg/m3/K \n ',ff.drhokg_dT);
        fprintf(fid,'dp/dT_rho:     %10.6f   kPa/K \n ',ff.dp_dT_kPaK);
        fprintf(fid,'Entropy:       %10.6f   kJ/K/kg\n ',ff.ss_JkgK/1000);
        %     fprintf(fid,'Entropy:       %10.6f   J/K/mol\n ',ff.ss_JmolK);
        fprintf(fid,'Enthalpy:      %10.3f   kJ/kg\n ',ff.hh_Jkg/1000);
        %     fprintf(fid,'Enthalpy:      %10.3f   J/mol\n ',ff.hh_Jmol);
        fprintf(fid,'cp:            %10.6f   kJ/K/kg\n ',ff.cp_JkgK/1000);
        fprintf(fid,'cv:            %10.6f   kJ/K/kg\n ',ff.cv_JkgK/1000);
        fprintf(fid,'SoS:           %10.3f   m/s\n ',ff.sos_ms);
        fprintf(fid,'vis:           %10.6f   mPa s\n ',ff.vis_Pas*1000);
        fprintf(fid,'TC:            %10.6f   W/m/K\n ',ff.lambda_WmK);
        fprintf(fid,'Kappa          %10.6f   \n ',ff.kappa);
    end
catch
    fprintf(fid,'计算失败');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% try ff = OilPropm('All','P',pres_kPa,'Q',1,MassFrac,GL,0,0);  fprintf(fid,'\n V: %4s, Den: %10.3f   kg/m3',ff.Phase,ff.rho_kgm3_all); catch,end
% try ff = OilPropm('All','P',pres_kPa,'Q',0,MassFrac,GL,0,0);  fprintf(fid,'\n L: %4s, Den: %10.3f   kg/m3',ff.Phase,ff.rho_kgm3_all); catch,end
% try ff = OilPropm('All','T',temp_K,'Q',1,MassFrac,GL,0,0);         fprintf(fid,'\n V: %4s, Den: %10.3f   kg/m3',ff.Phase,ff.rho_kgm3_all); catch,end
% try ff = OilPropm('All','T',temp_K,'Q',0,MassFrac,GL,0,0);  fprintf(fid,'\n L: %4s, Den: %10.3f   kg/m3 \n',ff.Phase,ff.rho_kgm3_all); catch,end
%

try
    if ncomp == 1
        D_ref = refpropm('D',input1,para1,input2,para2,Refrigerant{:});
        Cp_ref = refpropm('C',input1,para1,input2,para2,Refrigerant{:});
        Cv_ref = refpropm('O',input1,para1,input2,para2,Refrigerant{:});
        sos_ref = refpropm('A',input1,para1,input2,para2,Refrigerant{:});
        vis_ref = refpropm('V',input1,para1,input2,para2,Refrigerant{:});
        tc_ref = refpropm('L',input1,para1,input2,para2,Refrigerant{:});
        drhodP = refpropm('R',input1,para1,input2,para2,Refrigerant{:});
        drhodT = refpropm('W',input1,para1,input2,para2,Refrigerant{:});
        dpdT = refpropm('#',input1,para1,input2,para2,Refrigerant{:});

    elseif ncomp == 2
        D_ref = refpropm('D',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        Cp_ref = refpropm('C',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        Cv_ref = refpropm('C',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        sos_ref = refpropm('A',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        vis_ref = refpropm('V',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        tc_ref = refpropm('L',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        drhodP = refpropm('R',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        drhodT = refpropm('W',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
        dpdT = refpropm('#',input1,para1,input2,para2,Refrigerant{1},Refrigerant{2},MassFrac);
    elseif ncomp == 3
        D_ref = refpropm('D','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        Cp_ref = refpropm('C','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        Cv_ref = refpropm('C','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        sos_ref = refpropm('A','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        vis_ref = refpropm('V','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        tc_ref = refpropm('L','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        drhodP = refpropm('R','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        drhodT = refpropm('W','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
        dpdT = refpropm('#','T',temp_K,'P',pres_kPa,Refrigerant{1},Refrigerant{2},Refrigerant{3},MassFrac);
    end
    fprintf(fid,'\n\n----- Refprop Calculations --- \n ');
    fprintf(fid,'rho_kg:        %10.3f   kg/m3 \n ',D_ref);
    fprintf(fid,'cp:            %10.6f   kJ/K/kg\n ',Cp_ref/1000);
    fprintf(fid,'cv:            %10.6f   kJ/K/kg\n ',Cv_ref/1000);
    fprintf(fid,'SoS:           %10.3f   m/s\n ',sos_ref);
    %         fprintf(fid,'sr:            %10.3f   J/(kg K)\n ',sr);
    fprintf(fid,'vis:           %10.6f   mPa s\n ',vis_ref*1000);
    fprintf(fid,'TC:            %10.6f   W/m/K\n ',tc_ref);
    fprintf(fid,'d(rho)/dP:     %10.6f   kg/m3/kPa \n ',drhodP);
    fprintf(fid,'d(rho)/dT:     %10.6f   kg/m3/K \n ',drhodT);
    fprintf(fid,'dp/dT:         %10.6f   kPa/K \n ',dpdT);
    fprintf(fid,'------------------------------- \n ');
catch
    fprintf(fid,'\n\n--- Refprop Calculation fail --=  \n ');
end
fclose(fid);
text2put= fileread('tmp.txt');
delete('tmp.txt');
end