function ff = OilPropm(OutPut,InPut_1,Value_1,InPut_2,Value_2,MassFrac,GL,T_K_guess,p_kPa_guess)
% Equation Package for OilMixProp
% Developed by Xiaoxian Yang (xiaoxian.yang@mb.tu-chemnitz.de or runnizone@qq.com)
% Package is free to use; bugs exist; the developer is no responsible for the
% potential errors and the consequence due to them. 
% the unit system: K, mol, kg, J, m, kPa (molar mass: g/mol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%  Example to call this function %%%%%%%%%%%%%%%
%%% see FluidCalc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Avaiable input sets
% First input parameter:
% P   Pressure [kPa]
% T   Temperature [K]
%
% Second input parameter:
% P   Pressure [kPa]
% D   Density [kg/m3]
% Q   Quality (vapor fraction) (kg/kg)
% H   Enthalpy [J/kg]
% S   Entropy [J/kg/K]
% T   Temperature [K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Avaiable OutPut sets
%     A   Speed of sound [m/s]                       [liquid phase, gas phase]
%     C   Cp [J/(kg K)]                              [liquid phase, gas phase, all phase]
%     D   Density [kg/m^3]                           [liquid phase, gas phase, all phase]
%     H   Enthalpy [J/kg]                            [liquid phase, gas phase, all phase]
%     K   Ratio of specific heats (Cp/Cv) [-]        [liquid phase, gas phase, all phase]
%     L   Thermal conductivity [W/(m K)]             [liquid phase, gas phase]
%     O   Cv [J/(kg K)]                              [liquid phase, gas phase, all phase]
%     P   Pressure [kPa]                              
%     Q   Quality (vapor fraction) (kg/kg)        
%     S   Entropy [J/(kg K)]                         [liquid phase, gas phase, all phase]
%     T   Temperature [K]                            
%     U   Internal energy [J/kg]                     [liquid phase, gas phase, all phase]
%     V   Dynamic viscosity [Pa*s]                   [liquid phase, gas phase]
%     X   Liquid phase & gas phase comp.(mass frac)  [1st column, liquid phase; 2nd column, gas phase]
%     Z   Compressibility factor                     [liquid phase, gas phase]
%     #   dP/dT     (constant rho) [kPa/K]           [liquid phase, gas phase]
%     R   d(rho)/dP (constant T)   [kg/m^3/kPa]      [liquid phase, gas phase]
%     W   d(rho)/dT (constant p)   [kg/(m^3 K)]      [liquid phase, gas phase]
% to be developped
%     F   Fugacity [kPa] (returned as an array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 7 
        T_K_guess = 0;
        p_kPa_guess = 0;
    elseif nargin == 8
        p_kPa_guess = 0;
    end

    [mm_mix_gmol,Zi] = EOSmodel.MassF_2_MoleF(GL.MM_gmol,MassFrac);  
    ncomp = length(Zi);

    if strcmpi(InPut_1,'T')
        if strcmpi(InPut_2,'p')
            temp = Value_1;
            pres_Pa = Value_2 * 1000;
        elseif strcmpi(InPut_2,'D')
            temp = Value_1;
            v_molar_all = 1/Value_2 * mm_mix_gmol / 1000;   % might be in two phase region
            pres_Pa = EOSmodel.f_tv2p(GL.CubicEOS,temp,p_kPa_guess*1000,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            if strcmpi(OutPut,'P'), ff = pres_Pa/1000; return; end
        elseif strcmpi(InPut_2,'Q')
            temp = Value_1;
            FracV_mass = Value_2;
            if FracV_mass < 0 || FracV_mass > 1, error('Should be: 0 <= Q <=1' );  end
            if ncomp == 1
                if temp > GL.tempc
                    error('Input temperature higher than critical temperature')
                elseif temp < GL.temp_r
                    error('Input temperature lower than triple temperature')
                end
                [pres_Pa,v_molar_L,v_molar_V] = EOSmodel.p_LVE_pure(GL.CubicEOS, temp, GL.tempc, GL.presc, GL.Zc, GL.acentric, GL.R);
                if FracV_mass == 0
                    v_molar = v_molar_L;
                elseif FracV_mass == 1
                    v_molar = v_molar_V;
                else
                    error('Pure fluids: set Q = 0 for saturated liquid, then Q = 1 for saturate vapor, respectively')
                end
            else
                [pres_Pa, MoleF_Li, MoleF_Vi, ~, ~, v_molar] = ...
                    EOSmodel.f_tq2pv(GL.CubicEOS, FracV_mass, temp,Zi, p_kPa_guess*1000 ,GL.MM_gmol ,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
        elseif strcmpi(InPut_2,'H')
            temp = Value_1;
            hh_Jmol = Value_2 * mm_mix_gmol / 1000;
            [pres_Pa,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                EOSmodel.f_th2pv(GL.CubicEOS,hh_Jmol,temp,Zi,p_kPa_guess*1000,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        elseif strcmpi(InPut_2,'S')
            temp = Value_1;
            ss_JmolK = Value_2 * mm_mix_gmol / 1000;
            [pres_Pa,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                EOSmodel.f_ts2pv(GL.CubicEOS,ss_JmolK,temp,Zi,p_kPa_guess*1000,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        else
            error(['When 1st input is T, the 2nd input can only be P or Q'])
        end
    elseif strcmpi(InPut_1,'P')
        if strcmpi(InPut_2,'T')
            pres_Pa = Value_1 * 1000;
            temp = Value_2;
        elseif strcmpi(InPut_2,'D')
            pres_Pa = Value_1 * 1000;
            v_molar_all = 1/Value_2 * mm_mix_gmol / 1000;   % might be in two phase region
            temp = EOSmodel.f_pv2t(GL.CubicEOS,T_K_guess,pres_Pa,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
        elseif strcmpi(InPut_2,'Q')
            pres_Pa = Value_1 * 1000;
            FracV_mass = Value_2;
            if FracV_mass < 0 || FracV_mass > 1, error('Should be: 0 <= Q <=1' );  end
            if ncomp == 1
                if pres_Pa > GL.presc
                    error('Input pressure higher than critical pressure')
                elseif pres_Pa < GL.pres_r
                    error('Input pressure lower than triple pressure')
                end
                [temp,v_molar_L,v_molar_V] = EOSmodel.T_LVE_pure(GL.CubicEOS, pres_Pa, GL.tempc, GL.presc, GL.Zc, GL.acentric, GL.R);
                if FracV_mass == 0
                    v_molar = v_molar_L;
                elseif FracV_mass == 1
                    v_molar = v_molar_V;
                else
                    error('Pure fluids: set Q = 0 for saturated liquid, then Q = 1 for saturate vapor, respectively');
                end
            else
                if pres_Pa < max(GL.pres_r)
                    error('Input pressure lower than triple pressure')
                end
                [temp, MoleF_Li, MoleF_Vi, ~, ~,v_molar] = ...
                    EOSmodel.f_pq2tv(GL.CubicEOS, FracV_mass, pres_Pa,Zi, T_K_guess ,GL.MM_gmol,GL.temp_r, GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end

        elseif strcmpi(InPut_2,'H')
            pres_Pa = Value_1 * 1000;
            hh_Jmol = Value_2 * mm_mix_gmol / 1000;
            [temp,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                EOSmodel.f_ph2tv(GL.CubicEOS,hh_Jmol,pres_Pa,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        elseif strcmpi(InPut_2,'S')
            pres_Pa = Value_1 * 1000;
            ss_JmolK = Value_2 * mm_mix_gmol / 1000;
            [temp,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                EOSmodel.f_ps2tv(GL.CubicEOS,ss_JmolK,pres_Pa,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        else
            error(['When 1st input is T, the 2nd input can only be T Q H S'])
        end
    else
        error('The first input of OilPropm can only be T or p')
    end

    if strcmpi(InPut_2,'Q') || strcmpi(InPut_2,'S') || strcmpi(InPut_2,'H') 
        if strcmpi(OutPut,'T'), ff = temp; return; end
        if strcmpi(OutPut,'P'), ff = pres_Pa/1000; return; end
        if strcmpi(OutPut,'Q') 
           if length(v_molar) == 1
                Z = pres_Pa/GL.R/temp*v_molar; 
                if Z > 0.3, ff = 1; else, ff = 0; end
                return; 
            else
                [mm_L_gmol,~] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Li);  
                [mm_V_gmol,~] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Vi);  
                if ~strcmpi(InPut_2,'Q') 
                    ff = FracV_mole*mm_V_gmol / ((1-FracV_mole) * mm_L_gmol + FracV_mole*mm_V_gmol);
                else
                    ff = FracV_mass;
                end
                return; 
           end
        end
        if strcmpi(OutPut,'X')
            if length(v_molar) == 1
                ff = MassFrac; return; 
            else
                [~,MassF_Li] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Li);  
                [~,MassF_Vi] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Vi);  
                ff = [MassF_Li,MassF_Vi]; 
                return; 
            end
        end
        if strcmpi(OutPut,'D')
            if length(v_molar) == 1
                ff = mm_mix_gmol / 1000 / v_molar; return; 
            else
                [mm_L_gmol,~] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Li);  
                [mm_V_gmol,~] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Vi);  
                if ~exist('FracV_mass','var'), FracV_mass = FracV_mole*mm_V_gmol / ((1-FracV_mole) * mm_L_gmol + FracV_mole*mm_V_gmol); end
                ff = [mm_L_gmol / 1000 / v_molar(1), mm_V_gmol / 1000 / v_molar(2)]; 
                ff = [ff,1 / ((1-FracV_mass)/ff(1) + FracV_mass/ff(2))];
                return; 
            end
        end
        if strcmpi(OutPut,'Z'), ff = pres_Pa/GL.R/temp.*v_molar; return; end
    end


    %% call OilProp
    if ~exist('v_molar','var'), v_molar = 0; end
    ff = OilProp(OutPut,GL,temp,pres_Pa,v_molar,Zi);

    %% prepare output
    if strcmpi(OutPut,'D') 
        if ff.nphase == 1, ff = ff.rho_kgm3_all; else, ff = [ff.rho_kgm3,ff.rho_kgm3_all]; end
    elseif strcmpi(OutPut,'X')
        if ff.nphase == 1, ff = ff.MassF_Zi; else, ff = [ff.MassF_Li,ff.MassF_Vi]; end
    elseif strcmpi(OutPut,'Z')  
        ff = ff.Z; 
    elseif strcmpi(OutPut,'Q')  
        ff = ff.FracV_mass;
    elseif strcmpi(OutPut,'W')
        ff = ff.drhokg_dT; 
    elseif strcmpi(OutPut,'R')
        ff = ff.drhokg_dpkPa; 
    elseif strcmpi(OutPut,'#')
        ff = ff.dp_dT_kPaK; 
    elseif strcmpi(OutPut,'O')
        if ff.nphase == 1, ff = ff.cv_JkgK_all; else, ff = [ff.cv_JkgK,ff.cv_JkgK_all]; end
    elseif strcmpi(OutPut,'C')
        if ff.nphase == 1, ff = ff.cp_JkgK_all; else, ff = [ff.cp_JkgK,ff.cp_JkgK_all]; end
    elseif strcmpi(OutPut,'K')
        ff = ff.kappa; 
    elseif strcmpi(OutPut,'S')
        if ff.nphase == 1, ff = ff.ss_JkgK_all; else, ff = [ff.ss_JkgK,ff.ss_JkgK_all]; end
    elseif strcmpi(OutPut,'H')
        if ff.nphase == 1, ff = ff.hh_Jkg_all; else, ff = [ff.hh_Jkg,ff.hh_Jkg_all]; end
    elseif strcmpi(OutPut,'U')
        if ff.nphase == 1, ff = ff.uu_Jkg_all; else, ff = [ff.uu_Jkg,ff.uu_Jkg_all]; end
    elseif strcmpi(OutPut,'A')
        ff = ff.sos_ms; 
    elseif strcmpi(OutPut,'V')
        ff = ff.vis_Pas; 
    elseif strcmpi(OutPut,'L')
        ff = ff.lambda_WmK; 
    elseif strcmpi(OutPut,'All')
        return;
    else
        warning('Unknown OutPut, Please specify an OutPut (D, H, etc) to speed up the calculation.')
    end
end