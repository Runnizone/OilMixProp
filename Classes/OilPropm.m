function ff = OilPropm(OutPut,InPut_1,Value_1,InPut_2,Value_2,MassFrac,GL,T_K_guess,p_kPa_guess) %#codegen
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
%     All all the following properties               [a class]
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

    OutPut = upper(OutPut);
    [mm_mix_gmol,Zi] = EOSmodel.MassF_2_MoleF(GL.MM_gmol,MassFrac);  
    ncomp = length(Zi);

    if strcmpi(InPut_1,'T')
        if contains(OutPut,'T'), ff = Value_1; return; end
        if strcmpi(InPut_2,'P')
            T_K = Value_1;
            p_Pa = Value_2 * 1000;
            if contains(OutPut,'P'), ff = p_Pa/1000; return; end
        elseif strcmpi(InPut_2,'D')
            T_K = Value_1;
            v_molar_all = 1/Value_2 * mm_mix_gmol / 1000;   % might be in two phase region
            if contains(OutPut,'+')
                p_Pa = EOSmodel.f_tv2p_1phase(GL.CubicEOS,T_K,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            else
                p_Pa = EOSmodel.f_tv2p(GL.CubicEOS,T_K,p_kPa_guess*1000,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
            if contains(OutPut,'P'), ff = p_Pa/1000; return; end
        elseif strcmpi(InPut_2,'Q')
            T_K = Value_1;
            FracV_mass = Value_2;
            if FracV_mass < 0 || FracV_mass > 1, error('Should be: 0 <= Q <=1' );  end
            if ncomp == 1
                if T_K > GL.tempc
                    error('Input temperature higher than critical temperature')
                elseif T_K < GL.temp_r
                    error('Input temperature lower than triple temperature')
                end
                [p_Pa,v_molar_L,v_molar_V] = EOSmodel.p_LVE_pure(GL.CubicEOS, T_K, GL.tempc, GL.presc, GL.Zc, GL.acentric, GL.R);
                if FracV_mass == 0
                    v_molar = v_molar_L;
                elseif FracV_mass == 1
                    v_molar = v_molar_V;
                else
                    error('Pure fluids: set Q = 0 for saturated liquid, then Q = 1 for saturate vapor, respectively')
                end
            else
                [p_Pa, MoleF_Li, MoleF_Vi, ~, ~, v_molar] = ...
                    EOSmodel.f_tq2pv(GL.CubicEOS, FracV_mass, T_K,Zi, p_kPa_guess*1000 ,GL.MM_gmol ,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
        elseif strcmpi(InPut_2,'H')
            T_K = Value_1;
            hh_Jmol = Value_2 * mm_mix_gmol / 1000;
            if contains(OutPut,'+')
                [p_Pa,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_th2pv_1phase(GL.CubicEOS,hh_Jmol,T_K,Zi,p_kPa_guess*1000,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            else
                [p_Pa,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_th2pv(GL.CubicEOS,hh_Jmol,T_K,Zi,p_kPa_guess*1000,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            end
        elseif strcmpi(InPut_2,'S')
            T_K = Value_1;
            ss_JmolK = Value_2 * mm_mix_gmol / 1000;
            if contains(OutPut,'+')
                [p_Pa,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_ts2pv_1phase(GL.CubicEOS,ss_JmolK,T_K,Zi,p_kPa_guess*1000,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            else
                [p_Pa,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_ts2pv(GL.CubicEOS,ss_JmolK,T_K,Zi,p_kPa_guess*1000,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            end
        else
            error('When 1st input is T, the 2nd input can only be P, Q, D, H, S')
        end
    elseif strcmpi(InPut_1,'P')
        if contains(OutPut,'P'), ff = Value_1; return; end
        if strcmpi(InPut_2,'T')
            p_Pa = Value_1 * 1000;
            T_K = Value_2;
            if contains(OutPut,'T'), ff = T_K; return; end
        elseif strcmpi(InPut_2,'D')
            p_Pa = Value_1 * 1000;
            v_molar_all = 1/Value_2 * mm_mix_gmol / 1000;   % might be in two phase region
            if contains(OutPut,'+')
                T_K = EOSmodel.f_pv2t_1phase(GL.CubicEOS,T_K_guess,p_Pa,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            else
                T_K = EOSmodel.f_pv2t(GL.CubicEOS,T_K_guess,p_Pa,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
            if contains(OutPut,'T'), ff = T_K; return; end
        elseif strcmpi(InPut_2,'Q')
            p_Pa = Value_1 * 1000;
            FracV_mass = Value_2;
            if FracV_mass < 0 || FracV_mass > 1, error('Should be: 0 <= Q <=1' );  end
            if ncomp == 1
                if p_Pa > GL.presc
                    error('Input pressure higher than critical pressure')
                elseif p_Pa < GL.pres_r
                    error('Input pressure lower than triple pressure')
                end
                [T_K,v_molar_L,v_molar_V] = EOSmodel.T_LVE_pure(GL.CubicEOS, p_Pa, GL.tempc, GL.presc, GL.Zc, GL.acentric, GL.R);
                if FracV_mass == 0
                    v_molar = v_molar_L;
                elseif FracV_mass == 1
                    v_molar = v_molar_V;
                else
                    error('Pure fluids: set Q = 0 for saturated liquid, then Q = 1 for saturate vapor, respectively');
                end
            else
                if p_Pa < max(GL.pres_r)
                    error('Input pressure lower than triple pressure')
                end
                [T_K, MoleF_Li, MoleF_Vi, ~, ~,v_molar] = ...
                    EOSmodel.f_pq2tv(GL.CubicEOS, FracV_mass, p_Pa,Zi, T_K_guess ,GL.MM_gmol,GL.temp_r, GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
        elseif strcmpi(InPut_2,'H')
            p_Pa = Value_1 * 1000;
            hh_Jmol = Value_2 * mm_mix_gmol / 1000;
            if contains(OutPut,'+')
                [T_K,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_ph2tv_1phase(GL.CubicEOS,hh_Jmol,p_Pa,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            else
                [T_K,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_ph2tv(GL.CubicEOS,hh_Jmol,p_Pa,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            end
        elseif strcmpi(InPut_2,'S')
            p_Pa = Value_1 * 1000;
            ss_JmolK = Value_2 * mm_mix_gmol / 1000;
            if contains(OutPut,'+')
                [T_K,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_ps2tv_1phase(GL.CubicEOS,ss_JmolK,p_Pa,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            else
                [T_K,v_molar,FracV_mole, MoleF_Li, MoleF_Vi] = ...
                    EOSmodel.f_ps2tv(GL.CubicEOS,ss_JmolK,p_Pa,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
            end
        else
            error('When 1st input is P, the 2nd input can only be T Q D H S')
        end
    else
        error('The first input of OilPropm can only be T or p')
    end

    if strcmpi(InPut_2,'Q') || strcmpi(InPut_2,'S') || strcmpi(InPut_2,'H') 
        if contains(OutPut,'T'), ff = T_K; return; end
        if contains(OutPut,'P'), ff = p_Pa/1000; return; end
        if contains(OutPut,'Q') 
           if length(v_molar) == 1
                Z = p_Pa/GL.R/T_K*v_molar; 
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
        if contains(OutPut,'X')
            if length(v_molar) == 1
                ff = MassFrac; return; 
            else
                [~,MassF_Li] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Li);  
                [~,MassF_Vi] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Vi);  
                ff = [MassF_Li,MassF_Vi]; 
                return; 
            end
        end
        if contains(OutPut,'D')
            if length(v_molar) == 1
                ff = mm_mix_gmol / 1000 / v_molar; return; 
            else
                [mm_L_gmol,~] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Li);  
                [mm_V_gmol,~] = EOSmodel.MoleF_2_MassF(GL.MM_gmol,MoleF_Vi);  
                if ~exist('FracV_mass','var'), FracV_mass = FracV_mole*mm_V_gmol / ((1-FracV_mole) * mm_L_gmol + FracV_mole*mm_V_gmol); end
                ff = [mm_L_gmol / 1000 / v_molar(1), mm_V_gmol / 1000 / v_molar(2)]; 
                ff = [1 / ((1-FracV_mass)/ff(1) + FracV_mass/ff(2)), ff];
                return; 
            end
        end
        if contains(OutPut,'Z'), ff = p_Pa/GL.R/T_K.*v_molar; return; end
    end


    %% call OilProp
    if ~exist('v_molar','var'), v_molar = 0; end
    ff = OilProp(OutPut,GL,T_K,p_Pa,v_molar,Zi);
    ff.Phase = replace(ff.Phase,'X','L');
  
    %% prepare output
    if contains(OutPut,'D') 
        if ff.nphase == 1, ff = ff.rho_kgm3_all; else, ff = [ff.rho_kgm3_all,ff.rho_kgm3]; end
    elseif contains(OutPut,'X')
        if ff.nphase == 1, ff = ff.MassF_Zi; else, ff = [ff.MassF_Li,ff.MassF_Vi]; end
    elseif contains(OutPut,'Z')  
        ff = ff.Z; 
    elseif contains(OutPut,'Q')  
        ff = ff.FracV_mass;
    elseif contains(OutPut,'W')
        ff = ff.drhokg_dT; 
    elseif contains(OutPut,'R')
        ff = ff.drhokg_dpkPa; 
    elseif contains(OutPut,'#')
        ff = ff.dp_dT_kPaK; 
    elseif contains(OutPut,'O')
        if ff.nphase == 1, ff = ff.cv_JkgK_all; else, ff = [ff.cv_JkgK_all,ff.cv_JkgK]; end
    elseif contains(OutPut,'C')
        if ff.nphase == 1, ff = ff.cp_JkgK_all; else, ff = [ff.cp_JkgK_all,ff.cp_JkgK]; end
    elseif contains(OutPut,'K')
        ff = ff.kappa; 
    elseif contains(OutPut,'S')
        if ff.nphase == 1, ff = ff.ss_JkgK_all; else, ff = [ff.ss_JkgK_all,ff.ss_JkgK]; end
    elseif contains(OutPut,'H')
        if ff.nphase == 1, ff = ff.hh_Jkg_all; else, ff = [ff.hh_Jkg_all,ff.hh_Jkg]; end
    elseif contains(OutPut,'U')
        if ff.nphase == 1, ff = ff.uu_Jkg_all; else, ff = [ff.uu_Jkg_all,ff.uu_Jkg]; end
    elseif contains(OutPut,'A') && ~contains(OutPut,'ALL')
        ff = ff.sos_ms; 
    elseif contains(OutPut,'V')
        ff = ff.vis_Pas; 
    elseif contains(OutPut,'L') && ~contains(OutPut,'ALL')
        ff = ff.lambda_WmK; 
    elseif contains(OutPut,'ALL')
        return;
    else
        warning('Unknown OutPut, Please specify an OutPut (D, H, etc) to speed up the calculation.')
    end
end