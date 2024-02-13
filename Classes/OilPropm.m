function ff = OilPropm(parameter_1,value_1,parameter_2,value_2,MassFrac,GL,T_K_guess,p_pa_guess)
% Equation Package for OilProp
% Developed by Xiaoxian Yang (xiaoxian.yang@mb.tu-chemnitz.de or runnizone@qq.com)
% Package is free to use; bugs exist; the developer is no responsible for the
% potential errors and the consequence due to them. 
% the unit system: K, mol, kg, J, m, Pa (molar mass: g/mol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%  Example to call this function %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
%%% see FluidCalc.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% call the function %%%%%%%%%%%%%%%%%%%%%%%%%
% ff1 = OilPropm('P',pres,'T',temp,Zi,GL,T_K_guess,p_Pa_guess);
% % ff1 = OilPropm('T',temp,'Q',0,Zi,GL,T_K_guess,p_Pa_guess);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Avaiable input sets
% First input parameter:
% P   Pressure [Pa]
% T   Temperature [K]
%
% Second input parameter:
% P   Pressure [Pa]
% D   Density [kg/m3]
% Q   Quality (vapor fraction) (kg/kg)
% H   Enthalpy [J/kg]
% S   Entropy [J/kg/K]
% T   Temperature [K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 6 
        T_K_guess = 0;
        p_pa_guess = 0;
    elseif nargin == 7
        p_pa_guess = 0;
    end

    [mm_mix_gmol,Zi] = EOSmodel.MassF_2_MoleF(GL.MM_gmol,MassFrac);  
    ncomp = length(Zi);

    if strcmpi(parameter_1,'T')
        if strcmpi(parameter_2,'p')
            temp = value_1;
            pres = value_2;
        elseif strcmpi(parameter_2,'D')
            temp = value_1;
            v_molar_all = 1/value_2 * mm_mix_gmol / 1000;   % might be in two phase region
            pres = EOSmodel.f_tv2p(GL.CubicEOS,temp,p_pa_guess,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
        elseif strcmpi(parameter_2,'Q')
            temp = value_1;
            FracV = value_2;
            if FracV < 0 || FracV > 1, error('Should be: 0 <= Q <=1' );  end
            if ncomp == 1
                if temp > GL.tempc
                    error('Input temperature higher than critical temperature')
                elseif temp < GL.temp_r
                    error('Input temperature lower than triple temperature')
                end
                [pres,v_molar_L,v_molar_V] = EOSmodel.p_LVE_pure(GL.CubicEOS, temp, GL.tempc, GL.presc, GL.Zc, GL.acentric, GL.R);
                if FracV == 0
                    v_molar = v_molar_L;
                elseif FracV == 1
                    v_molar = v_molar_V;
                else
                    error('Pure fluids: set Q = 0 for saturated liquid, then Q = 1 for saturate vapor, respectively')
                end
            else
%                 [pres, ~, ~, ~, ~, v_molar] = EOSmodel.f_tq2pv(GL.CubicEOS, FracV, temp,Zi, p_pa_guess ,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
                [pres, ~, ~, ~, ~, v_molar] = EOSmodel.f_tq2pv(GL.CubicEOS, FracV, temp,Zi, p_pa_guess ,GL.MM_gmol ,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
        elseif strcmpi(parameter_2,'H')
            temp = value_1;
            hh_Jmol = value_2 * mm_mix_gmol / 1000;
            [pres,v_molar] = EOSmodel.f_th2pv(GL.CubicEOS,hh_Jmol,temp,Zi,p_pa_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        elseif strcmpi(parameter_2,'S')
            temp = value_1;
            ss_JmolK = value_2 * mm_mix_gmol / 1000;
            [pres,v_molar] = EOSmodel.f_ts2pv(GL.CubicEOS,ss_JmolK,temp,Zi,p_pa_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        else
            error(['When 1st input is T, the 2nd input can only be P or Q'])
        end
    elseif strcmpi(parameter_1,'P')
        if strcmpi(parameter_2,'T')
            pres = value_1;
            temp = value_2;
        elseif strcmpi(parameter_2,'D')
            pres = value_1;
            v_molar_all = 1/value_2 * mm_mix_gmol / 1000;   % might be in two phase region
            temp = EOSmodel.f_pv2t(GL.CubicEOS,T_K_guess,pres,v_molar_all,Zi,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
        elseif strcmpi(parameter_2,'Q')
            pres = value_1;
            FracV = value_2;
            if FracV < 0 || FracV > 1, error('Should be: 0 <= Q <=1' );  end
            if ncomp == 1
                if pres > GL.presc
                    error('Input pressure higher than critical pressure')
                elseif pres < GL.pres_r
                    error('Input pressure lower than triple pressure')
                end
                [temp,v_molar_L,v_molar_V] = EOSmodel.T_LVE_pure(GL.CubicEOS, pres, GL.tempc, GL.presc, GL.Zc, GL.acentric, GL.R);
                if FracV == 0
                    v_molar = v_molar_L;
                elseif FracV == 1
                    v_molar = v_molar_V;
                else
                    error('Pure fluids: set Q = 0 for saturated liquid, then Q = 1 for saturate vapor, respectively');
                end
            else
                if pres < max(GL.pres_r)
                    error('Input pressure lower than triple pressure')
                end
%                 [temp, ~, ~, ~, ~,v_molar] = EOSmodel.f_pq2tv(GL.CubicEOS, FracV, pres,Zi, T_K_guess ,GL.temp_r, GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
                [temp, ~, ~, ~, ~,v_molar] = EOSmodel.f_pq2tv(GL.CubicEOS, FracV, pres,Zi, T_K_guess ,GL.MM_gmol,GL.temp_r, GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.R);
            end
        elseif strcmpi(parameter_2,'H')
            pres = value_1;
            hh_Jmol = value_2 * mm_mix_gmol / 1000;
            [temp,v_molar] = EOSmodel.f_ph2tv(GL.CubicEOS,hh_Jmol,pres,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        elseif strcmpi(parameter_2,'S')
            pres = value_1;
            ss_JmolK = value_2 * mm_mix_gmol / 1000;
            [temp,v_molar] = EOSmodel.f_ps2tv(GL.CubicEOS,ss_JmolK,pres,Zi,T_K_guess,GL.tempc,GL.presc,GL.Zc,GL.acentric,GL.kij,GL.k0,GL.k1,GL.T_ref,GL.v_ref,GL.R);
        else
            error(['When 1st input is T, the 2nd input can only be T Q H S'])
        end
    else
        error('The first input of OilPropm can only be T or p')
    end
    if ~exist('v_molar','var'), v_molar = 0; end
    ff = OilProp(GL,temp,pres,v_molar,Zi);
end