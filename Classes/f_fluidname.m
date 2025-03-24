function  f_fluidname(INPUT)
% different type of inputs 
% All: all fluids
% REFPROP: all fluids avaiable in REFPROP 10.0
% Others: all fluids except for those in REFPROP 10.0
% oils: all oils which are fitted using OilMixProp
% CASID: the fluid name with the given CASID

    %% fluid constants
    filename_ref = ['Classes\Fluid_Constants.txt'];
    [Num,Material_ref,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,CASID_ref] = ...
        textread(filename_ref,'%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','headerlines',1,'delimiter',';');
    N_ref = length(Num);

    %% fluid constants more than NIST
    filename_e = ['Classes\Fluid_Constants_ext.txt'];
    [Num_e,Material_e,~,CASID_e,~,~,~,~,~,~,~,~,~] = textread(filename_e,'%f%s%s%s%f%f%f%f%f%f%f%f%f','headerlines',1,'delimiter',';');
    N_e = length(Num_e);


    %% fluid constants for oils
    filename_Oil = ['Classes\Fluid_Constants_Fitted.txt'];
    [Material_Oil,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        textread('Classes/Fluid_Constants_Fitted.txt','%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'delimiter',';');         
    N_Oil = length(Material_Oil);

    Lprint = 0;
    icount = 0;

    if strcmpi(INPUT,'All') || strcmpi(INPUT,'REFPROP')
        disp('************** Fluids in REFPROP 10.0 **************')
        fprintf('No.                Name              CAS RN\n');
        for i = 1:N_ref
            icount = icount + 1;
            fprintf('%3i%20s%20s\n',icount,Material_ref{i},CASID_ref{i});
        end
        Lprint = 1;
    end

    if strcmpi(INPUT,'All') || strcmpi(INPUT,'Others') || strcmpi(INPUT,'Other')
        disp('***************************** More fluids *****************************')
        fprintf('  No.                                    Name         CAS RN\n');
        for i = 1:N_e
            icount = icount + 1;
            fprintf('%5i%40s%15s\n',icount,Material_e{i},CASID_e{i});
        end
        Lprint = 1;
    end

    if strcmpi(INPUT,'All') || strcmpi(INPUT,'Oils') || strcmpi(INPUT,'Oil')
        disp('************** Oils in OilMixProp 1.0 **************')
        fprintf('No.                Name              CAS RN\n');
        for i = 1:N_Oil
            icount = icount + 1;
            fprintf('%3i%20s%20s\n',icount,Material_Oil{i},'No Data');
        end
        Lprint = 1;
    end

    if Lprint == 0
        INPUT = regexprep(INPUT,'-','');
        for i = 1:N_ref
            CASID = regexprep(CASID_ref{i},'-','');  
            if strcmpi(INPUT,CASID)
                icount = icount + 1;
                fprintf('No.                Name              CAS RN\n');
                fprintf('%3i%20s%20s\n',icount,Material_ref{i},CASID_ref{i});
            end
        end
        Lprint = 1;
    end

    if Lprint == 0
        disp('invalid input!!!!')
        disp('Here are valid inputs:');
        disp('All: all fluids');
        disp('REFPROP: all fluids avaiable in REFPROP 10.0');
        disp('Others: all fluids except for those in REFPROP 10.0');
        disp('oils: all oils which are fitted using OilMixProp');
        disp('CASID: the fluid name with the given CASID');
    end


end

