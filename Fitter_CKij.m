%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']);  options = optimset('Display',  'off');   warning('off');  
SSS = dbstack();  thisfile = SSS(1).file;  LL = length(thisfile);   thisfilename = thisfile(1:LL-2);
linewidth = 1; fontsize = 10;  markersize = 4;
AllEOS = {'PR','SRK','PTV','YFR'};     CubicEOS = AllEOS{4}; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the fluid to be fitted
OilNum = 2;            % oil to be study, see the fluid definition below
Lplot = 1;             % plot and save the figure? 1 yes, 0 not
Lsave2database = 1;    % save parameters to database (Classes/Bin_kij_fit.txt)? 1 yes, 0 not
% Oil mixture details
if OilNum == 1
    Material = {'Emkarate RL32','R1233zde'};     % make sure oil is the second component
    FitIndex_VP = [1,6,9];                        % Index of data used for fit, only one point needed
elseif OilNum == 2
    Material = {'PEC8','CO2'};     % make sure oil is the second component
    FitIndex_VP = [6];                        % Index of data used for fit, only one point needed
elseif OilNum == 3

else
    error(['not defined yet'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% parameter preperation
GL = GetGlobals(CubicEOS,Material);  % obtain fluid constants
Material_mix = [Material{1},'_',Material{2}];
VPfile = ['ExpData/',Material_mix,'/VaporPressure.txt'];
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
    n0 = [-0.03];    lb = [-1];    ub = [1];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
if Lplot && LvpExist == 1
    figure(1);clf; hold on; box on;
    plot(T_K_all,ReDev_all*100,'o','color','k','markersize',markersize,'linewidth',linewidth-0.4); 
    plot(T_K,ReDev*100,'o','markerfacecolor','r','markeredgecolor','r','markersize',markersize,'linewidth',linewidth-0.4); 
    plot([min(T_K_all)-10,max(T_K_all)+10],[0,0],'k-');   
    text(min(T_K_all)-5,3,['\itk\rm_i_j = ',num2str(kij_fit)])
    % axis([min(T_K_all),max(T_K_all)+10,-2.0,2.0])
    ylabel(['100\cdot(\itp\rm_b_u_b_,_e_x_p\it',char(hex2dec('2212')),'\itp\rm_b_u_b_,_f_i_t)/\itp\rm_b_u_b_,_f_i_t'],'rotation',90,'fontsize',fontsize,'FontName','Arial','linewidth',linewidth); 
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition',[0 0 18 8]);
    print(gcf,'-dtiff','-r200',['Figures/',thisfilename,'_',Material_mix,'.tiff']); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the parameters to database
if Lsave2database
    copyfile Classes/Bin_kij_fit.txt 'Classes/Bin_kij_fit_2.txt'
%     try
        [Material_1,Material_2,kij_all] = textread('Classes/Bin_kij_fit_2.txt','%s%s%f','headerlines',1,'delimiter',';');
        nM = length(Material_1);
        fid_out = fopen(['./Classes/Bin_kij_fit.txt'],'w');
        fprintf(fid_out,['          Material_1;          Material_2; Cubic_kij\n']);
        Lfound = 0;
        if nM == 0
            fprintf(fid_out,'%20s;%20s;%10.4f\n',Material{1},Material{2},kij_fit);
        end
        for im = 1:nM
            if (strcmpi(Material{1},Material_1{im}) && strcmpi(Material{2},Material_2{im})) ...
                || (strcmpi(Material{1},Material_2{im}) && strcmpi(Material{2},Material_1{im})) 
                Lfound = 1;
                fprintf(fid_out,'%20s;%20s;%10.4f\n',Material{1},Material{2},kij_fit);
            else
                fprintf(fid_out,'%20s;%20s;%10.4f\n',Material_1{im},Material_2{im},kij_all(im));
            end
            if im == nM && Lfound == 0
                fprintf(fid_out,'%20s;%20s;%10.4f\n',Material{1},Material{2},kij_fit);
            end
        end
        fclose(fid_out);
        delete 'Classes/Bin_kij_fit_2.txt'
%     catch
%         copyfile 'Classes/Bin_kij_fit_2.txt' Classes/Bin_kij_fit.txt 
%         delete 'Classes/Bin_kij_fit_2.txt'
%         error(['Error occurs'])
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%