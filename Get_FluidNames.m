%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  AllEOS = {'PR','SRK','PTV','YFR'};      warning('off'); 
% rmpath('C:\Xiaoxian\Github\OilMixProp\Classes')

% different type of inputs 
% All: all fluids
% REFPROP: all fluids avaiable in REFPROP 10.0
% Others: all fluids except for those in REFPROP 10.0
% oils: all oils which are fitted using OilMixProp
% CASID: the fluid name with the given CASID


% The adopted name of a fluid in coding package
% CAS_RN = '106-99-0';
% f_fluidname(CAS_RN); 

% 
% % fluids avaiable in REFPROP
% f_fluidname('REFPROP') 


% the current existing oils 
f_fluidname('oils') 


% % complete list of fluids
% f_fluidname('all') 