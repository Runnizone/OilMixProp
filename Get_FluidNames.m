%%%%%%%%%%%%%%%%%%%%%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'/Classes']); format short;  AllEOS = {'PR','SRK','PTV','YFR'};      warning('off'); 
% rmpath('C:\Xiaoxian\Github\OilMixProp\Classes')

% different type of inputs 
% All: all fluids
% REFPROP: all fluids avaiable in REFPROP 10.0
% Others: all fluids except for those in REFPROP 10.0
% oils: all oils which are fitted using OilMixProp
% CASID: the fluid name with the given CASID

% INPUT = 'All';
% [Fluids, CASID] = fluidname(INPUT);


INPUT = '106-99-0';




