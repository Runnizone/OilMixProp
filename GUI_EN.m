clear;close;clc;  

Language  = 'DE';   % 'EN' or 'DE'
Executable = 1;     % 1 yes   or    0 no

if Executable
    BasePath =  'C:\Program Files (x86)\Chemnitz University of Technology\OilMixProp\application';
    if ~exist(BasePath,'dir')
        BasePath = 'C:\Program Files\Chemnitz University of Technology\OilMixProp\application';
    end
else
    BasePath = '.\';
    path(path,[pwd,'\Classes']);
    rmpath('C:\Xiaoxian\Github\OilMixProp\Classes');
    rmpath('C:\Xiaoxian\Github\CubicEoS\Classes');
    % rmpath('C:\Xiaoxian\Github\OilMixProp_MatLabExe\Classes');
end


format short; warning('off');
% inputs = {BasePath,Language};
oilmixprop_v1(BasePath,Language);
