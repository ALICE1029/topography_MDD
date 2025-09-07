function prepare_vars_for_norm_models(t,type,network,y, sex,age, covars)
% Inputs:

% y = IDP
% age = age
%/ sex = biological attribute of sex
% /site = scanner/site
% neuroimaging covariates (matrix)

% Outputs:
% variables prepared for normative modeling 

covars_norm=zeros(size(covars));
for i=1:size(covars,2)
    covars_norm(:,i)=inormal(covars(:,i));
end
%save variables for cross-sectional normative models
dir=strcat('/home/cxpang/matlab/code/8.brain_age/variables_for_normative_modeling/',t,'/');
mkdir(dir);
outstr=strcat(type,'_',num2str(network),'_',['lo_vars_for_R_cross']);
save([dir,outstr], 'y','sex','age','covars','-v6')     
end