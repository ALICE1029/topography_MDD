addpath(genpath('/home/cxpang/matlab/myPLS-master'));
addpath(genpath('/home/cxpang/matlab/Collaborative_Brain_Decomposition-master/lib/NIfTI_20140122'))
addpath(('/home/cxpang/matlab/code'));
method='NGSR'
addpath('/home/cxpang/matlab/code/8.brain_age/function')
%% loading
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub2_remove65.mat','loading_sub2')
load('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hamd_sub_sub2_remove65.mat','age','sex','mfd')
des = [age,sex,mfd];
for i=1:287276
[loading_sub2(:,i), b, stats] = regress_out(loading_sub2(:,i), des);
end
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub2_regress_remove65.mat','loading_sub2')
%% 
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub1_remove65.mat','loading_sub1')
load('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hamd_sub_sub1_remove65.mat','age','sex','mfd')
des = [age,sex,mfd];
for i=1:287276
[loading_sub1(:,i), b, stats] = regress_out(loading_sub1(:,i), des);
end
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub1_regress_remove65.mat','loading_sub1')
%% ---------- Input data ----------
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub2_regress_remove65.mat','loading_sub2')
X0 = loading_sub2;
load('/home/cxpang/matlab/code/8.brain_age/using_mat/hamd_sub_sub2_remove65.mat')
Y0behav = hamd; % sex - 0=female/1=male
%regress hamd
load('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hamd_sub_sub2_remove65.mat','age','sex','site')
des = [age,sex];
for i=1:17
[Y0behav(:,i), b, stats] = regress_out(Y0behav(:,i), des);
end
hamd=Y0behav;


%%
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub1_regress_remove65.mat','loading_sub1')
X0 = loading_sub1;
load('/home/cxpang/matlab/code/8.brain_age/using_mat/hamd_sub_sub1_remove65.mat')
Y0behav = hamd; % sex - 0=female/1=male
%regress hamd
load('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hamd_sub_sub1_remove65.mat','age','sex')
des = [age,sex];
for i=1:17
[Y0behav(:,i), b, stats] = regress_out(Y0behav(:,i), des);
end
%%
% consistent with 10f cv
[my_corr,p]=corr(X0,Y0behav,'type','Pearson');
pos = find(any(p<0.05, 2));
%save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/sub1/p0.05/all_mask.mat','pos')
for i=1:287276
    tmp(i)=i;
end
pos=setdiff(tmp,pos);
X0(:,pos)=[];

input.grouping = zeros(size(hamd,1),1);
pls_opts.nPerms = 100;
pls_opts.nBootstraps =1;
% if you want to regress nuisance variables, do it here

% --- brain data ---
% Matrix X0 is typically a matrix with brain imaging data,
% of size subjects (rows) x imaging features (columns)
input.brain_data=X0;

% --- behavior data ---
input.behav_data = Y0behav;

% --- grouping data ---
input.group_names={'group 1'}; 

% --- Names of the behavior data ---
input.behav_names = {'hamd'}; 

% --- Names of the imaging data ---

for ii = 1:20; input.img_names{ii,1} = ['img ' num2str(ii)]; end


%% ---------- Options for PLS ----------

% --- Permutations & Bootstrapping ---


% --- Data normalization options ---
pls_opts.normalization_img =1;
pls_opts.normalization_behav = 1;

% --- PLS grouping option ---
pls_opts.grouped_PLS = 0; 

% --- Permutations grouping option ---

pls_opts.grouped_perm = 0;

% --- Bootstrapping grouping option ---

pls_opts.grouped_boot = 0;

% --- Mode for bootstrapping procrustes transform ---

pls_opts.boot_procrustes_mod = 2;

% --- Save bootstrap resampling data? ---

pls_opts.save_boot_resampling=1;

% --- Type of behavioral analysis ---

pls_opts.behav_type = 'behavior';

%% ---------- Options for result saving and plotting ----------
% --- path where to save the results ---
save_opts.output_path = './example_results';

% --- prefix of all results files ---
% this is also the default prefix of the toolbox if you don't define
% anything
save_opts.prefix = sprintf('myPLS_%s_norm%d-%d',pls_opts.behav_type,...
    pls_opts.normalization_img, pls_opts.normalization_behav);

% --- Plotting grouping option ---
% 0: Plots ignoring grouping
% 1: Plots considering grouping
save_opts.grouped_plots = 1;

% --- Significance level for latent components ---
save_opts.alpha = 0.1; % for the sake of the example data



% uncomment the following to see example for barplot figures:
s.img_type = 'barPlot';
save_opts.fig_pos_img = [440   606   560   192];

% --- Brain visualization thresholds ---
% (thresholds for bootstrap scores and loadings, only required if imagingType='volume' or imagingType='corrMat')
save_opts.BSR_thres = [-2.3 2.3]; % negative and positive threshold for visualization of bootstrap ratios
save_opts.load_thres = [-0.4 0.4]; % negative and positive threshold for visualization of loadings

% --- Brain mask ---
% (gray matter mask, only required if imagingType='volume')
save_opts.mask_file = 'example_mask.nii'; % filename of binary mask that will constrain analysis

% --- Structural template file for visualization ---
% (structural volume for background, only required if imagingType='volume')
save_opts.struct_file = 'example_struct.nii';

% --- Orientation of volumes in slice plots ---
% (only required if imagingType='volume')
save_opts.volume_orientation = 'axial'; %'axial','coronal','sagittal'

% --- Bar plot options ---
save_opts.plot_boot_samples = 1; % binary variable indicating if bootstrap samples should be plotted in bar plots
save_opts.errorbar_mode = 'CI'; % 'std' = plotting standard deviations; 'CI' = plotting 95% confidence intervals
save_opts.hl_stable = 1; % binary variable indicating if stable bootstrap scores should be highlighted

% --- Customized figure size for behavior bar plots ---
save_opts.fig_pos_behav = [440   606   320   192];
%%
res = myPLS_analysis(input,pls_opts);
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/bootstrap/sub1_noboot.mat','res','save_opts')
%% boot
res = myPLS_analysis(input,pls_opts);
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/bootstrap/sub1.mat','res','save_opts')
myPLS_plot_results(res,save_opts);
%%
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/bootstrap/sub1.mat','res','save_opts')
boot=res.boot_results.LC_behav_loadings_boot;
boot=boot(:,1,:);
boot_tmp=squeeze(boot);

l=res.boot_results.LC_img_loadings_lB(:,1);
r=res.boot_results.LC_img_loadings_uB(:,1);
weight=res.LC_img_loadings(:,1);
k=1
tmp=[]
for i=1:length(l)
    if(l(i)<0&r(i)>0)
        tmp(k)=i;
        k=k+1;
    end
end
weight(tmp)=0;

l=res.boot_results.LC_behav_loadings_lB(:,1);
r=res.boot_results.LC_behav_loadings_uB(:,1);
k=1
tmp=[]
for i=1:17
    if(l(i)<0&r(i)>0)
        tmp(k)=i;
        k=k+1;
    end
end
%%
method='NGSR'
voxel_num=45892;
network_num=18;
xweight=zeros(287276,1);
load('/home/cxpang/matlab/code/8.brain_age/sub2_pos_0.05.mat','pos')
xweight(pos)=weight;
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/hamd/loading_pos_all.mat'),'pos')
x_tmp=zeros(voxel_num*network_num,1);
tmp=[];
for i=1:voxel_num*network_num
    tmp(i)=i;
end
tmp2=setdiff(tmp,pos);
x_tmp(tmp2)=xweight;
addpath(genpath('/home/cxpang/matlab/Collaborative_Brain_Decomposition-master/lib/NIfTI_20140122'))
analyse=reshape(x_tmp,[voxel_num,network_num]);
%%
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamdsub_weight/analyse_sub1_boot_remove65.mat','analyse')
%one map
one=abs(analyse);
all=mean(one,2);
hdr_mask = spm_vol('/HeLabData2/cxpang/fmri/mask/GMMask_3mm.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
vol = zeros(hdr_mask.dim);
vol(ind) = all;% z value nii
[L, num] = bwlabeln(vol,26);
n = 0;
cluster_size=50;
for x = 1:num
    theCurrentCluster = L == x;
   len(x)= length(find(theCurrentCluster));
    if length(find(theCurrentCluster)) <= cluster_size
        n = n + 1;
       vol (logical(theCurrentCluster)) = 0;
    
    end
end
outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/hamdsub1_weight_remove65');
saveFig = 0;
new_vol=vol(ind);

maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
maskNii = load_untouch_nii(maskName);
if ~exist(outDir,'dir')
    mkdir(outDir);
end
kNii = maskNii;
kNii.img(maskNii.img~=0) =new_vol;
outName = ['/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamdsub1_weight_remove65/allweight_bo_thre50.nii.gz'];
save_untouch_nii(kNii,outName);

%num
weight_neg=analyse;
weight_neg(find(analyse>0))=0;
weight_neg_sum=zeros(18,1);
weight_pos_sum=zeros(18,1);
for i=1:18
    weight_neg_sum(i)=sum (weight_neg(:,i));
end
weight_pos=analyse;
weight_pos(find(analyse<0))=0;
for i=1:18
    weight_pos_sum(i)=sum (weight_pos(:,i));
end
weight_sum=[weight_pos_sum,abs(weight_neg_sum)];

save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamdsub1_weight_remove65/sum.mat','weight_sum')
