load CSU
%%
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
method='NGSR'
voxel_num=45892;
network_num=18;
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
addpath(genpath('/home/cxpang/matlab/SVR'));
%%
load('/HeLabData2/cxpang/treatment_NGSR/CSU/result/regress_variable.mat','x')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/deltahamd.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/covariate.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/mfd.mat','mfd')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/mfd_co.mat','mfd_base','mfd_after')
%% topography changes
load('/home/cxpang/matlab/treatment/CSU_baseline.mat')
lo1=x(:,pos)';
load('/home/cxpang/matlab/treatment/CSU_after.mat','pos')
lo2=x(:,pos)';
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/group.mat','group')

for i=1:826056
    lo1(:,i)=regress_out(lo1(:,i), mfd_base);
    lo2(:,i)=regress_out(lo2(:,i), mfd_after);
end
%% generate change_sub1/2 folders, 18folder with sub_num subfolders, for
% paired ttest with GRF corrected 
lo1_sub1=lo1(find(group==1),:);
lo2_sub1=lo2(find(group==1),:);
lo1_sub2=lo1(find(group==2),:);
lo2_sub2=lo2(find(group==2),:);
lo1_sub1_cell=cell(18,1);
lo2_sub1_cell=cell(18,1);
lo1_sub2_cell=cell(25,1);
lo2_sub2_cell=cell(25,1);
for i=1:18
    tmp=reshape(lo1_sub1(i,:),[voxel_num,network_num]);
    lo1_sub1_cell{i}=tmp;
    tmp=reshape(lo2_sub1(i,:),[voxel_num,network_num]);
    lo2_sub1_cell{i}=tmp;
end
for i=1:25
    tmp=reshape(lo1_sub2(i,:),[voxel_num,network_num]);
    lo1_sub2_cell{i}=tmp;
    tmp=reshape(lo2_sub2(i,:),[voxel_num,network_num]);
    lo2_sub2_cell{i}=tmp;
end
for i=1:18
    i
    tmp=lo1_sub1_cell{i};
    tmp2=lo2_sub1_cell{i};
    for j=1:18%  network num
        
        outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/change_sub1_base/',num2str(j));
        saveFig = 0;
        maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
        maskNii = load_untouch_nii(maskName);
        if ~exist(outDir,'dir')
            mkdir(outDir);
        end
        kNii = maskNii;
        kNii.img(maskNii.img~=0) =tmp(:,j);
        outName = [outDir,filesep,strcat('t',num2str(i),'.nii.gz')];
        save_untouch_nii(kNii,outName);
        outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/change_sub1_after/',num2str(j));
        saveFig = 0;
        maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
        maskNii = load_untouch_nii(maskName);
        if ~exist(outDir,'dir')
            mkdir(outDir);
        end
        kNii = maskNii;
        kNii.img(maskNii.img~=0) = tmp2(:,j);
        outName = [outDir,filesep,strcat('t',num2str(i),'.nii.gz')];
        save_untouch_nii(kNii,outName);
    end
    
    
end

for i=1:25 
    i
    tmp=lo1_sub2_cell{i};
    tmp2=lo2_sub2_cell{i};
    for j=1:18%  network num
        
        outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/change_sub2_base/',num2str(j));
        saveFig = 0;
        maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
        maskNii = load_untouch_nii(maskName);
        if ~exist(outDir,'dir')
            mkdir(outDir);
        end
        kNii = maskNii;
        kNii.img(maskNii.img~=0) =tmp(:,j);
        outName = [outDir,filesep,strcat('t',num2str(i),'.nii.gz')];
        save_untouch_nii(kNii,outName);
        outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/change_sub2_after/',num2str(j));
        saveFig = 0;
        maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
        maskNii = load_untouch_nii(maskName);
        if ~exist(outDir,'dir')
            mkdir(outDir);
        end
        kNii = maskNii;
        kNii.img(maskNii.img~=0) = tmp2(:,j);
        outName = [outDir,filesep,strcat('t',num2str(i),'.nii.gz')];
        save_untouch_nii(kNii,outName);
    end
    
    
end
