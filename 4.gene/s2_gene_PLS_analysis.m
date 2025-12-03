%% Gene Association Analysis
%% define work dir
clear all
data_dir = 'D:\study\sub2\liu\';
script_dir = 'D:\study\sub2\code\code\gene\';
sub='2'
sign='abs'
%% load parcellation,Z-maps,gene
hdr_par = spm_vol([script_dir,'rglasser360MNI.nii']);
vol_par =spm_read_vols(hdr_par);

V1=spm_vol([data_dir,'corr_sub2_',sub,'_remove65sex\t_combine_',sign,'.nii.gz']);
Y1=spm_read_vols(V1);

correlation = zeros(360,1);
overlap_rate = zeros(360,1);
for i=1:360
    correlation(i,1)=mean(Y1(vol_par==i));
  overlap_rate(i,1)=length(find(Y1(vol_par==i)))/length(Y1(vol_par==i));%有值的/该脑区voxel数目
end
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat'],'parcelExpression')
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat'],'probeInformation')

%% remove missing roi
temp1=find(overlap_rate<0.5);
temp2=find(isnan(parcelExpression(:,2)));
missingdata_regions=union(temp1,temp2);
region_ind=setdiff(parcelExpression(:,1),missingdata_regions);

group_express=parcelExpression(region_ind,2:end);
gene_name = probeInformation.GeneSymbol;

GENEdata=group_express;
MRIdata=correlation(region_ind);

%% PLS_calculation
Y = zscore(MRIdata);
dim =10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(GENEdata,Y,dim,'CV',dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);

%align PLS components with desired direction%
R1 = corr([XS(:,1),XS(:,2),XS(:,3)],MRIdata);
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end
%% calculate corrected weight gene
gene_name = probeInformation.GeneSymbol;
geneindex=1:size(GENEdata,2);
genes = gene_name;
X=GENEdata;
Y=zscore(MRIdata);
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
a=XS(:,1);
b=YS(:,1);
tmp=[a,b];
filename = strcat('D:\study\sub2\scatterplot\',sign,'_gene_sub',sub,'.txt');
dlmwrite(filename, tmp, 'delimiter', '\t');
[R1,p1]=corr(XS(:,1),YS(:,1));
[R1,p1]=corr(XS(:,1),MRIdata);
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
PLS1_ROIscores_280=XS(:,1);
save([data_dir,'\PLS1_ROIscore_sub',sub,sign,'.mat'],'PLS1_ROIscores_280');

PLS1weights = zeros(10027,1000);
%PLS2weights = zeros(10027,10000);
bootnum=1000;
parfor i=1:bootnum
    i
    myresample = randsample(size(X,1),length(region_ind));
    %res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVARb,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run
       
end

PLS1sw = std(PLS1weights');
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
% fid1 = fopen([data_dir,'\PLS1_geneWeights_sub',sub,sign,'.csv'],'w');
% for i=1:length(genes)
%   fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
% end
% fclose(fid1);
% new=PLS1(find(Z1>3|Z1<-3))
% fid1 = fopen([data_dir,'\PLS1_geneWeights_sub',sub,sign,'.csv'],'w');
% for i=1:length(new)
%   fprintf(fid1,'%s\n', new{i});
% end
% fclose(fid1);

% % permutation test using brainsmash to generate the surrogate maps for t value
load('D:\study\sub2\gene\sub2\surrogate_g1.mat')
hdr_mask = spm_vol('D:\study\sub2\liu\GMMask_3mm.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
for j = 1:1000
    disp(j);
    gradients_sur = zeros(360,1);
    Y1 = zeros(hdr_par.dim);    
    Y1(ind) = surrogate_maps(j,:);
    for i=1:360
        gradients_sur(i,1)=mean(Y1(vol_par==i));
    end
    MRIdata_s = gradients_sur(region_ind);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(GENEdata,MRIdata_s,dim);%Percentage of variance explained by the model, returned as a numeric matrix. PCTVAR is a 2-by-ncomp matrix, where ncomp is the number of PLS components. The first row of PCTVAR contains the percentage of variance explained in X by each PLS component, and the second row contains the percentage of variance explained in Y.
    PCTVARrand(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);    
end
p_single = zeros(1,10);
for l=1:dim
    p_single(l)=length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))/1000;
end
p_cum = length(find(Rsq>=Rsquared))/1000;
myStats=[PCTVAR; p_single];



% csvwrite([data_dir,'\PLS_stats.csv'],myStats);
% 
%% Draw variance explanation

py = plot(sort(myStats(2,:)','descend'),'.-','LineWidth',2);
py.Color = [115 130 184]/255;
py.MarkerSize = 20;
hold on
plot(0:11,0.1*ones(1,12),'--','Color',[226,115,134]./255,'LineWidth',1);
hold off

xlabel('PLS Component');
ylabel('Explained variance');
set(gca,'XLim',[0,11]);
set(gca,'YLim',[0,0.25],'YTick',0:0.1:0.25);
%t1 = text(1,0.20,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
%t2 = text(2,0.36,'**','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
set(gca,'LineWidth',0.5);

set(gca,'FontName','Arial','FontSize',14);
box off

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 9 5.6]);
print(gcf,[figure_dir,'PLS_Variance_main.tif'],'-dtiff','-r1000')

%% calculate corrected weight gene
gene_name = probeInformation.GeneSymbol;
geneindex=1:size(GENEdata,2);
genes = gene_name;
X=GENEdata;
Y=zscore(MRIdata);
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
a=XS(:,1);
b=YS(:,1);
tmp=[a,b];
filename = strcat('D:\study\sub2\scatterplot\',sign,'_gene_sub',sub,'.txt');
dlmwrite(filename, tmp, 'delimiter', '\t');
[R1,p1]=corr(XS(:,1),YS(:,1));
[R1,p1]=corr(XS(:,1),MRIdata);
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
PLS1_ROIscores_280=XS(:,1);
save([data_dir,'\PLS1_ROIscore_sub',sub,sign,'.mat'],'PLS1_ROIscores_280');
%csvwrite([data_dir,'\PLS1_ROIscores_sub',sub,sign,'.csv'],XS(:,1));
PLS1_score=XS(:,1);


PLS1weights = zeros(10027,1000);

bootnum=1000;
parfor i=1:bootnum
    i
    myresample = randsample(size(X,1),length(region_ind));
    %res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run
end

PLS1sw = std(PLS1weights');
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
fid1 = fopen([data_dir,'\PLS1_geneWeights_sub',sub,sign,'.csv'],'w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
end
fclose(fid1);


%% generate 360weighted image for BrainNet Viewer
sub='1';
sign='abs'
hdr_par = spm_vol([script_dir,'rglasser360MNI.nii']);
vol_par =spm_read_vols(hdr_par);

V1=spm_vol([data_dir,'corr_sub2_',sub,'_remove65sexmask\t_combine_',sign,'.nii.gz']);
Y1=spm_read_vols(V1);

correlation = zeros(360,1);
overlap_rate = zeros(360,1);
for i=1:360
    correlation(i,1)=mean(Y1(vol_par==i));
  overlap_rate(i,1)=length(find(Y1(vol_par==i)))/length(Y1(vol_par==i));%有值的/该脑区voxel数目
end
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat'],'parcelExpression')
load([script_dir,'100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat'],'probeInformation')
temp1=find(overlap_rate<0.5);
temp2=find(isnan(parcelExpression(:,2)));
missingdata_regions=union(temp1,temp2);
region_ind=setdiff(parcelExpression(:,1),missingdata_regions);
load([data_dir,'\PLS1_ROIscore_sub',sub,sign,'.mat'],'PLS1_ROIscores_280');
gii1 = gifti([script_dir,'Glasser180_210P_L.label.gii']);
gii2 = gifti([script_dir,'Glasser180_210P_R.label.gii']);
ParcelLabel = double([gii1.cdata;gii2.cdata]);

PLS1_weight = zeros(length(ParcelLabel),1);
for i = 1:length(region_ind)
    PLS1_weight(ParcelLabel==region_ind(i)) = PLS1_ROIscores_280(i);
end
save([data_dir,'\PLS1_weight',sub,sign,'.txt'],'PLS1_weight','-ascii');
%% correlation between PLS score and G1 difference
load PLS1_ROIscore.mat
load PLS2_ROIscore.mat
load region_ind.mat
V=spm_vol('4mm_Glasser360.nii');
Y = spm_read_vols(V);

V1 = spm_vol('g1_Z_T2.nii');
Y1 = spm_read_vols(V1);

gradients = zeros(360,1);

for i=1:360
    gradients(i,1)=mean(Y1(Y==i));
end
corr_real = zeros(2,1);
corr_real(1) = corr(gradients(region_ind),PLS1_ROIscores_280);
corr_real(2) = corr(gradients(region_ind),PLS2_ROIscores_280);


load surrogate_g1.mat
hdr_mask = spm_vol('Reslice_group_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
corr_surr = zeros(2,10000);
gradients_surr = zeros(360,1);
for j = 1:10000
    disp(j);
    Y_s = zeros(V.dim);
    Y_s(ind) = surrogate_g1(j,:);
    for i = 1:360
        gradients_surr(i) = mean(Y_s(Y==i));
    end
    corr_surr(1,j) = corr(gradients_surr(region_ind),PLS1_ROIscores_280);
    corr_surr(2,j) = corr(gradients_surr(region_ind),PLS2_ROIscores_280);
end
p = zeros(2,1);
p(1) = length(find(corr_surr(1,:)>corr_real(1)))/10000;
p(2) = length(find(corr_surr(2,:)>corr_real(1)))/10000;
save corr_surr.mat corr_surr corr_real p

