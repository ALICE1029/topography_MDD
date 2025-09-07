load dida
voxel_num=45892;
network_num=18;
method='NGSR';
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load('/home/cxpang/matlab/code/8.brain_age/using_mat/index_all_remove65sex.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%age,group,mfd,sex

%% 
hc=find(group==1);
mdd=find(group==2);
sub1_index=mdd(find(Idx==1));
sub2_index=mdd(find(Idx==2));

group_sub1=zeros(length(hc)+length(sub1_index),1);
for i=1:length(hc)
    group_sub1(i)=1;
end
for i=length(hc)+1:length(hc)+length(sub1_index)
    group_sub1(i)=2;
end
group_sub2=zeros(length(hc)+length(sub2_index),1);
for i=1:length(hc)
    group_sub2(i)=1;
end
for i=length(hc)+1:length(hc)+length(sub2_index)
    group_sub2(i)=2;
end
sex_hc=sex(hc);
sex_mdd=sex(mdd);
sex_sub1=sex(sub1_index);
sex_sub2=sex(sub2_index);
sex_all_sub1=[sex_hc;sex_sub1];
sex_all_sub2=[sex_hc;sex_sub2];
age_hc=age(hc);
age_sub1=age(sub1_index);
age_sub2=age(sub2_index);
mfd_hc=mfd(hc);
mfd_sub1=mfd(sub1_index);
mfd_sub2=mfd(sub2_index);
site_hc=site(hc);
site_sub1=site(sub1_index);
site_sub2=site(sub2_index);
%% 2. age
all_age = [age_sub1; age_sub2];
g1 = repmat({'subtype1'},length(sub1_index),1);
g2 = repmat({'subtype2'},length(sub2_index),1);
[h,p]=ttest2(age_sub1,age_sub2)
g = [g1; g2];
boxplot(all_age,g)
%% sex
sex_1_1=length(find(sex_sub1==1))/length(sex_sub1)%female
sex_1_2=length(find(sex_sub1==2))/length(sex_sub1)%male
sex_2_1=length(find(sex_sub2==1))/length(sex_sub2)%female
sex_2_2=length(find(sex_sub2==2))/length(sex_sub2)%male
%spss
%% 3. hamd
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_hamd.txt');
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd.mat'))
%regress age sex 
des = [sex(mdd_hamd),age(mdd_hamd),site(mdd_hamd)];
[hamd, b, stats] = regress_out(hamd, des);

hamd_all=[mdd_hamd,hamd];
sub_hamd_1_index=intersect(mdd_hamd,sub1_index);
sub_hamd_2_index=intersect(mdd_hamd,sub2_index);
for i=1:length(sub_hamd_1_index)
    hamd_sub1(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_1_index(i)),2);
end
for i=1:length(sub_hamd_2_index)
    hamd_sub2(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_2_index(i)),2);
end
[t,p]=ttest2(hamd_sub1,hamd_sub2);

all_hamd = [hamd_sub1,hamd_sub2];
g1 = repmat({'subtype1'},length(sub_hamd_1_index),1);
g2 = repmat({'subtype2'},length(sub_hamd_2_index),1);
g = [g1; g2];
boxplot(all_hamd,g)

%% 4. subscale
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_sub.txt');
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd_sub.mat'))
%regress age sex 
des = [sex(mdd_hamd),age(mdd_hamd),site(mdd_hamd)];
for i=1:17
[hamd(:,i), b, stats] = regress_out(hamd(:,i), des);
end
for j=1:17
    hamd_all=[mdd_hamd,hamd(:,j)]
    sub_hamd_1_index=intersect(sub1_index,mdd_hamd);
    sub_hamd_2_index=intersect(mdd_hamd,sub2_index);
    hamd_sub_sub1=zeros(length(sub_hamd_1_index),1);
    hamd_sub_sub2=zeros(length(sub_hamd_2_index),1);
    for i=1:length(sub_hamd_1_index)
        hamd_sub_sub1(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_1_index(i)),2);
    end
    for i=1:length(sub_hamd_2_index)
        hamd_sub_sub2(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_2_index(i)),2);
    end
    [t(j),p(j)]=ttest2(hamd_sub_sub1,hamd_sub_sub2)
    all_hamd = [hamd_sub_sub1;hamd_sub_sub2];
%     g1 = repmat({'subtype1'},length(sub_hamd_1_index),1);
%     g2 = repmat({'subtype2'},length(sub_hamd_2_index),1);
%     g = [g1; g2];
%     boxplot(all_hamd,g)
    
end
q=mafdr(p,'BHFDR', true);
%% medicated
load('/HeLabData2/cxpang/DIDA/using_mat/clinical_pos_remove65.mat')

sub_hamd_1_med_0_index=intersect(med_0,sub1_index);
sub_hamd_1_med_1_index=intersect(med_1,sub1_index);
med_1_sub=length(sub_hamd_1_med_0_index)/(length(sub_hamd_1_med_0_index)+length(sub_hamd_1_med_1_index))

sub_hamd_2_med_0_index=intersect(med_0,sub2_index);
sub_hamd_2_med_1_index=intersect(med_1,sub2_index);
med_2sub=length(sub_hamd_2_med_0_index)/(length(sub_hamd_2_med_0_index)+length(sub_hamd_2_med_1_index))

%% first episode
sub_hamd_1_first_0_index=intersect(first_0,sub1_index);
sub_hamd_1_first_1_index=intersect(first_1,sub1_index);
first_1_sub=length(sub_hamd_1_first_0_index)/(length(sub_hamd_1_first_0_index)+length(sub_hamd_1_first_1_index))

sub_hamd_2_first_0_index=intersect(first_0,sub2_index);
sub_hamd_2_first_1_index=intersect(first_1,sub2_index);
first_2sub=length(sub_hamd_2_first_0_index)/(length(sub_hamd_2_first_0_index)+length(sub_hamd_2_first_1_index))
%% onset age
save(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/onset_pos_remove65.mat'),'onset_age','onset_age_pos')
on_all=[onset_age_pos,onset_age];
sub_on_1_index=intersect(onset_age_pos,sub1_index);
sub_on_2_index=intersect(onset_age_pos,sub2_index);
for i=1:length(sub_on_1_index)
    on_sub1(i)=on_all(find(on_all(:,1)==sub_on_1_index(i)),2);
end
for i=1:length(sub_on_2_index)
    on_sub2(i)=on_all(find(on_all(:,1)==sub_on_2_index(i)),2);
end
[t,p]=ttest2(on_sub1,on_sub2);

all_on = [on_sub1,on_sub2];
g1 = repmat({'subtype1'},length(sub_on_1_index),1);
g2 = repmat({'subtype2'},length(sub_on_2_index),1);
g = [g1; g2];
boxplot(all_on,g)

%% duration
load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/dur_pos_remove65.mat'),'duration_pos','duration')
%regress age sex 
%des = [sex(duration_pos),age(duration_pos)];

%
du_all=[duration_pos,duration];
sub_du_1_index=intersect(duration_pos,sub1_index);
sub_du_2_index=intersect(duration_pos,sub2_index);
for i=1:length(sub_du_1_index)
    du_sub1(i)=du_all(find(du_all(:,1)==sub_du_1_index(i)),2);
end
for i=1:length(sub_du_2_index)
    du_sub2(i)=du_all(find(du_all(:,1)==sub_du_2_index(i)),2);
end
[t,p]=ttest2(du_sub1,du_sub2)

all_du = [du_sub1,du_sub2];
g1 = repmat({'subtype1'},length(sub_du_1_index),1);
g2 = repmat({'subtype2'},length(sub_du_2_index),1);
g = [g1; g2];
boxplot(all_du,g)
