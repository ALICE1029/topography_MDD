clear all
close all
inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\vali\');
outpath=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\vali_image\');
%% validation:loso
for s=1:10
    s
yfitF=nan(1000,5);
yfitM=nan(1000,5);
load(strcat(inpath2,...
    'GAMLSS_','ot1npsighc_',num2str(s),'_predicted_sex0_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load(strcat(inpath2,...
    'GAMLSS_','ot1npsighc_',num2str(s),'_predicted_sex1_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile

fit=yfit(:,3);
hold on
% gap>0
load(strcat(inpath2,...
    'GAMLSS_','ot1npsigsub1_',num2str(s),'_predicted_sex0_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load(strcat(inpath2,...
    'GAMLSS_','ot1npsigsub1_',num2str(s),'_predicted_sex1_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
fit_pos=yfit_pos(:,3);
loso_sub1{s}=yfit_pos(:,3);% for compr
plot(fake_age_pos,(yfit_pos(:,3)),'b','linewidth',2.5); %plot 50th centile]
set(gca, 'FontSize', 16);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub1_abs', num2str(s), '.png'); % 定义输出文件名  
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式      
    
yfitF=nan(1000,5);
yfitM=nan(1000,5);
load(strcat(inpath2,...
    'GAMLSS_','ot2npsighc_',num2str(s),'_predicted_sex0_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load(strcat(inpath2,...
    'GAMLSS_','ot2npsighc_',num2str(s),'_predicted_sex1_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);
hold on
% gap>0
load(strcat(inpath2,...
    'GAMLSS_','ot2npsigsub2_',num2str(s),'_predicted_sex0_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load(strcat(inpath2,...
    'GAMLSS_','ot2npsigsub2_',num2str(s),'_predicted_sex1_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
fit_pos=yfit_pos(:,3);
loso_sub2{s}=yfit_pos(:,3);
plot(fake_age_pos,(yfit_pos(:,3)),'r','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 16);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub2_abs', num2str(s), '.png'); % 定义输出文件名  
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式      
    
end
%% mfd
yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','omfdt1npsighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','omfdt1npsighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);

hold on
% gap>0
load([inpath2,...
    'GAMLSS_','omfdt1npsigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','omfdt1npsigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
mfd_sub1=yfit_pos(:,3);
fake_age_pos=(faf_age+fam_age)/2;
plot(fake_age_pos,(yfit_pos(:,3)),'b','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 14);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub1_abs_mfd.png'); % 定义输出文件名  
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式  


yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','omfdt2npsighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','omfdt2npsighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);

hold on
% gap>0
load([inpath2,...
    'GAMLSS_','omfdt2npsigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','omfdt2npsigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
fit_pos=yfit_pos(:,3);
mfd_sub2=yfit_pos(:,3);
plot(fake_age_pos,(yfit_pos(:,3)),'r','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 14);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub2_abs_mfd.png'); % 定义输出文件名  
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式  
%% sex
for s=1:2
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','osext1npsighc_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
    hold on
    plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
    hold on
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','osext1npsigsub1_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hold on
    sex_sub1{s}=yfit(:,3);
    plot(fake_age,(yfit(:,3)),'b','linewidth',2.5); %plot 50th centile
    set(gca, 'FontSize', 14);
    outputFileName = strcat(outpath, 'GAMLSS_plot_sub1_abs_sex',num2str(s),'.png'); % 定义输出文件名
    saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式
 
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','osext2npsighc_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
     hold o  
    plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
    hold on
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','osext2npsigsub2_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hold on
      sex_sub2{s}=yfit(:,3);
    plot(fake_age,(yfit(:,3)),'r','linewidth',2.5); %plot 50th centile
    set(gca, 'FontSize', 14);
    outputFileName = strcat(outpath, 'GAMLSS_plot_sub2_abs_sex',num2str(s),'.png'); % 定义输出文件名
    saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式
   
end

%% srpbs
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\vali\');
yfitF=nan(1000,5);
yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t1asighcsr_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1asighcsr_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);

hold on
% gap>0
load([inpath2,...
    'GAMLSS_','t1asigsub1sr_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1asigsub1sr_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
plot(fake_age_pos,(yfit_pos(:,3)),'b','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 14);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub1_abs_srpbs.png'); % 定义输出文件名  
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式  


yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t2asighcsr_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asighcsr_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);

hold on
% gap>0
load([inpath2,...
    'GAMLSS_','t2asigsub2sr_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asigsub2sr_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
fit_pos=yfit_pos(:,3);
plot(fake_age_pos,(yfit_pos(:,3)),'r','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 14);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub2_abs_srpbs.png'); % 定义输出文件名  
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式  
%% bootstrapyfitF=nan(1000,5);
clear all
close all

inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\boot_abs\');
for i=1:100
    i
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_boot(:,i)=(faf_age+fam_age)/2;
    fit_boot(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_sub1_boot(:,i)=(faf_age+fam_age)/2;
    fit_sub1_boot(:,i)=yfit(:,3);
end
% 95%ci
for i=1:1000
    mean_normal = mean(fit_boot(i,:));
    std_normal = std(fit_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_hc(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
    
    
    mean_normal = mean(fit_sub1_boot(i,:));
    std_normal = std(fit_sub1_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_sub1_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_sub1(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t1psighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

figure; % Ensure it's a new figure
set(gcf, 'Position', [100, 100, 800, 400]);
%gap<0
load([inpath2,...
    'GAMLSS_','t1npsigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1npsigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
%
x = [fake_age', fliplr(fake_age')];
y1=ci_normal_hc(:,1);
y2=ci_normal_hc(:,2);
y = [(y1)', fliplr(y2')];
fill(x, y, [0.6 0.6 0.6],'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile

hold on
% plot(fake_age,ci_normal_sub1(:,1),'b--','linewidth',1.5); %plot 50th centile
% hold on
% plot(fake_age,ci_normal_sub1(:,2),'b--','linewidth',1.5); %plot 50th centile
x = [fake_age', fliplr(fake_age')];
y = [ci_normal_sub1(:,1)', fliplr(ci_normal_sub1(:,2)')];
fill(x, y, [0.7, 0.7, 1], 'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age_neg,(yfit_neg(:,3)),'b','linewidth',2.5); %plot 50th centile
% *
% 定义第一个横坐标范围和纵坐标
x_stars1 = 11:21;
y_stars1 = 6;

% 绘制灰色横线
plot(x_stars1, ones(size(x_stars1)) * y_stars1, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);

% 在横线中间绘制星号
hold on;
plot(mean(x_stars1), y_stars1 + 0.2, '*', 'MarkerSize', 8, 'Color', [0.5, 0.5, 0.5]);  % 绘制星号

% 绘制竖直虚线到 x 轴
% plot([min(x_stars1), min(x_stars1)], [0, y_stars1], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);  % 左边竖直虚线
% plot([max(x_stars1), max(x_stars1)], [0, y_stars1], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);  % 右边竖直虚线

% 定义第二个横坐标范围和纵坐标
x_stars2 = 54:65;
y_stars2 = 7;

% 绘制第二条灰色横线
plot(x_stars2, ones(size(x_stars2)) * y_stars2, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);

% 在横线中间绘制星号
plot(mean(x_stars2), y_stars2 + 0.2, '*', 'MarkerSize', 8, 'Color', [0.5, 0.5, 0.5]);  % 绘制星号

% % 绘制第二条竖直虚线到 x 轴
% plot([min(x_stars2), min(x_stars2)], [0, y_stars2], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);  % 左边竖直虚线
% plot([max(x_stars2), max(x_stars2)], [0, y_stars2], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);  % 右边竖直虚线

% 去除边框和设置字体大小
box off;
set(gca, 'FontSize', 14);



box off
set(gca, 'FontSize', 14);

%% sub 2 pos
clear all
close all

inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\boot_abs\');
for i=1:100%i boot times

yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath3,...
    'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath3,...
    'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age_boot(:,i)=(faf_age+fam_age)/2;
fit_boot(:,i)=yfit(:,3);
% sub 1
yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath3,...
    'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath3,...
    'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age_sub1_boot(:,i)=(faf_age+fam_age)/2;
fit_sub1_boot(:,i)=yfit(:,3);
end
% 95%ci
for i=1:1000
mean_normal = mean(fit_boot(i,:));
std_normal = std(fit_boot(i,:));

confidence_level = 0.95; % 95%

df = (length(fit_boot(i,:))-1);
critical_value = tinv((1 + confidence_level) / 2, df);

ci_normal_hc(i,:) = [mean_normal+ - critical_value * std_normal , ...
 mean_normal + critical_value * std_normal];


mean_normal = mean(fit_sub1_boot(i,:));
std_normal = std(fit_sub1_boot(i,:));

confidence_level = 0.95; % 95%

df = (length(fit_sub1_boot(i,:))-1);
critical_value = tinv((1 + confidence_level) / 2, df);

ci_normal_sub1(i,:) = [mean_normal+ - critical_value * std_normal , ...
 mean_normal + critical_value * std_normal];
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

% plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
% hold on
%gap<0
load([inpath2,...
    'GAMLSS_','t2npsigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2npsigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
% 
figure; % Ensure it's a new figure
set(gcf, 'Position', [100, 100, 800, 400]);
x = [fake_age', fliplr(fake_age')];
y = [ci_normal_sub1(:,1)', fliplr(ci_normal_sub1(:,2)')];
fill(x, y, [1, 0.7, 0.7], 'EdgeColor', 'none','facealpha',0.5); % 灰色阴影
hold on
x = [fake_age', fliplr(fake_age')];
y1=ci_normal_hc(:,1);
y2=ci_normal_hc(:,2);
y = [(y1)', fliplr(y2')];%y的转置在里面，之前犯了先flip再转置的错误
fill(x, y, [0.5 0.5 0.5],'EdgeColor', 'none','facealpha',0.5); 
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
plot(fake_age_neg,(yfit_neg(:,3)),'r','linewidth',2.5); %plot 50th centile
hold on
% plot(fake_age,ci_normal_sub1(:,1),'b--','linewidth',1.5); %plot 50th centile
% hold on
% plot(fake_age,ci_normal_sub1(:,2),'b--','linewidth',1.5); %plot 50th centile

box off
set(gca, 'FontSize', 14);
x_stars1 = 13:50;
y_stars1 = 14;

% 绘制灰色横线
plot(x_stars1, ones(size(x_stars1)) * y_stars1, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);

% 在横线中间绘制星号
hold on;
plot(mean(x_stars1), y_stars1 + 0.3, '*', 'MarkerSize', 8, 'Color', [0.5, 0.5, 0.5]);  % 绘制星号
%% permut
% sub 1 pos
clear all
close all

inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\permut_abs_new\');
for i=1:1000
    i
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_permut(:,i)=(faf_age+fam_age)/2;
    fit_permut(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fit_sub1_permut(:,i)=yfit(:,3);
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t1asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
real_hc=yfit(:,3);

%gap<0
load([inpath2,...
    'GAMLSS_','t1psigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1psigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
real_sub1_pos=yfit_neg(:,3);
real_diff_pos_sub2=(real_sub1_pos-real_hc);
permut_diff_pos=[];
for i=1:1000
    permut_diff_pos(i,:)=(fit_sub1_permut(:,i)-fit_permut(:,i));
end
% for i=1:1000
% sig_sub_pos(i)=length(find( real_diff_pos_sub1(i)< permut_diff_pos(i,:)))/1000;
% end
for i=11:64
    pos=find(faf_age>=i&faf_age<i+1);
    real_sub1_pos1(i-10)=mean(real_diff_pos_sub1(pos));
    permut_diff_pos1=mean(permut_diff_pos(:,pos),2);

    sig_sub1_pos(i-10)=length(find( real_sub1_pos1(i-10)>= permut_diff_pos1))/1000;
end
 find(sig_sub1_pos<0.05)
q=mafdr(sig_sub1_pos,'BHFDR', true);

%%
clear all
close all

inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\permut_abs_new\');
for i=1:1000
    i
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_permut(:,i)=(faf_age+fam_age)/2;
    fit_permut(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fit_sub1_permut(:,i)=yfit(:,3);
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
real_hc=yfit(:,3);

%gap<0
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
real_sub1_pos=yfit_neg(:,3);
real_diff_pos_sub1=(real_sub1_pos-real_hc);
permut_diff_pos=[];
for i=1:1000
    permut_diff_pos(i,:)=(fit_sub1_permut(:,i)-fit_permut(:,i));
end
% for i=1:1000
% sig_sub_pos(i)=length(find( real_diff_pos_sub1(i)< permut_diff_pos(i,:)))/1000;
% end
for i=11:64
    pos=find(faf_age>=i&faf_age<i+1);
    real_sub1_pos1(i-10)=mean(real_diff_pos_sub1(pos));
    permut_diff_pos1=mean(permut_diff_pos(:,pos),2);

    sig_sub1_pos(i-10)=length(find( real_sub1_pos1(i-10)>=permut_diff_pos1))/1000;
end
 find(sig_sub1_pos<0.05)
q=mafdr(sig_sub1_pos,'BHFDR', true);
%%
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated males at site i (assumes males = 1)
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
plot(fake_age_neg,(yfit_neg(:,3)),'k','linewidth',2.5); %plot 50th centile
hold on
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated males at site i (assumes males = 1)
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
%
plot(fake_age_neg,(yfit_neg(:,3)),'r','linewidth',2.5); %plot 50th centile
load([inpath2,...
    'GAMLSS_','t2asigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated males at site i (assumes males = 1)
load([inpath2,...
    'GAMLSS_','t2asigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
%
plot(fake_age_neg,(yfit_neg(:,3)),'b','linewidth',2.5); %plot 50th centile
%%
for i=1:10
    [lr1(i),lp1(i)]=corr(loso_sub1{i},real_sub1);
    [lr2(i),lp2(i)]=corr(loso_sub2{i},real_sub2);
end
[mr1,mp1]=corr(mfd_sub1,real_sub1);
[mr2,mp2]=corr(mfd_sub2,real_sub2);
for i=1:2
    [sr1(i),sp1(i)]=corr(sex_sub1{i},real_sub1);
    [sr2(i),sp2(i)]=corr(sex_sub2{i},real_sub2);
end
%%
%% validation:loso
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\indi\');
for s=1:3514
    s
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','a2sighc_',num2str(s),'__predicted_sex0_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load(strcat(inpath2,...
        'GAMLSS_','a2sighc_',num2str(s),'__predicted_sex1_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age=(faf_age+fam_age)/2;
    

   
    fit(s,:)=yfit(:,3);

    % gap>0
    load(strcat(inpath2,...
        'GAMLSS_','a2sigsub2_',num2str(s),'__predicted_sex0_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load(strcat(inpath2,...
        'GAMLSS_','a2sigsub2_',num2str(s),'__predicted_sex1_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit_pos=(yfitF+yfitM)/2;
    fake_age_pos=(faf_age+fam_age)/2;
    fit_pos(s,:)=yfit_pos(:,3);
end
yfit_new=pca(fit');
[coeff, score, latent, ~, explained, ~] = pca(fit_pos');
explained_variance = latent / sum(latent);

yfit_pos_new=pca(fit_pos');
yfit_pca = fit'* yfit_new(:,:);
yfit_pca=yfit_pca(:,1);
yfit_pos_pca = fit_pos'* yfit_pos_new(:,:);
yfit_pos_pca=yfit_pos_pca(:,1);
  plot(fake_age,(yfit_pca),'k','linewidth',2.5); %plot 50th centile
  hold on
  plot(fake_age,(yfit_pos_pca),'r','linewidth',2.5); %plot 50th centile