load('D:\study\sub2\brain_age\brain_span\age_cortical_am_hip.mat','age_new_day')
load('D:\study\sub2\brain_age\brain_span\column_pos_cortical_am_hip.mat','column_pos')
data = readmatrix('D:\study\sub2\brain_age\brain_span\expression_matrix_my.csv');
data_new=data(:,column_pos);
filename = 'D:\study\sub2\brain_age\brain_span\columns_metadata.csv';
opts = detectImportOptions(filename);
opts.SelectedVariableNames = opts.VariableNames(7); 
strData = readmatrix(filename, opts, 'OutputType', 'string');
brain_new=strData(column_pos);
rowsToDelete = all(data_new == 0, 2);
data_new(rowsToDelete, :) = [];
filename = 'D:\study\sub2\brain_age\brain_span\rows_metadata.csv';
opts = detectImportOptions(filename);
opts.SelectedVariableNames = opts.VariableNames(4);  % 选择第一列
strData = readmatrix(filename, opts, 'OutputType', 'string');
gene_new=strData;
gene_new(rowsToDelete)=[];
%%

condi='sub1'
domi_pos=find(brain_new=='OFC'|brain_new=='ITC'|brain_new=='STC');
domi_notpos=find(~(brain_new=='OFC'|brain_new=='ITC'|brain_new=='STC'));

%condi='sub2'
%domi_pos=find(brain_new=='M1C'|brain_new=='V1C'|brain_new=='VFC'|brain_new=='HIP'|brain_new=='AMY'|brain_new=='DFC');
%domi_notpos=find(~(brain_new=='M1C'|brain_new=='V1C'|brain_new=='VFC'|brain_new=='HIP'|brain_new=='AMY'|brain_new=='DFC'));

%%
load('D:\study\sub2\brain_age\brain_span\neurodev_proc.mat');
field = fieldnames(neurodev_proc); % four neurodevelopment processes
for i = 1:length(field)
    i
    proc = field{i};
    gene_proc = getfield(neurodev_proc, proc);
    [isMember, indexInB] = ismember(gene_new, gene_proc);
    genepos1 = find(isMember);
    gene_remain_pos=find(isMember==0);
    gene_domi1=data_new(genepos1,domi_pos);
    gene_nodomi1=data_new(genepos1,domi_notpos);
    Data = zscore( [ gene_domi1'; gene_nodomi1' ] );%region(domi/nodomi)*gene(developmental process gene)
    [ ~, PC ] = pca( Data, 'Centered', false );
    if corr( PC( :, 1 ), mean( Data, 2 ) ) < 0
        PC = -PC;
    end
    PC = PC(:,1);%region*1
    PC_normalized = (PC- min(PC)) ./ (max(PC) - min(PC));
    % for trajactory
    domAge = log2(age_new_day(domi_pos));%是domi 脑区的那些年龄
    notdomAge = log2(age_new_day(domi_notpos));%不是domi 脑区的那些年龄
    domPC_minmax = PC_normalized(1:size(gene_domi1,2));%  Data = zscore( [ gene_domi1'; gene_nodomi1' ] );%region*gene
    notdomPC_minmax  = PC_normalized((size(gene_domi1,2)+1):end);
    dataTable = table(domAge, domPC_minmax, 'VariableNames', {'age', 'y'});
    filename = strcat('D:\study\sub2\brain_age\brain_span\',proc,'_',condi,'_domi_cortical.xlsx');
    writetable(dataTable, filename);
    dataTable = table(notdomAge, notdomPC_minmax, 'VariableNames', {'age', 'y'});
    filename = strcat('D:\study\sub2\brain_age\brain_span\',proc,'_',condi,'_nodomi_cortical.xlsx');
    writetable(dataTable, filename);

    % for significant analyze
    domPC = PC(1:length(domAge));
    notdomPC = PC((length(domAge)+1):end);
    Index_dom_6_14 = find(domAge >= log2(4281) & domAge <= log2(14866));
    Index_notdom_6_14 = find(notdomAge >=log2(4281) & notdomAge <=log2(14866));%here
    domPC_6_14 = mean(domPC(Index_dom_6_14));
    notdomPC_6_14 = mean(notdomPC(Index_notdom_6_14));
    diff_PC(i) = domPC_6_14 - notdomPC_6_14;
    nperm = 1000; % permutation number
    for j = 1 : nperm
        j
        myresample = randsample(length(gene_remain_pos),length(genepos1));
        myresample=gene_remain_pos(myresample);
        domGene_perm = data_new(myresample,domi_pos);
        notdomGene_perm = data_new(myresample,domi_notpos);
        Dataperm = zscore( [ domGene_perm'; notdomGene_perm' ] );
        [ ~, PCperm ] = pca( Dataperm, 'Centered', false );
        if corr( PCperm( :, 1 ), mean( Dataperm, 2 ) ) < 0
            PCperm = -PCperm;
        end
        domPCperm = PCperm(1:length(domAge));
        notdomPCperm = PCperm((length(domAge) + 1):end);
        domPCperm_6_14 = mean(domPCperm(Index_dom_6_14));
        notdomPCperm_6_14 = mean(notdomPCperm(Index_notdom_6_14));
        diff_PCperm(j,i) = domPCperm_6_14 - notdomPCperm_6_14;
    end
    if diff_PC(i) > 0
        p_perm(i) = length(find(diff_PCperm(:,i) > diff_PC(i))) / nperm;
    else
        p_perm(i) = length(find(diff_PCperm(:,i) < diff_PC(i))) / nperm;
    end
end
Myelin=diff_PCperm(:,1);
Axon=diff_PCperm(:,2);
Dendrite=diff_PCperm(:,3);
Synapse=diff_PCperm(:,6);

%save('D:\study\sub2\draw\permut_all.mat','Myelin','Axon','Dendrite','Synapse')
% save('D:\study\sub2\draw\permut_sub1_cortical.mat','Myelin','Axon','Dendrite','Synapse')
save('D:\study\sub2\draw\permut_sub1.mat','Myelin','Axon','Dendrite','Synapse')
