function my_func_initialization_woLoadSrc(rand_num,path,filename,list_cell,maskFile,prepData,outDir,resId,numUsed,K,pS,pL,spaR,vxI,ard,iterNum,random_number)
list_cell_new=cell(numUsed,1);
for i=1:numUsed
    list_cell_new{i}=list_cell{rand_num(i)};
end
list_cell=list_cell_new;
sbjData=prepareFuncData_vol_func(path,filename,list_cell,maskFile);
if isdeployed
    numUsed = str2double(numUsed);
    K = str2double(K);
    pS = str2double(pS);
    pL = str2double(pL);
    spaR = str2double(spaR);
    vxI = str2double(vxI);
    ard = str2double(ard);
    iterNum = str2double(iterNum);
end

% parameter setting
options = [];
options.maxIter = iterNum;
options.error = 1e-6;
options.nRepeat = 1;
options.minIter = 100; %iterNum;
options.Converge = 1;
options.meanFitRatio = 0.1;

options.S1 = pS;
options.L = pL;
options.robust = 0; % robust NMF

neiR = spaR;
vxlInfo = vxI;

if ard~=0
    options.ard = 1;
end

% output name  INCLUDE PARAMETER INFORMATION
resDir = [outDir,filesep,resId,'boostrap_comp',num2str(K)];
if ~exist(resDir,'dir')
    mkdir(resDir);
end

% load preparation data, it contains gNb
load(prepData);

% load data, it contains sbjData
disp('load data...');

nSmp = size(sbjData{1},2);
totFea = 0;
mFea = zeros(numUsed,1);

%concentrate all subjects:obtain catx
for si=1:length(sbjData)
    mFea(si) = size(sbjData{si},1);
    totFea = totFea + mFea(si);
end
catX = zeros(totFea,nSmp,'single');

sInd = 1;
for si=1:length(sbjData)
%     fprintf('  sbj%d->',si);
%     if mod(si,15)==0
%         fprintf('\n');
%     end
    
    origSbjData = sbjData{si};
    origSbjData = dataPrepro(origSbjData,'vp','vmax');%delete negative value+normalization

    nanSbj = isnan(origSbjData);
    if sum(nanSbj(:))>0
        error([' nan exists: ','sbj',num2str(si)]);
    end

    eInd = sInd + mFea(si) - 1;
    catX(sInd:eInd,:) = origSbjData;
    sInd = eInd + 1;
end
clear sbjData;

fprintf('\n');


load('/HeLabData2/cxpang/fmri/bnu/s0/test_CreatePrepData') % containing gNb
vNum = size(gNb,1);  
catW = sparse(vNum,vNum);
for vi=1:vNum
    for ni=1:length(gNb{vi})
        nei = gNb{vi}(ni);
        if vxlInfo~=0%
            if vi<nei
                corrVal = (1+corr(catX(:,vi),catX(:,nei)))/2;
            else
                continue;
            end
        else
            corrVal = 1;
        end
        if isnan(corrVal)
            corrVal = 0;
        end
        catW(vi,nei) = corrVal;
        catW(nei,vi) = corrVal;
    end
end
disp('snmf...');
iU = [];
iV = [];
[initU, initV] = mNMF_sp(catX, K, catW, options, iU, iV);
save([resDir,filesep,strcat('init',num2str(random_number),'.mat')],'initV','initU','-v7.3');
end


% if isdeployed
%     exit;
% else
%     disp('Done!');
% end


