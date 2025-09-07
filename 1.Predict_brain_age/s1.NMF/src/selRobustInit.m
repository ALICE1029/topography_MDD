function initV = selRobustInit(inFile,K,outDir,network_num)
% select robust initialization
% inFile is a text file including paths to all the candidate initializations, each line for one initialization

%ncutPath = [gdPath,filesep,'Ncut_9'];
%addpath(genpath(ncutPath));

%
fid = fopen(inFile, 'r');
initList = textscan(fid, '%s');
initList = initList{1};
fclose(fid);

repNum = length(initList);
disp('select best V...');

% load results
resSet = [];
resSet2 = [];
for ri=1:repNum
    resName = initList{ri};
    load(resName);

    resSet = [resSet, initV];
    resSet2=[resSet2, initU];
    clear initV;
end
resSet = resSet';

% clustering by ncut
corrVal = corr(resSet');
corrVal(isnan(corrVal)) = -1;
nDis = 1 - corrVal;
triuInd = triu(ones(size(nDis)),1);
nDisVec = nDis(triuInd==1);

nW = exp(-nDis.^2 ./ (median(nDisVec).^2));
nW(isnan(nW)) = 0;

sumW = sum(nW,1);
sumW(sumW==0) = 1;
D = diag(sumW);
L = sqrt(inv(D))*nW*sqrt(inv(D));  
L = (L+L')/2;
opts.disp = 0;
[Ev,~] = eigs(double(L),K,'LA',opts);
normvect = sqrt(diag(Ev*Ev'));
normvect(normvect==0.0) = 1;
Ev = diag(normvect) \ Ev;

[EvDiscrete,~] = discretisation(Ev);
EvDiscrete = full(EvDiscrete);
[~, C] = max(EvDiscrete,[],2);

% get centroid
initV = zeros(size(resSet,2),K);
for ki=1:K            
    % % typical point
    if sum(C==ki)>1
        pos=find(C==ki);
        candSet = resSet(C==ki,:)';
        corrW = abs(corr(candSet));
        corrW(isnan(corrW)) = 0;
        [mVal, mInd] = max(sum(corrW));
        times=pos(mInd)/network_num+1;
         resSet2=[resSet2, initU];
          initU(:,ki) = resSet2(:,pos(mInd));
        initV(:,ki) = candSet(:,mInd);
    elseif sum(C==ki)==1
        initV(:,ki) = resSet(C==ki,:);
         initU(:,ki) =resSet(C==ki,:);
    end
    
    
    
end

initV = initV ./ max(eps,repmat(max(initV),size(initV,1),1));

if ~exist(outDir,'dir')
    mkdir(outDir);
end
save([outDir,filesep,'init.mat'],'initV','initU');
