function deployFuncInit_vol(random_number,rand_num,every_time_sub,path,filename,list_cell,maskFile,prepDataFile,outDir,spaR,vxI,ard,iterNum,tNum,alpha,beta,resId)%,rand_num)

if isdeployed
    spaR = str2double(spaR);
    vxI = str2double(vxI);
    ard = str2double(ard);
    iterNum = str2double(iterNum);
    tNum = str2double(tNum);
    alpha = str2double(alpha);
    beta = str2double(beta);
end

if ~exist(prepDataFile,'file')
    maskNii = load_untouch_nii(maskFile);
    maskMat = int32(maskNii.img~=0);
    gNb = constructW_vol(maskMat,spaR);
    
    save(prepDataFile,'gNb','-v7.3');
else
    load(prepDataFile); % containing gNb
end

nmVec = zeros(length(gNb),1);
for gni=1:length(gNb)
    nmVec(gni) = length(gNb{gni});
end
nM = median(nmVec);
Opt.mode = 'qsub'; % qsub session
Opt.max_queued =30;% number
Opt.flag_verbose = true;
Opt.flag_pause = false;
Opt.flag_update = false;
Opt.time_between_checks=30;
Opt.path_logs=fullfile(pwd, 'gsrlog')
Opt.qsub_options = '-q  he_queue.q ';
for j=18 % this step can change to many networks to select a best number of network
    
    
    numUsed = every_time_sub;
    %  pL = round((beta*tNum*numUsed)/(network_num*nM));
    % my_func_initialization_woLoadSrc(sbjData,prepDataFile,outDir,resId,numUsed,network_num,pS,pL,spaR,vxI,ard,iterNum)
    
    for i=1:random_number
        pS = round((alpha*tNum*numUsed)/j);
        pL = round((beta*tNum*numUsed)/(j*nM));
       % my_func_initialization_woLoadSrc(rand_num(i,:),path,filename,list_cell,maskFile,prepDataFile,outDir,resId,numUsed,j,pS,pL,spaR,vxI,ard,iterNum,i);
       
        Pl.(sprintf('Task%d_%d', j,i)).command=...
            'my_func_initialization_woLoadSrc(opt.rand,opt.path,opt.filename,opt.list_cell,opt.maskFile,opt.prefile,opt.outdir,opt.resid,opt.numused,opt.k,opt.ps,opt.pl,opt.spar,opt.vi,opt.Ard,opt.iternum,opt.num)';
        Pl.(sprintf('Task%d_%d', j,i)).opt.rand=rand_num(i,:);
        Pl.(sprintf('Task%d_%d', j,i)).opt.path=path;
        Pl.(sprintf('Task%d_%d', j,i)).opt.filename=filename;
        Pl.(sprintf('Task%d_%d', j,i)).opt.list_cell=list_cell;
        Pl.(sprintf('Task%d_%d', j,i)).opt.maskFile=maskFile;
        Pl.(sprintf('Task%d_%d', j,i)).opt.prefile=prepDataFile;
        Pl.(sprintf('Task%d_%d', j,i)).opt.outdir=outDir;
        Pl.(sprintf('Task%d_%d', j,i)).opt.resid=resId;
        Pl.(sprintf('Task%d_%d', j,i)).opt.numused=numUsed;
        Pl.(sprintf('Task%d_%d', j,i)).opt.k=j;
        Pl.(sprintf('Task%d_%d', j,i)).opt.ps=pS;
        Pl.(sprintf('Task%d_%d', j,i)).opt.pl=pL;
        Pl.(sprintf('Task%d_%d', j,i)).opt.spar=spaR;
        Pl.(sprintf('Task%d_%d', j,i)).opt.vi=vxI;
        Pl.(sprintf('Task%d_%d', j,i)).opt.Ard=ard;
        Pl.(sprintf('Task%d_%d', j,i)).opt.iternum=iterNum;
        Pl.(sprintf('Task%d_%d', j,i)).opt.num=i;
        
    end
end
psom_run_pipeline(Pl, Opt);