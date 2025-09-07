function deployFuncMvnmfL21p1_func_vol(path,file_name,list_cell,maskFile,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn)


if isdeployed
    K = str2double(K);
    alphaS21 = str2double(alphaS21);
    alphaL = str2double(alphaL);
    vxI = str2double(vxI);
    spaR = str2double(spaR);
    ard = str2double(ard);
    eta = str2double(eta);
    iterNum = str2double(iterNum);
    calcGrp = str2double(calcGrp);
    parforOn = str2double(parforOn);
end
maskNii = load_untouch_nii(maskFile);
maskMat = maskNii.img~=0;
sbjNum = length(list_cell);



Opt.mode = 'qsub'; % qsub session
Opt.max_queued =50;% number
Opt.flag_verbose = true;
Opt.flag_pause = false;
Opt.flag_update = false;
Opt.time_between_checks=30;
Opt.path_logs=fullfile(pwd, 'nmff_log')%%change
Opt.qsub_options = '-q  he_queue.q '; 

for si=1:sbjNum%%
   %sbjData = prepareFuncData_vol_func_single(maskMat,list_cell{si});
    %^sbjNum = length(sbjData);
   % display('reading')
  
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
     Pl.(sprintf('Task%d', si)).command=...
    'func_mMvNMF4fmri_l21p1_ard_woSrcLoad(opt.path,opt.file_name,opt.si,opt.mask,opt.prefile,opt.outdir,opt.resid,opt.initName,opt.K,opt.alphaS21,opt.alphaL,opt.spar,opt.vi,opt.Ard,opt.eta,opt.iternum,opt.calcGrp,opt.parforOn)';
     Pl.(sprintf('Task%d', si)).opt.path=path;
      Pl.(sprintf('Task%d', si)).opt.file_name=file_name;
    Pl.(sprintf('Task%d', si)).opt.si=list_cell{si};
    Pl.(sprintf('Task%d', si)).opt.mask=maskMat;
    Pl.(sprintf('Task%d', si)).opt.prefile=prepDataFile;
    Pl.(sprintf('Task%d', si)).opt.outdir=outDir;
    Pl.(sprintf('Task%d', si)).opt.resid=resId;
    Pl.(sprintf('Task%d', si)).opt.initName=initName;
    Pl.(sprintf('Task%d', si)).opt.K=K;
    Pl.(sprintf('Task%d', si)).opt.alphaS21=alphaS21;
    Pl.(sprintf('Task%d', si)).opt.alphaL=alphaL;
    Pl.(sprintf('Task%d', si)).opt.spar=spaR;
    Pl.(sprintf('Task%d', si)).opt.vi=vxI;
    Pl.(sprintf('Task%d', si)).opt.Ard=ard;
    Pl.(sprintf('Task%d', si)).opt.eta=eta;
    Pl.(sprintf('Task%d', si)).opt.iternum=iterNum;
    Pl.(sprintf('Task%d', si)).opt.calcGrp=calcGrp;
    Pl.(sprintf('Task%d', si)).opt.parforOn=parforOn;
    
end

psom_run_pipeline(Pl, Opt);

if isdeployed
    exit;
else
    disp('Done!');
end
