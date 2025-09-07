function vali_deploy(path,sbj_num,maskFile,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn)


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
sbjNum = sbj_num;



Opt.mode = 'qsub'; % qsub session
Opt.max_queued =100;% number
Opt.flag_verbose = true;
Opt.flag_pause = false;
Opt.flag_update = false;
Opt.time_between_checks=30;
Opt.path_logs=fullfile(pwd, 'nmf2_log')%%change
Opt.qsub_options = '-q  public_queue.q '; 

for si=1:sbjNum%%
    si
   %sbjData = prepareFuncData_vol_func_single(maskMat,list_cell{si});
    %^sbjNum = length(sbjData);
   % display('reading')
   if(si<10)
       path_new=strcat(path,'sub-000',num2str(si),'/');
       file_name=strcat('/xbcNGSdswransub-000',num2str(si),'.nii');%for dida dataset
   elseif(si<100)
       path_new=strcat(path,'sub-00',num2str(si),'/');
       file_name=strcat('/xbcNGSdswransub-00',num2str(si),'.nii');%for dida dataset
   elseif(si<1000)
       path_new=strcat(path,'sub-0',num2str(si),'/');
       file_name=strcat('/xbcNGSdswransub-0',num2str(si),'.nii');%for dida dataset
   else
       path_new=strcat(path,'sub-',num2str(si),'/');
       file_name=strcat('/xbcNGSdswransub-',num2str(si),'.nii');%for dida dataset
   end
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    % vali_func(path,file_name,si,maskMat,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,spaR,vxI,ard,eta,iterNum,calcGrp,parforOn);
     Pl.(sprintf('Task%d', si)).command=...
    'vali_func(opt.path,opt.file_name,opt.si,opt.mask,opt.prefile,opt.outdir,opt.resid,opt.initName,opt.K,opt.alphaS21,opt.alphaL,opt.spar,opt.vi,opt.Ard,opt.eta,opt.iternum,opt.calcGrp,opt.parforOn)';
     Pl.(sprintf('Task%d', si)).opt.path=path_new;
      Pl.(sprintf('Task%d', si)).opt.file_name=file_name;
    Pl.(sprintf('Task%d', si)).opt.si=si;
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
