
% get sbj list
%sbjListFile = '/home/cxpang/matlab/hcp_sbjLst.txt';
%disp('Get subject list...');
%fid = fopen(sbjListFile,'r');
%sbjList = textscan(fid,'%s');
%fclose(fid);
load hcplist
for si=1:length(data)

    
fileName1=strcat('/HeLabData2/cxpang/fmri/',num2str(data(si)),'_rfMRI_REST1_LR_Atlas_MSMAll.L.fs5.func.gii');
    
fileName2=strcat('/HeLabData2/cxpang/fmri/',num2str(data(si)),'_rfMRI_REST1_LR_Atlas_MSMAll.R.fs5.func.gii');
        Job_Name = ['idemo_' num2str(si)];
        pipeline.(Job_Name).command = 'my_cifti(opt.para1,opt.para2,opt.para3)';
        pipeline.(Job_Name).opt.para1 = fileName1;
        pipeline.(Job_Name).opt.para2 = fileName2;
        pipeline.(Job_Name).opt.para3 = si;
        disp([num2str(si),'. ',fileName]);
        %cii = ciftiopen(fileName,wbPath);   % contain cii.cdata: vNum
        %Xtnum
        %cii=ft_read_cifti(fileName);
end
Pipeline_opt.mode = 'qsub';
Pipeline_opt.qsub_options = '-q all.q,basic.q';
Pipeline_opt.mode_pipeline_manager = 'batch';
Pipeline_opt.max_queued = 1000;
Pipeline_opt.flag_verbose = 1;
Pipeline_opt.flag_pause = 0;
Pipeline_opt.path_logs = ['/home/cxpang'];

psom_run_pipeline(pipeline, Pipeline_opt);