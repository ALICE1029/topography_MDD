%read cifti and save as mat


function cii=my_cifti(gii_l,gii_r,med_l,med_r,id)
disp('reading data');
tmp=gifti(gii_l);
fun_l=tmp.cdata;
tmp=gifti(gii_r);
fun_r=tmp.cdata;
fun_l(med_l,:) = [];
fun_r(med_r,:) = [];
fun=[fun_l;fun_r];
cii=fun;

 %concentrate
save(strcat('/HeLabData2/cxpang/fmri/bnu/func_mat/func_',num2str(id),'.mat'),'fun','-v7.3')
%save(strcat('/HeLabData2/cxpang/fmri/bnu/func_mat/func_',(id),'.mat'),'fun','-v7.3')
end