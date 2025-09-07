maskfile ='/HeLabData2/cxpang/fmri/bnu/mask/GMMask_3mm.nii';
hdr_mask = spm_vol(maskfile);
[vol_mask,XYZmm] = spm_read_vols(hdr_mask);
vol_mask=vol_mask(:);
for i=1:271633
    if(vol_mask(i)==0)
        XYZmm(:,i)=zeros(3,1);
    end
end
XYZmm(:,all(XYZmm==0,1)) = [];

