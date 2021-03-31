function [counts_ccf,areas_ccf,counts_fp,areas_fp,counts_ccf2fp,areas_ccf2fp]=AreaCounts_ccf2fp(P,atlas_resolution,bregma,CCFtoFPtable, FPtable, chon_images_loc,area_acr)

% pulls out areas associated with user defined area grouping in CCFv3 and FP labels and corresponding cell counts
% area_acr contains list of strings for grouping areas

if nargin<12
    area_acr={'AUD','MO','SS','ACA','CLA','TEa','RSPagl','RSPd','RSPv','ECT','PRE','POST','PAR','ENTl','ENTm','VISpm','VISam','Anterior area','VISrl','VISal','Lateral visual area','VISli','VISpor','VISpl'};

% use AP, DV, and ML coordinates to get the point's position in the atlas
    curr_objectPoints= NaN(length(P.ap),3);
    curr_objectPoints(:,1)= round(-P.ap/atlas_resolution + bregma(1));
    curr_objectPoints(:,2)= round(P.dv/atlas_resolution + bregma(2));
    curr_objectPoints(:,3)= round(P.ml/atlas_resolution + bregma(3));
    
 for point = 1:size(curr_objectPoints,1)
     
    [Ann_FP, Name_FP, Acr_FP] = CCF_to_FP(curr_objectPoints(point,1), curr_objectPoints(point,2), curr_objectPoints(point,3), CCFtoFPtable, FPtable, chon_images_loc);
        ann_FP{point} = Ann_FP;
        name_FP{point} = Name_FP;
        acr_FP{point} = Acr_FP;    
 end
 

    roi_ccf_fp = table(P.name,P.acr, name_FP', acr_FP', P.ap, P.dv, P.ml, P.ann, ann_FP',...
         'VariableNames', {'CCF_name', 'CCF_acronym', 'FP_name', 'FP_acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex', 'fpIndex'});
%      disp(roi_ccf_fp)

Areas_ccf=roi_ccf_fp.CCF_acronym;
Areas_fp=roi_ccf_fp.FP_acronym;

AreaCounts_ccf=categorical(Areas_ccf);
AreaCounts_fp=categorical(Areas_fp);

% figure; 
% subplot(1,2,1)
% histogram(AreaCounts_ccf,'DisplayOrder', 'ascend','Orientation','horizontal');%'DisplayOrder', 'ascend',
% title('Allen CCFv3 labels')
% 
% subplot(1,2,2)
% histogram(AreaCounts_fp,'DisplayOrder', 'ascend','Orientation','horizontal');
% title('Franklin & Paxinos labels (Chon et al. 2020)')

AC_ccf=nominal(AreaCounts_ccf);
AC_fp=nominal(AreaCounts_fp);

[counts_ccf,areas_ccf] = summary(AC_ccf);
[counts_fp,areas_fp] = summary(AC_fp);

for i=1:length(area_acr)
    if strcmp(area_acr{i},'Anterior area')|| strcmp(area_acr{i},'Lateral visual area')
        allen_idx(:,i)=contains(roi_ccf_fp.CCF_name,area_acr{i});
    else
        allen_idx(:,i)=contains(roi_ccf_fp.CCF_acronym,area_acr{i});
    end
end

for i=1:length(area_acr)
allen_labels{i}= roi_ccf_fp.CCF_acronym(allen_idx(:,i));
fp_labels{i}= roi_ccf_fp.FP_acronym(allen_idx(:,i));
end

for i=1:length(area_acr)
[counts_ccf2fp{i},areas_ccf2fp{i}]=summary(nominal(categorical(fp_labels{i})));
end
end