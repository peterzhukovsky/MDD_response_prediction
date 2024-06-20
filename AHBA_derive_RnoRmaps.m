%% AHBA analysis of responder vs non responder differences
%first run rsfc_pred_CANBIND.m fully 
%canbind
for i=1:360; mdl=fitlm([demo_CANBIND.response, demo_CANBIND.meanfd, demo_CANBIND.AGE, demo_CANBIND.SEX], demo_CANBIND.rsfmri(:,i));pval(i)=mdl.Coefficients.pValue(2);tstatCAN(i)=mdl.Coefficients.tStat(2);end
min(pval); find(pval<0.005)
cort_rois.Var2(bhfdr(pval)<0.05)'  
%embarc 1
clin_ct(strcmp(clin_ct.w8_responder, ''), :)=[];
for i=1:360; mdl=fitlm([strcmp(clin_ct.w8_responder, 'Yes'), clin_ct.meanfd, clin_ct.age, clin_ct.sex], clin_ct.rsfmri(:,i));pval(i)=mdl.Coefficients.pValue(2);tstatEMB1(i)=mdl.Coefficients.tStat(2);end
min(pval); find(pval<0.005)
cort_rois.Var2(bhfdr(pval)<0.1)'  
%anova1(clin_ct.meanfd, strcmp(clin_ct.w8_responder, 'Yes'))

%embarc 2
clin_stage2(isnan(clin_stage2.w16_score_17),:)=[]; 

for i=1:360; mdl=fitlm([(clin_stage2.w16_score_17./clin_stage2.w0_score_17<0.5), clin_stage2.meanfd, clin_stage2.age, clin_stage2.sex], clin_stage2.rsfmri(:,i));pval(i)=mdl.Coefficients.pValue(2);tstatEMB2(i)=mdl.Coefficients.tStat(2);end
min(pval); find(pval<0.01)
cort_rois.Var2(bhfdr(pval)<0.1)'  

c3vs1_tstat_old=[tstatCAN', tstatEMB1', tstatEMB2'];

[r,p]=corr(c3vs1_tstat_old)