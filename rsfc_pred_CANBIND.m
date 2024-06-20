clear
%% 1. %% %% import canbind data
cd C:\Users\peter\Documents\CANBIND
demo=readtable('data\tabular\Clinical_DEMO_MDD_BSLN\CBN01_DEMO_DATA_Z3_01_V01_TRTMT.csv');
SHAPS=readtable('data\tabular\Clinical_SHAPS_MDD_BSLN\CBN01_SHAPS_DATA_Z3_01_V01_TRTMT.csv'); NEO=readtable('data\tabular\Clinical_NEOFFI_MDD_BSLN\CBN01_NEOFFI_DATA_Z3_01_V01_TRTMT.csv');
MADRS_BL=readtable('data\tabular\Clinical_MADRS_MDD_BSLN\CBN01_MADRS_DATA_Z3_01_V01_TRTMT.csv');
MADRS=readtable('data\tabular\Clinical_MADRS_MDD_PHASE1\CBN01_MADRS_DATA_Z3_02_V01_TRTMT.csv');
MADRS(strcmp(MADRS.EVENTNAME,'Week 2') & strcmp(MADRS.VISITSTATUS, 'END'),:)=[]; 
ID=MADRS.SUBJLABEL(strcmp(MADRS.EVENTNAME,'Week 4') & strcmp(MADRS.VISITSTATUS, 'END'),:); for i=1:length(ID); MADRS(strcmp(MADRS.SUBJLABEL, ID{i}),:)=[]; end
%MADRS(strcmp(MADRS.EVENTNAME,'Week 6') & strcmp(MADRS.VISITSTATUS, 'END'),:)=[];
BMI=readtable('data\tabular\Clinical_BMI_MDD_BSLN\CBN01_BMI_DATA_Z3_01_V01_TRTMT.csv');
GAD7=readtable('data\tabular\Clinical_GAD7_MDD_BSLN\CBN01_GAD7_DATA_Z3_01_V01_TRTMT.csv');
subjects_madrs=unique(MADRS.SUBJLABEL);
MADRS_FU=table();
for i=1:length(subjects_madrs)
    s=subjects_madrs{i};
    d=MADRS(strcmp(MADRS.SUBJLABEL, s),:);
    MADRS_FU.SUBJLABEL(i)=subjects_madrs(i);
    MADRS_FU.MADRS_TOT_final_score(i)=d.MADRS_TOT_PRO_RATED(length(d.SUBJLABEL));
    MADRS_FU.MADRS_2wks(i)=d.MADRS_TOT_PRO_RATED(1);
end
demo_CANBIND=innerjoin(demo, MADRS_BL);demo_CANBIND=innerjoin(demo_CANBIND, MADRS_FU);demo_CANBIND=innerjoin(demo_CANBIND, SHAPS, 'Keys','SUBJLABEL'); demo_CANBIND=innerjoin(demo_CANBIND, NEO, 'Keys','SUBJLABEL'); demo_CANBIND=innerjoin(demo_CANBIND, BMI, 'Keys','SUBJLABEL');demo_CANBIND=innerjoin(demo_CANBIND, GAD7, 'Keys','SUBJLABEL');
demo_CANBIND.response=(demo_CANBIND.MADRS_TOT_final_score./demo_CANBIND.MADRS_TOT_PRO_RATED<0.5); 
demo_CANBIND.remission=demo_CANBIND.MADRS_TOT_final_score<=9; 
demo_CANBIND(isnan(demo_CANBIND.AGE), :)=[];
nansum(demo.SEX==2)
cd C:\Users\peter\Documents\CANBIND\reports\
load("globalFC_canbind.mat");  %load('globalsubFC_embarc.mat'); l=reshape(globalsubFC(:,4,:), [157 180]);r=reshape(globalsubFC(:,9,:), [157 180]); globalfFC=[l,r];
      roi=[  179 180 359 360  ]; for i=1:157; globalfFC(i,:)=nanmean(Fnetmats_all(i,roi,:)); end % 221   252   270 241 %   68    90    98   108   165   221   240   241   245   252   270   345
      %roi=[ 179 180 359 360  ]; %108 114 %0 35 64 66 76 91 93   115   146   161   165   180   197   216   241   253   272   283   292   307 322   335   360
       %not bad for placebo:  [ 61 64 241 244 ] : 64 66 76      % OR 61 180 241 360      % OR 61 179 180 241 359 360        %  62 179 180 242 359 360     % 61 62 179 180 241  242 359 360            35    64    66    91    93   108   115   161   180  197   216   272   283   307   335   360
tmp=subjects;       
subjects=table();globalfFC(globalfFC(:,1)==0,:)=NaN;
subjects.rsfmri=globalfFC; subjects.meanfd=mean_fd';
subjects.SUBJLABEL=tmp;
subcort_rois={'L_accumbens', 'L_amygdala', 'L_caudate', 'L_hippocampus', 'L_putamen','R_accumbens', 'R_amygdala', 'R_caudate', 'R_hippocampus', 'R_putamen'};
subjects.SUBJLABEL=replace(subjects.SUBJLABEL,'sub-','');

cort_rois=readtable('C:\Users\peter\Documents\GABA\HCPMMP1_on_MNI152_ICBM2009a_nlin.txt');
%
demo_CANBIND.SUBJLABEL=replace(demo_CANBIND.SUBJLABEL,'_','');
demo_CANBIND=innerjoin(demo_CANBIND, subjects);
demo_CANBIND.site(1:8)=1;demo_CANBIND.site(9:37)=2;demo_CANBIND.site(38:48)=3;demo_CANBIND.site(49:73)=4; demo_CANBIND.site(74:122)=5;demo_CANBIND.site(123:length(demo_CANBIND.SUBJLABEL))=6;
demo_CANBIND.SHAPS=nansum(demo_CANBIND{:,48:61}')'; demo_CANBIND.SHAPS(demo_CANBIND.SHAPS==0)=NaN; %demo_CANBIND.SHAPS(isnan(demo_CANBIND.SHAPS))=nanmean(demo_CANBIND.SHAPS);
%% 2. %%%% %%%% %% import embarc data 
cd C:\Users\peter\Documents\EMBARC\reports
load('globalFC_embarc.mat');  % load('globalsubFC_embarc.mat');l=reshape(globalsubFC(:,4,:), [330 180]);r=reshape(globalsubFC(:,9,:), [330 180]); globalfFC=[l,r];
       for i=1:330; globalfFC(i,:)=nanmean(Fnetmats_all(i,roi,:)); end % 221   252   270 241 %   68    90    98   108   165   221   240   241   245   252   270   345
%108 114 FOP3 FOP4; 240   241 24dv p32pr         L_d23ab_ROI   L_s32_ROI
%0 35 64 66 76 91 93   115   146   161   165   180   197   216   241   253   272   283   292   307 322   335   360
%0    90   108   114   126   210   221   229   240   252   341
tmp=subjects;
subjects=table();globalfFC(globalfFC(:,1)==0,:)=NaN;
subjects.rsfmri=globalfFC; subjects.meanfd=mean_fd';
subjects.ProjectSpecificID=tmp;
subcort_rois={'L_accumbens', 'L_amygdala', 'L_caudate', 'L_hippocampus', 'L_putamen','R_accumbens', 'R_amygdala', 'R_caudate', 'R_hippocampus', 'R_putamen'};

clindata=readtable('C:\Users\peter\Documents\EMBARC\data\clinical_data_fromChristianW.xlsx');clindata(clindata.delete==1, :)=[];
cort_rois=readtable('C:\Users\peter\Documents\GABA\HCPMMP1_on_MNI152_ICBM2009a_nlin.txt');
%cort_rois(181:360,:)=[];

clin_ct=innerjoin(subjects, clindata);clin_ct.remission_step1=clin_ct.w8_score_17<8;clin_ct.remission_step1=double(clin_ct.remission_step1);clin_ct.remission_step1(isnan(clin_ct.w8_score_17))=NaN;
clin_ct.remission_step2=clin_ct.w16_score_17<8;clin_ct.remission_step2=double(clin_ct.remission_step2);clin_ct.remission_step2(isnan(clin_ct.w16_score_17))=NaN;
clin_ct.strf_14(clin_ct.strf_14<13)=NaN;  clin_ct.strf_14(clin_ct.strf_14>55)=NaN;  clin_ct.strf_waist(clin_ct.strf_waist>70)=NaN;  
beta=regress(clin_ct.strf_14, [ones(length(clin_ct.strf_waist),1), clin_ct.strf_waist]); imputed_bmi=beta'.*[ones(length(clin_ct.strf_waist),1), clin_ct.strf_waist]; imputed_bmi=sum(imputed_bmi')';
clin_ct.strf_14(isnan(clin_ct.strf_14))=imputed_bmi(isnan(clin_ct.strf_14));% clin_ct.strf_14(isnan(clin_ct.strf_14))=nanmean(clin_ct.strf_14);
clin_ct.education=clin_ct.demo_educa_years;
%clin_ct(strcmp(clin_ct.Stage1TX, 'PLA') & strcmp(clin_ct.Stage2TX, 'PLA'),:)=[]; clin_ct(strcmp(clin_ct.Stage1TX, 'PLA') & strcmp(clin_ct.Stage2TX, ''),:)=[]; clin_ct(strcmp(clin_ct.Stage1TX, 'PLA') & strcmp(clin_ct.Stage2TX, 'BUP'),:)=[];
%clin_ct.remission_step1(strcmp(clin_ct.Stage1TX, 'PLA') & strcmp(clin_ct.Stage2TX, 'SER'),:)=clin_ct.remission_step2(strcmp(clin_ct.Stage1TX, 'PLA') & strcmp(clin_ct.Stage2TX, 'SER'),:);
%clin_ct.w8_responder(isnan(clin_ct.w8_score_17) & clin_ct.w6_score_17./clin_ct.w0_score_17<0.5)={'Yes'}; clin_ct.w8_responder(isnan(clin_ct.w8_score_17) & clin_ct.w6_score_17./clin_ct.w0_score_17>0.5)={'No'};
%%%%%%%%%%%%% quick combat for EMBARC
clin_ct(isnan(clin_ct.rsfmri(:,1)),:)=[];
batch = clin_ct.siteNum';
dat = clin_ct.rsfmri';
age = clin_ct.age; sex = clin_ct.sex;
mod = [age sex];
full_harmonized = combat(dat, batch, mod, 1); full_harmonized=full_harmonized';
%clin_ct.rsfmri=full_harmonized;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clin_stage2=clin_ct(strcmp(clin_ct.Stage1TX, 'PLA') & strcmp(clin_ct.Stage2TX, 'SER'),:);
clin_pla=clin_ct(strcmp(clin_ct.Stage1TX, 'PLA'),:);
clin_ct(strcmp(clin_ct.Stage1TX, 'PLA'),:)=[];
%clin_ct(~strcmp(clin_ct.Stage2TX, 'BUP'),:)=[];
nansum(clin_ct.remission_step2)/length(clin_ct.remission_step2)
nansum(clin_ct.remission_step1)/length(clin_ct.remission_step2)
sum(strcmp(clin_ct.w8_responder,'Yes'))/length(clin_ct.remission_step2)
stai=zscore(clin_ct.stai_pre_final_score); bmi=clin_ct.strf_14; masq=clin_ct.masq2_score_gd+clin_ct.masq2_score_ad+clin_ct.masq2_score_aa; shaps=clin_ct.shaps_total_continuous;sex=clin_ct.sex;employ=clin_ct.demo_employ_status;clin_ct(:, 118:1191)=[];clin_ct.sex=sex;clin_ct.shaps=shaps;clin_ct.employ=employ;clin_ct.masq=masq;clin_ct.bmi=bmi;clin_ct.stai=stai;
stai=zscore(clin_stage2.stai_pre_final_score); bmi=clin_stage2.strf_14; masq=clin_stage2.masq2_score_gd+clin_stage2.masq2_score_ad+clin_stage2.masq2_score_aa; shaps=clin_stage2.shaps_total_continuous;sex=clin_stage2.sex;employ=clin_stage2.demo_employ_status;clin_stage2(:, 118:1191)=[];clin_stage2.sex=sex;clin_stage2.shaps=shaps;clin_stage2.employ=employ;clin_stage2.masq=masq;clin_stage2.bmi=bmi;clin_stage2.stai=stai;
stai=zscore(clin_pla.stai_pre_final_score); bmi=clin_pla.strf_14; masq=clin_pla.masq2_score_gd+clin_pla.masq2_score_ad+clin_pla.masq2_score_aa; shaps=clin_pla.shaps_total_continuous;sex=clin_pla.sex;employ=clin_pla.demo_employ_status;clin_pla(:, 118:1191)=[];clin_pla.sex=sex;clin_pla.shaps=shaps;clin_pla.employ=employ;clin_pla.masq=masq; clin_pla.bmi=bmi;clin_pla.stai=stai;
cd C:\Users\peter\Documents\CANBIND\reports\

%%  %% %% %%%% %%%% 3. ComBat %% %% %%%% %%%% 
batch = demo_CANBIND.site';
dat = demo_CANBIND.rsfmri';
age = demo_CANBIND.AGE; sex = demo_CANBIND.SEX;sex = dummyvar(sex);
mod = [age sex(:,2)];
full_harmonized = combat(dat, batch, mod, 1); full_harmonized=full_harmonized';
%demo_CANBIND.rsfmri=full_harmonized;

clin_ct(isnan(clin_ct.rsfmri(:,1)),:)=[];
batch = [demo_CANBIND.site', clin_ct.siteNum'+6, clin_stage2.siteNum'+6];
dat = [demo_CANBIND.rsfmri;  clin_ct.rsfmri; clin_stage2.rsfmri]';
age = [demo_CANBIND.AGE; clin_ct.age; clin_stage2.age]; sex = [demo_CANBIND.SEX==1; clin_ct.sex; clin_stage2.sex];
mod = [age sex];
full_harmonized = combat(dat, batch, mod, 1); full_harmonized=full_harmonized';

%    demo_CANBIND.rsfmri=full_harmonized(1:145,:);
%    clin_ct.rsfmri=full_harmonized(146:146+140,:);
%    clin_stage2.rsfmri=full_harmonized(287:359,:);

%% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %%
%% 4. RUN THE TRAINING MODEL IN CANBIND %% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %%
%% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %%    
Xlogist=[demo_CANBIND.rsfmri,  demo_CANBIND.SHAPS, demo_CANBIND.AGE, demo_CANBIND.SEX, madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1 , demo_CANBIND.BL_BMI  ];%  , demo_CANBIND.lh_hippo, demo_CANBIND.rh_hippo];%    ];%  ];%, demo_CANBIND.GAD7_Tot];% 
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 
Ylogist=demo_CANBIND.response==1; 
%Ylogist=demo_CANBIND.remission;
%Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 
%       Ylogist(isnan(Xlogist(:,1)) )=[]; Xlogist(isnan(Xlogist(:,1) ),:)=[]; Xlogist(:,roi)=[];varnames(roi)=[];

[B,FitInfo] = lassoglm(Xlogist,Ylogist,'binomial','CV',10, 'PredictorNames', varnames  ,'alpha', 0.3); %0.01 with the base only
%lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = glmval(coef,Xlogist,'logit');[X,Y,T,AUC] = perfcurve(Ylogist,yhat, 1);AUC

%Xlogist=Xlogist(:,B(:,idxLambda1SE)~=0);varnames=varnames(B(:,idxLambda1SE)~=0); %pasimonious model test > set alpha to very low 0.001
save('test_coef.mat','coef')
test_can_in_embarc

%test_embarc_in_can

%% figure out the rois
%for i=1:360; mdl=fitlm(demo_CANBIND.response, demo_CANBIND.rsfmri(:,i));pval(i)=mdl.Coefficients.pValue(2);end
%min(pval); find(pval<0.005)
%cort_rois.Var2(pval<0.05)'   


%find(strcmp(cort_rois.Var2,'L_9m_ROI'))
roinames=varnames(coef(2:length(coef))>0); roinames=replace(roinames,'L_','lh.L_');  roinames=replace(roinames,'R_','rh.R_'); roinames=replace(roinames,'_ROI','_ROI.label');

roinames=varnames(coef(2:length(coef))<0);roinames=replace(roinames,'L_','lh.L_');  roinames=replace(roinames,'R_','rh.R_'); roinames=replace(roinames,'_ROI','_ROI.label');
