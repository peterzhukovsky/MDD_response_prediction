function [output] = elasticnet_embarc(Xlogist,Ylogist)

varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];% 

[B,FitInfo] = lassoglm(Xlogist(:, 361:366),Ylogist,'binomial','CV',10, 'PredictorNames', varnames  ,'alpha', 0.001); %     
%lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = glmval(coef,Xlogist(:, 361:366),'logit');[X,Y,T,AUC] = perfcurve(Ylogist,yhat, 1);AUC

%% CANBIND test
load('clin_canbind_data.mat');
XTest=[demo_CANBIND.SHAPS_Tot, demo_CANBIND.AGE, demo_CANBIND.SEX,...
    madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1, demo_CANBIND.BL_BMI ];%   , demo_CANBIND.lh_hippo, demo_CANBIND.rh_hippo];%    
varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
%varnames(ix==0)=[];XTest(:, ix==0)=[];

YTest=demo_CANBIND.response==1; 
yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUC_clinCANBIND] = perfcurve(YTest,yhat, 1);AUC_clinCANBIND %AUC_test(randfold)=AUC;

%figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')

%% EMBARC stage 2 test
load('clin_embarc_data.mat');
XTest=[clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17,...
    clin_stage2.employ<=2, clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %
clin_stage2.w16_score_17(isnan(clin_stage2.w16_score_17))=clin_stage2.w10_score_17(isnan(clin_stage2.w16_score_17));
YTest=(clin_stage2.w16_score_17./clin_stage2.w0_score_17<0.5); %Ytest=clin_ct.remission_step2;
%XTest(:, ix==0)=[];
%YTest=clin_stage2.remission_step2;
XTest(isnan(clin_stage2.w16_score_17),:)=[]; YTest(isnan(clin_stage2.w16_score_17))=[]; 
XTest(isnan(YTest),:)=[]; YTest(isnan(YTest))=[]; 
%        YTest(isnan(XTest(:,1)) )=[]; XTest(isnan(XTest(:,1) ),:)=[]; XTest(:,roi)=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCclin2EMB] = perfcurve(YTest,yhat, 1);AUCclin2EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color','r'); plot([0,1],[0,1],'k'); 

%% EMBARC placebo test
XTest=[clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,clin_pla.employ<=2, clin_pla.bmi];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
YTest=strcmp(clin_pla.w8_responder,'Yes');
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];
yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCclin1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCclin1PLA_EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 













%% Clinical model eval
cort_rois=readtable('C:\Users\peter\Documents\GABA\HCPMMP1_on_MNI152_ICBM2009a_nlin.txt');
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];% 
load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef'); 
coef(1)=[];
ix=coef~=0;
varnames(ix==0)=[];Xlogist(:, ix==0)=[];

%%

[B,FitInfo] = lassoglm(Xlogist,Ylogist,'binomial','CV',10, 'PredictorNames', varnames  ,'alpha', 0.001); %     
%lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = glmval(coef,Xlogist,'logit');[X,Y,T,AUC] = perfcurve(Ylogist,yhat, 1);AUC

%% CANBIND test
load('clin_canbind_data.mat');
XTest=[demo_CANBIND.rsfmri,  demo_CANBIND.SHAPS_Tot, demo_CANBIND.AGE, demo_CANBIND.SEX,...
    madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1, demo_CANBIND.BL_BMI ];%   , demo_CANBIND.lh_hippo, demo_CANBIND.rh_hippo];%    
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
varnames(ix==0)=[];XTest(:, ix==0)=[];

YTest=demo_CANBIND.response==1; 
yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUC_CANBIND] = perfcurve(YTest,yhat, 1);AUC_CANBIND %AUC_test(randfold)=AUC;

%figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')

%% EMBARC stage 2 test
load('clin_embarc_data.mat');
XTest=[clin_stage2.rsfmri, clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17,...
    clin_stage2.employ<=2, clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %
clin_stage2.w16_score_17(isnan(clin_stage2.w16_score_17))=clin_stage2.w10_score_17(isnan(clin_stage2.w16_score_17));
YTest=(clin_stage2.w16_score_17./clin_stage2.w0_score_17<0.5); %Ytest=clin_ct.remission_step2;
XTest(:, ix==0)=[];

%YTest=clin_stage2.remission_step2;
XTest(isnan(clin_stage2.w16_score_17),:)=[]; YTest(isnan(clin_stage2.w16_score_17))=[]; 
XTest(isnan(YTest),:)=[]; YTest(isnan(YTest))=[]; 
%        YTest(isnan(XTest(:,1)) )=[]; XTest(isnan(XTest(:,1) ),:)=[]; XTest(:,roi)=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep2EMB] = perfcurve(YTest,yhat, 1);AUCstep2EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color','r'); plot([0,1],[0,1],'k'); 

%% EMBARC placebo test
XTest=[clin_pla.rsfmri, clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,clin_pla.employ<=2, clin_pla.bmi];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
YTest=strcmp(clin_pla.w8_responder,'Yes');
%Ylogist=clin_ct.remission_step1;
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];
XTest(:, ix==0)=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCstep1PLA_EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 

cd C:\Users\peter\Documents\CANBIND\reports
output=[coef', AUC, AUC_CANBIND, AUCstep2EMB, AUCstep1PLA_EMB, AUC_clinCANBIND, AUCclin2EMB, AUCclin1PLA_EMB];