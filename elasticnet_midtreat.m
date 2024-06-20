function [output] = elasticnet_midtreat(Xlogist,Ylogist)


%% TEST CLINICAL MODEL IN THE SAME BOOTSTRAP 

varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'}, {'MADRSchange'}, {'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 


[B,FitInfo] = lassoglm(Xlogist(:,361:367),Ylogist,'binomial','CV',10, 'PredictorNames', varnames  ,'alpha', 0.01); %0.01 with the base only
%lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
%yhat = glmval(coef,Xlogist,'logit');[X,Y,T,AUC] = perfcurve(Ylogist,yhat, 1);AUC
%Xlogist=Xlogist(:,B(:,idxLambda1SE)~=0);varnames=varnames(B(:,idxLambda1SE)~=0); %pasimonious model test > set alpha to very low 0.001
%save('test_coef.mat','coef')

%% test in embarc
load('C:\Users\peter\Documents\CANBIND\reports\clin_embarc_data.mat');


Xlogist=[clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, (clin_ct.w2_score_17-clin_ct.w0_score_17)./clin_ct.w0_score_17,  clin_ct.employ<=2       , clin_ct.bmi   ];%    ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  ];% 
Ylogist=strcmp(clin_ct.w8_responder,'Yes'); 
%Ylogist=clin_ct.remission_step1;
Xlogist(strcmp(clin_ct.w8_responder,''),:)=[];Ylogist(strcmp(clin_ct.w8_responder,''))=[];
Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 
%       Ylogist(isnan(Xlogist(:,1)) )=[]; Xlogist(isnan(Xlogist(:,1) ),:)=[]; Xlogist(:,roi)=[];varnames(roi)=[];
%load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef'); 
yhat = glmval(coef,Xlogist,'logit'); [X,Y,T,AUCclin1EMB] = perfcurve(Ylogist,yhat, 1);AUCclin1EMB %AUC_test(randfold)=AUC;

%figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')

%%% embarc step 2
XTest=[clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17, (clin_stage2.w10_score_17-clin_stage2.w0_score_17)./clin_stage2.w0_score_17, clin_stage2.employ<=2    , clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %  ];% 
clin_stage2.w16_score_17(isnan(clin_stage2.w16_score_17))=clin_stage2.w12_score_17(isnan(clin_stage2.w16_score_17));
YTest=(clin_stage2.w16_score_17./clin_stage2.w0_score_17<0.5);
%YTest=clin_stage2.remission_step2;
XTest(isnan(clin_stage2.w16_score_17),:)=[]; YTest(isnan(clin_stage2.w16_score_17))=[]; 
XTest(isnan(YTest),:)=[]; YTest(isnan(YTest))=[]; 
%        YTest(isnan(XTest(:,1)) )=[]; XTest(isnan(XTest(:,1) ),:)=[]; XTest(:,roi)=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCclin2EMB] = perfcurve(YTest,yhat, 1);AUCclin2EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color','r'); plot([0,1],[0,1],'k'); 

cd C:\Users\peter\Documents\CANBIND\reports

%%% embarc placebo step 1
XTest=[clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,  (clin_pla.w2_score_17-clin_pla.w0_score_17)./clin_pla.w0_score_17, clin_pla.employ<=2   , clin_pla.bmi   ];%  ,clin_pla.lh_hippo,clin_pla.rh_hippo  ];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
YTest=strcmp(clin_pla.w8_responder,'Yes');
%Ylogist=clin_ct.remission_step1;
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCclin1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCclin1PLA_EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 

cd C:\Users\peter\Documents\CANBIND\reports

output=[coef', AUCclin1EMB, AUCclin2EMB, AUCclin1PLA_EMB];
%output=[coef', AUC, AUC_CANBIND, AUCstep2EMB, AUCstep1PLA_EMB, AUC_clinCANBIND, AUCclin2EMB, AUCclin1PLA_EMB];




