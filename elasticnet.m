function [output] = elasticnet(Xlogist,Ylogist)


%% TEST CLINICAL MODEL IN THE SAME BOOTSTRAP 

varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   

[B,FitInfo] = lassoglm(Xlogist(:,361:366),Ylogist,'binomial','CV',10, 'PredictorNames', varnames  ,'alpha', 0.01); %0.01 with the base only
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
Xtest=[clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17,  clin_ct.employ<=2       , clin_ct.bmi   ];%    ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  ];% 
%Ylogist=clin_ct.remission_step1;
Ytest=strcmp(clin_ct.w8_responder,'Yes'); 
Xtest(strcmp(clin_ct.w8_responder,''),:)=[];Ytest(strcmp(clin_ct.w8_responder,''))=[];
Xtest(isnan(Ytest),:)=[]; Ytest(isnan(Ytest))=[]; 
%       Ylogist(isnan(Xlogist(:,1)) )=[]; Xlogist(isnan(Xlogist(:,1) ),:)=[]; Xlogist(:,roi)=[];varnames(roi)=[];
%load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef'); 
yhat = glmval(coef,Xtest,'logit'); [X,Y,T,AUCclin1EMB] = perfcurve(Ytest,yhat, 1);AUCclin1EMB %AUC_test(randfold)=AUC;

%figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')

%%% embarc step 2
XTest=[clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17, clin_stage2.employ<=2    , clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %  ];% 
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
XTest=[clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17, clin_pla.employ<=2   , clin_pla.bmi   ];% ,clin_pla.lh_hippo,clin_pla.rh_hippo  ];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
YTest=strcmp(clin_pla.w8_responder,'Yes');
%YTest=clin_ct.remission_step1;
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCclin1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCclin1PLA_EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 





























%%
cort_rois=readtable('C:\Users\peter\Documents\GABA\HCPMMP1_on_MNI152_ICBM2009a_nlin.txt');
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   

%%%%% unhash these 4 lines for the compact model
%load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef');
%coef(1)=[];
%ix=coef~=0;
%varnames(ix==0)=[];Xlogist(:, ix==0)=[];

[B,FitInfo] = lassoglm(Xlogist,Ylogist,'binomial','CV',10, 'PredictorNames', varnames,'alpha',0.3); %0.3 0.01 with the base only
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = glmval(coef,Xlogist,'logit');[X,Y,T,AUC] = perfcurve(Ylogist,yhat, 1);AUC
save('C:\Users\peter\Documents\CANBIND\reports\bootstrap\test_coef.mat','coef');

load('C:\Users\peter\Documents\CANBIND\reports\clin_embarc_data.mat');

%% 

XTest=[clin_ct.rsfmri, clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, clin_ct.employ<=2       , clin_ct.bmi   ];%    ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  ];% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    XTest=[clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, clin_ct.employ<=2       , clin_ct.bmi   ];%    ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  ];% 
Ytest=strcmp(clin_ct.w8_responder,'Yes'); 
XTest(strcmp(clin_ct.w8_responder,''),:)=[];Ytest(strcmp(clin_ct.w8_responder,''))=[];
XTest(isnan(Ytest),:)=[]; Ytest(isnan(Ytest))=[]; 
%XTest(:, ix==0)=[]; %%%%%%%%%%%%%%%%%%%%%% UNHASH line 

%       Ylogist(isnan(Xlogist(:,1)) )=[]; Xlogist(isnan(Xlogist(:,1) ),:)=[]; Xlogist(:,roi)=[];varnames(roi)=[];
load('C:\Users\peter\Documents\CANBIND\reports\bootstrap\test_coef.mat','coef');
yhat = glmval(coef,XTest,'logit'); [X,Y,T,AUCstep1EMB] = perfcurve(Ytest,yhat, 1);AUCstep1EMB %AUC_test(randfold)=AUC;

%figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')

%%
XTest=[clin_stage2.rsfmri, clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17,clin_stage2.employ<=2    , clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %  ];% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    XTest=[clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17,clin_stage2.employ<=2    , clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %  ];% 
clin_stage2.w16_score_17(isnan(clin_stage2.w16_score_17))=clin_stage2.w12_score_17(isnan(clin_stage2.w16_score_17));
YTest=(clin_stage2.w16_score_17./clin_stage2.w0_score_17<0.5);
%YTest=clin_stage2.remission_step2;
XTest(isnan(clin_stage2.w16_score_17),:)=[]; YTest(isnan(clin_stage2.w16_score_17))=[]; 
XTest(isnan(YTest),:)=[]; YTest(isnan(YTest))=[]; 
%        YTest(isnan(XTest(:,1)) )=[]; XTest(isnan(XTest(:,1) ),:)=[]; XTest(:,roi)=[];
% XTest(:, ix==0)=[]; %%%%%%%%%%%%%%%%%%%%%% UNHASH line 
yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep2EMB] = perfcurve(YTest,yhat, 1);AUCstep2EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color','r'); plot([0,1],[0,1],'k'); 

cd C:\Users\peter\Documents\CANBIND\reports

%%
XTest=[clin_pla.rsfmri, clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,clin_pla.employ<=2   , clin_pla.bmi   ];%  ,clin_pla.lh_hippo,clin_pla.rh_hippo  ];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   XTest=[clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,clin_pla.employ<=2   , clin_pla.bmi   ];%  
YTest=strcmp(clin_pla.w8_responder,'Yes');
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];
% XTest(:, ix==0)=[]; %%%%%%%%%%%%%%%%%%%%%% UNHASH line 
yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCstep1PLA_EMB %AUC_test(randfold)=AUC;
%figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 
output=[coef', AUC, AUCstep1EMB, AUCstep2EMB, AUCstep1PLA_EMB,  AUCclin1EMB, AUCclin2EMB, AUCclin1PLA_EMB];

