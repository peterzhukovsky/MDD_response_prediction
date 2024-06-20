%% 4. RUN THE TRAINING MODEL IN EMBARC 1 SSRI %% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %%
%% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %% %% %% %%%% %%%% %%%% %%%% %%%% %%%% %%
load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef'); load('C:\Users\peter\Documents\CANBIND\reports\summary_tests\dACC_testing.mat' ,'coef')
coef(1)=[];
ix=coef~=0;
Xlogist=[clin_ct.rsfmri, clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, clin_ct.employ<=2, clin_ct.bmi];%       ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  
Ylogist=strcmp(clin_ct.w8_responder,'Yes'); 
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
Xlogist(strcmp(clin_ct.w8_responder,''),:)=[];Ylogist(strcmp(clin_ct.w8_responder,''))=[];
Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 
varnames(ix==0)=[];Xlogist(:, ix==0)=[];
AUC=0.5;
while AUC==0.5
[B,FitInfo] = lassoglm(Xlogist,Ylogist,'binomial','CV',10, 'PredictorNames', varnames      ,'alpha', 0.001); %    0.001   
%lassoPlot(B,FitInfo,'PlotType','CV'); legend('show','Location','best') % show legend
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
MinModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinDeviance)~=0)
idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];
yhat = glmval(coef,Xlogist,'logit');[X,Y,T,AUC] = perfcurve(Ylogist,yhat, 1);AUC
end
%% CANBIND test
XTest=[demo_CANBIND.rsfmri,  demo_CANBIND.SHAPS_Tot, demo_CANBIND.AGE, demo_CANBIND.SEX,...
    madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1, demo_CANBIND.BL_BMI ];%   , demo_CANBIND.lh_hippo, demo_CANBIND.rh_hippo];%    
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
varnames(ix==0)=[];XTest(:, ix==0)=[];

YTest=demo_CANBIND.response==1; 
yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUC_CANBIND] = perfcurve(YTest,yhat, 1);AUC_CANBIND %AUC_test(randfold)=AUC;

figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')
tmp=(mean([1-X, Y]')); J=Y-X; if sum(1-X >0.55 & Y>0.55)>0; J(1-X<0.55 | Y<0.55)=0;  end; J(1-X<0.5 | Y<0.5)=0;  ixx=J==max(J); bAcc(1)=tmp(ixx); Sensitivity(1)=1-X(ixx);Specificity(1)=Y(ixx);

%% EMBARC stage 2 test
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
figure(2);hold on;plot(X,Y,'Color','r'); plot([0,1],[0,1],'k'); 

tmp=(mean([1-X, Y]')); J=Y-X; if sum(1-X >0.55 & Y>0.55)>0; J(1-X<0.55 | Y<0.55)=0;  end;   J(1-X<0.5 | Y<0.5)=0; ixx=J==max(J); bAcc(2)=tmp(ixx);Sensitivity(2)=1-X(ixx);Specificity(2)=Y(ixx);

%% EMBARC placebo test
XTest=[clin_pla.rsfmri, clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,  clin_pla.employ<=2, clin_pla.bmi];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
YTest=strcmp(clin_pla.w8_responder,'Yes');
%Ylogist=clin_ct.remission_step1;
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];
XTest(:, ix==0)=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCstep1PLA_EMB %AUC_test(randfold)=AUC;
figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 

tmp=(mean([1-X, Y]')); J=Y-X; if sum(1-X >0.55 & Y>0.55)>0; J(1-X<0.55 | Y<0.55)=0;  end;  J(1-X<0.5 | Y<0.5)=0; ixx=J==max(J); bAcc(3)=tmp(ixx); Sensitivity(3)=1-X(ixx);Specificity(3)=Y(ixx);

cd C:\Users\peter\Documents\CANBIND\reports

%tmp=(mean([1-X, Y]')); J=Y-X; ixx=J==max(J); bAcc(3)=tmp(ixx); Sensitivity(3)=1-X(ixx);Specificity(3)=Y(ixx);


