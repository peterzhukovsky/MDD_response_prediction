%% test can in embarc
%% importing data and merging sheets
%% embarc step 1

Xlogist=[clin_ct.rsfmri, clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, clin_ct.employ<=2       , clin_ct.bmi   ];%    ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  ];% 
%varnames=[cort_rois.Var2',{'shaps'}, {'age'}, {'sex'},{'blmadrs'},{'employ'}];%  
Ylogist=strcmp(clin_ct.w8_responder,'Yes'); 
%Ylogist=clin_ct.remission_step1;
Xlogist(strcmp(clin_ct.w8_responder,''),:)=[];Ylogist(strcmp(clin_ct.w8_responder,''))=[];
Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 
%       Ylogist(isnan(Xlogist(:,1)) )=[]; Xlogist(isnan(Xlogist(:,1) ),:)=[]; Xlogist(:,roi)=[];varnames(roi)=[];
load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef'); 
yhat = glmval(coef,Xlogist,'logit'); [X,Y,T,AUCstep1EMB] = perfcurve(Ylogist,yhat, 1);AUCstep1EMB %AUC_test(randfold)=AUC;

figure(2);hold on;plot(X,Y,'Color',[76/255 0 153/255]); plot([0,1],[0,1],'k')
tmp=(mean([1-X, Y]')); J=Y-X; if sum(1-X >0.55 & Y>0.55)>0; J(1-X<0.55 | Y<0.55)=0;  end; J(1-X<0.5 | Y<0.5)=0;  ix=J==max(J); bAcc(1)=tmp(ix); Sensitivity(1)=1-X(ix);Specificity(1)=Y(ix);
%% embarc step 2
XTest=[clin_stage2.rsfmri, clin_stage2.shaps, clin_stage2.age, clin_stage2.sex, clin_stage2.w0_score_17,clin_stage2.employ<=2    , clin_stage2.bmi  ];%    ,clin_stage2.lh_hippo,clin_stage2.rh_hippo];  %  ];% 
clin_stage2.w16_score_17(isnan(clin_stage2.w16_score_17))=clin_stage2.w12_score_17(isnan(clin_stage2.w16_score_17));
YTest=(clin_stage2.w16_score_17./clin_stage2.w0_score_17<0.5);
%YTest=clin_stage2.remission_step2;
XTest(isnan(clin_stage2.w16_score_17),:)=[]; YTest(isnan(clin_stage2.w16_score_17))=[]; 
XTest(isnan(YTest),:)=[]; YTest(isnan(YTest))=[]; 
%        YTest(isnan(XTest(:,1)) )=[]; XTest(isnan(XTest(:,1) ),:)=[]; XTest(:,roi)=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep2EMB] = perfcurve(YTest,yhat, 1);AUCstep2EMB %AUC_test(randfold)=AUC;
figure(2);hold on;plot(X,Y,'Color','r'); plot([0,1],[0,1],'k'); 

cd C:\Users\peter\Documents\CANBIND\reports
tmp=(mean([1-X, Y]')); J=Y-X; if sum(1-X >0.55 & Y>0.55)>0; J(1-X<0.55 | Y<0.55)=0;  end;   J(1-X<0.5 | Y<0.5)=0; ix=J==max(J); bAcc(2)=tmp(ix);Sensitivity(2)=1-X(ix);Specificity(2)=Y(ix);

%% embarc placebo step 1
XTest=[clin_pla.rsfmri, clin_pla.shaps, clin_pla.age, clin_pla.sex, clin_pla.w0_score_17,clin_pla.employ<=2   , clin_pla.bmi   ];%  ,clin_pla.lh_hippo,clin_pla.rh_hippo  ];%    ,clin_stage2.lh_inferiortemporal_thickness,clin_stage2.rh_inferiortemporal_thickness ];  
YTest=strcmp(clin_pla.w8_responder,'Yes');
%Ylogist=clin_ct.remission_step1;
XTest(strcmp(clin_pla.w8_responder,''),:)=[];YTest(strcmp(clin_pla.w8_responder,''))=[];

yhat = glmval(coef,XTest,'logit');
[X,Y,T,AUCstep1PLA_EMB] = perfcurve(YTest,yhat, 1);AUCstep1PLA_EMB %AUC_test(randfold)=AUC;
figure(2);hold on;plot(X,Y,'Color',[0/255 204/255 0/255]); plot([0,1],[0,1],'k'); 
tmp=(mean([1-X, Y]')); J=Y-X; if sum(1-X >0.55 & Y>0.55)>0; J(1-X<0.55 | Y<0.55)=0;  end;  J(1-X<0.5 | Y<0.5)=0; ix=J==max(J); bAcc(3)=tmp(ix); Sensitivity(3)=1-X(ix);Specificity(3)=Y(ix);

