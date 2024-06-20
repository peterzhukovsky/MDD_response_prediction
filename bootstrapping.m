%% bootstrapping CANBIND > EMBARC baseline response
bootfun = @(x, y)elasticnet(x, y);
%[ci, bootstat] = bootci(3, {bootfun, Xlogist, Ylogist}, 'Options', statset('UseParallel', true))

%[ci, bootstat] = bootci(2, bootfun, Xlogist, Ylogist)

Xlogist=[ demo_CANBIND.SHAPS, demo_CANBIND.AGE, demo_CANBIND.SEX, madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1 , demo_CANBIND.BL_BMI  ]; varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 
Xlogist=[demo_CANBIND.rsfmri,  demo_CANBIND.SHAPS, demo_CANBIND.AGE, demo_CANBIND.SEX, madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1 , demo_CANBIND.BL_BMI  ]; varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 

Ylogist=demo_CANBIND.response==1; 

[bootstat_coef, bootids] = bootstrp(500, bootfun, Xlogist, Ylogist);
figure; %sum(bootstat_coef==0)./10<0.11

mean(bootstat_coef)
std(bootstat_coef)

%% visuals- combine three bootstraps into 1
load("bootstrap\coefficients_50boot_366predictors_alpha0.3_dACC.mat");
tmp=bootstat_coef;
load("bootstrap\coefficients_50boot_366predictors_alpha0.4_dACC.mat");
tmp=[bootstat_coef;tmp];
load("bootstrap\coefficients_50boot_366predictors_alpha0.4_dACC_noCOMBAT.mat");
bootstat_coef=[bootstat_coef;tmp];
mean(bootstat_coef)
sortrows(bootstat_coef(:,370)-bootstat_coef(:,371))

bootstat_coef(:,61)=bootstat_coef(:,55)-bootstat_coef(:,58);
bootstat_coef(:,62)=bootstat_coef(:,56)-bootstat_coef(:,59);
bootstat_coef(:,63)=bootstat_coef(:,57)-bootstat_coef(:,60);
sum(bootstat_coef<0.5)/length(bootstat_coef(:,1))
sum(bootstat_coef<0)/length(bootstat_coef(:,1))

%% clear the 0.5 models
bootstat_coef(bootstat_coef(:,8)==0.5,:)=[];
bootstat_coef(bootstat_coef(:,56)==0.5,:)=[];
sortrows(bootstat_coef(:,49)-bootstat_coef(:,52))

%% cut the extra variables in the dacc model - 50 predictor model reverse from EMBARC TO CAN
load('C:\Users\peter\Documents\CANBIND\reports\test_coef.mat','coef'); 
coef(1)=[];
ix=coef~=0;
Xlogist=[clin_ct.rsfmri, clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, clin_ct.employ<=2, clin_ct.bmi];%       ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  
Ylogist=strcmp(clin_ct.w8_responder,'Yes'); 
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
Xlogist(strcmp(clin_ct.w8_responder,''),:)=[];Ylogist(strcmp(clin_ct.w8_responder,''))=[];
Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 
%varnames(ix==0)=[];Xlogist(:, ix==0)=[];

bootfun = @(x, y)elasticnet_embarc(x, y);
[bootstat_coef, bootids] = bootstrp(1000, bootfun, Xlogist, Ylogist);

%sum(bootstat_coef==0)./10<0.11

mean(bootstat_coef)
std(bootstat_coef)

%% bootstrapping mid treatment models from CANBIND to EMBARC

Xlogist=[demo_CANBIND.rsfmri,  demo_CANBIND.SHAPS, demo_CANBIND.AGE, demo_CANBIND.SEX, madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', ...
    (madrs2hdrs17(demo_CANBIND.MADRS_2wks)-madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED))'./madrs2hdrs17(demo_CANBIND.MADRS_TOT_PRO_RATED)', demo_CANBIND.EMPLOY_STATUS==1 , demo_CANBIND.BL_BMI  ];%  , demo_CANBIND.lh_hippo, demo_CANBIND.rh_hippo];%    ];%  ];%, demo_CANBIND.GAD7_Tot];% 
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'}, {'MADRSchange'}, {'employ'}   ,{'bmi'}  ];%   ,{'lHippo'},{'rHippo'}  ];%            % ,{'gad'}     ,{'linfTemp'},{'rinfTemp'}]; 
Ylogist=demo_CANBIND.response==1; 

bootfun = @(x, y)elasticnet_midtreat(x, y);
[bootstat_coef, bootids] = bootstrp(1000, bootfun, Xlogist, Ylogist);

%sum(bootstat_coef==0)./10<0.11

mean(bootstat_coef)
std(bootstat_coef)

bootstat_coef(bootstat_coef(:,8)==0.5,:)=[];


%% bootstrapping mid treatment models from EMBARC to CANBIND 

Xlogist=[clin_ct.rsfmri, clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, (clin_ct.w2_score_17-clin_ct.w0_score_17)./clin_ct.w0_score_17, clin_ct.employ<=2, clin_ct.bmi];%       ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  
Ylogist=strcmp(clin_ct.w8_responder,'Yes'); 
varnames=[cort_rois.Var2',{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'MADRSchange'}, {'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
Xlogist(strcmp(clin_ct.w8_responder,''),:)=[];Ylogist(strcmp(clin_ct.w8_responder,''))=[];
Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 

bootfun = @(x, y)elasticnet_midtreat_embarc(x, y);
[bootstat_coef, bootids] = bootstrp(1000, bootfun, Xlogist, Ylogist);

%sum(bootstat_coef==0)./10<0.11

mean(bootstat_coef)
std(bootstat_coef)

bootstat_coef(bootstat_coef(:,10)==0.5,:)=[];



























%% Clinical 5 predictor model reverse from EMBARC TO CAN

bootfun = @(x, y)elasticnet_embarc_baseclin(x, y);

Xlogist=[ clin_ct.shaps, clin_ct.age, clin_ct.sex, clin_ct.w0_score_17, clin_ct.employ<=2, clin_ct.bmi];%       ,clin_ct.lh_hippo,clin_ct.rh_hippo]; %  
Ylogist=strcmp(clin_ct.w8_responder,'Yes'); 
varnames=[{'shaps'},{'age'}, {'sex'}, {'blmadrs'},{'employ'}, {'bmi'}];%   {'lHippo'},{'rHippo'}  ];%    ,{'linfTemp'},{'rinfTemp'}]; 
Xlogist(strcmp(clin_ct.w8_responder,''),:)=[];Ylogist(strcmp(clin_ct.w8_responder,''))=[];
Xlogist(isnan(Ylogist),:)=[]; Ylogist(isnan(Ylogist))=[]; 

[bootstat_coef, bootids] = bootstrp(500, bootfun, Xlogist, Ylogist);
figure; %sum(bootstat_coef==0)./10<0.11
