%% scripts for pet image processing -  R noR
clear
cd C:\Users\peter\Documents\GABA\PET\'5HT and DA'\
load('C:\Users\peter\Documents\CANBIND\reports\AHBA\dACC_outputs.mat');
c3vs1_tstat=[tstatCAN', tstatEMB1', tstatEMB2'];

ht1a_left=readtable('lh.5HT1a_way_hc36_savli.sum.txt');ht1a_left.Mean(120)=NaN;
ht1a_right=readtable('rh.5HT1a_way_hc36_savli.sum.txt');ht1a_right.Mean(120)=NaN;

ht1b_left=readtable('lh.5HT1b_p943_hc22_savli.sum.txt');ht1a_left.Mean(120)=NaN;
ht1b_right=readtable('rh.5HT1b_p943_hc22_savli.sum.txt');ht1a_right.Mean(120)=NaN;


ht2a_left=readtable('lh.5HT2a_alt_hc19_savli.sum.txt'); ht2a_left.Mean(120)=NaN;
ht2a_right=readtable('rh.5HT2a_alt_hc19_savli.sum.txt'); ht2a_right.Mean(120)=NaN;
ht2a_left2=readtable('lh.5HT2a_mdl_hc3_talbot.sum.txt'); ht2a_left2.Mean(120)=NaN;
ht2a_right2=readtable('rh.5HT2a_mdl_hc3_talbot.sum.txt'); ht2a_right2.Mean(120)=NaN;


htt_left=readtable('lh.5HTT_dasb_hc30_savli.sum.txt');%htt_left.Mean(120)=NaN;
htt_right=readtable('rh.5HTT_dasb_hc30_savli.sum.txt');%htt_right.Mean(120)=NaN;
htt_left2=readtable('lh.5HTT_madam_hc10_fazio.sum.txt'); %htt_left2.Mean(120)=NaN;
htt_right2=readtable('rh.5HTT_madam_hc10_fazio.sum.txt');% htt_right2.Mean(120)=NaN;

cd C:\Users\peter\Documents\CANBIND\reports
ht2a_left=mean([(ht2a_left.Mean), (ht2a_left2.Mean)]')';
ht2a_right=mean([(ht2a_right.Mean), (ht2a_right2.Mean)]')';
htt_left=mean([(htt_left.Mean), (htt_left2.Mean)]')';
htt_right=mean([(htt_right.Mean), (htt_right2.Mean)]')';

%
[r,p]=corr(ht1b_left.Mean, c3vs1_tstat(1:180,:), 'rows', 'pairwise')
imagesc(r([1,3,5,2,4])'); colormap bone; colorbar
[r,p]=corr(ht1b_right.Mean, c3vs1_tstat(181:360,:), 'rows', 'pairwise')
imagesc(r([6,8,10,7,9])'); colormap bone; colorbar

figure(1); scatter(ht1b_left.Mean, c3vs1_tstat(1:180,1), 'filled', 'k'); lsline; %hold on; ix=bhfdr(c3vs1_pval(:,10))<0.05;
figure(2); scatter(ht1b_right.Mean, c3vs1_tstat(181:360,1), 'filled', 'k'); lsline; %hold on; ix=bhfdr(c3vs1_pval(:,10))<0.05;


%'L_accumbens', 'L_amygdala', 'L_caudate', 'L_hippocampus', 'L_putamen
[r,p]=corr(ht1a_left.Mean, c3vs1_tstat(1:180,:), 'rows', 'pairwise')
imagesc(r([1,3,5,2,4])'); colormap bone; colorbar
[r,p]=corr(ht1a_right.Mean, c3vs1_tstat(181:360,:), 'rows', 'pairwise')
imagesc(r([6,8,10,7,9])'); colormap bone; colorbar

figure(1); scatter(ht1a_left.Mean, c3vs1_tstat(1:180,1), 'filled', 'k'); lsline; %hold on; ix=bhfdr(c3vs1_pval(:,10))<0.05;
figure(2); scatter(ht1a_right.Mean, c3vs1_tstat(181:360,1), 'filled', 'k'); lsline; %hold on; ix=bhfdr(c3vs1_pval(:,10))<0.05;


[r,p]=corr(ht2a_left,  c3vs1_tstat(1:180,:), 'rows','pairwise')
figure; imagesc(r([1,3,5,2,4])'); colormap bone; colorbar
[r,p]=corr(ht2a_right, c3vs1_tstat(181:360,:), 'rows','pairwise')
imagesc(r([6,8,10,7,9])'); colormap bone; colorbar


figure(1); hold off; scatter(gaba_right.Mean, c3vs1_tstat(:,10), 'filled', 'k'); lsline; hold on; ix=bhfdr(c3vs1_pval(:,10))<0.05;
scatter(gaba_right.Mean(ix), c3vs1_tstat(ix,10), 40,'r');

figure(2); hold off; scatter(gaba_left.Mean, c3vs1_tstat(:,5), 'filled', 'k'); lsline; hold on; ix=bhfdr(c3vs1_pval(:,5))<0.05;
scatter(gaba_left.Mean(ix), c3vs1_tstat(ix,5), 40,'r');


[r,p]=corr(htt_left,  c3vs1_tstat(1:180,:), 'rows','pairwise')
figure; imagesc(r([1,3,5,2,4])'); colormap bone; colorbar
[r,p]=corr(htt_right, c3vs1_tstat(181:360,:), 'rows','pairwise')
imagesc(r([6,8,10,7,9])'); colormap bone; colorbar


figure(1); hold off; scatter(kappa_right.Mean, c3vs1_tstat(:,10), 'filled', 'k'); lsline; hold on; ix=bhfdr(c3vs1_pval(:,10))<0.05;
scatter(kappa_right.Mean(ix), c3vs1_tstat(ix,10), 40,'r');

figure(2); hold off; scatter(kappa_left.Mean, c3vs1_tstat(:,5), 'filled', 'k'); lsline; hold on; ix=bhfdr(c3vs1_pval(:,5))<0.05;
scatter(kappa_left.Mean(ix), c3vs1_tstat(ix,5), 40,'r');


