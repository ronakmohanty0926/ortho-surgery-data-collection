 clc;
clear;
A = readmatrix("YB3_Users_Sig_Arc.xlsx", 'Sheet',2);
idxTrial = 2;
X = A(:,idxTrial);
% nGauss = 2;
% options = statset('MaxIter',500);
% gmm = fitgmdist(X, nGauss, 'Options',options);
% intv = 0:0.025:1;
% intv = intv';
% pdfgmm = pdf(gmm,intv);
figure
% hold on
% plot(intv,pdfgmm,'LineWidth',2);
h = histogram(X,'Normalization','cdf');

% hold off