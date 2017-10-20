%STAT550_HW2_5_18.m

addpath(genpath('/Users/ADV/Documents/'))
clear all
%col1: history; col2: verbal; col3: science
X = xlsread('Table5_2.xlsx');
xbar = mean(X)';
S = cov(X);
S_inv=inv(S);
[n,p] = size(X);
mu0 = [500;50;30];

T2 = n*(xbar-mu0)'*S_inv*(xbar-mu0);
c2 = (n-1)*p/(n-p)*2.7;

[vec,val] = eig(S);

history = X(:,1);
verbal = X(:,2);
science = X(:,3);

% figure(1); clf;
% subplot(131)
% qqplot(history,verbal);
% xlabel('History quantiles')
% ylabel('Verbal quantiles')
% title('History vs. Verbal QQ plot')
% 
% figure(2); clf;
% subplot(132)
% qqplot(history,science);
% xlabel('History quantiles')
% ylabel('Science quantiles')
% title('History vs. Science QQ plot')
% 
% figure(3); clf;
% subplot(133)
% qqplot(verbal,science);
% xlabel('Verbal quantiles')
% ylabel('Science quantiles')
% title('Verbal vs. Science QQ plot')

figure(2); clf;
subplot(131)
scatter(history,verbal)
xlabel('History scores')
ylabel('Verbal scores')
title('History vs. Verbal scatter plot')
lsline

%figure(5); clf;
subplot(132)
scatter(history,science)
xlabel('History scores')
ylabel('Science scores')
title('History vs. Science scatter plot')
lsline

%figure(6); clf;
subplot(133)
scatter(verbal,science)
xlabel('Verbal scores')
ylabel('Science scores')
title('Verbal vs. Science scatter plot')
lsline


figure(3); clf;
subplot(131)
qqplot(history);
ylabel('History quantiles')
title('History QQ plot')

subplot(132)
qqplot(verbal);
ylabel('Verbal quantiles')
title('Verbal QQ plot')

subplot(133)
qqplot(science);
ylabel('Science quantiles')
title('Science QQ plot')