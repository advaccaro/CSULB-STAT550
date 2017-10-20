%STAT550_HW2_5_20.m
addpath(genpath('/Users/ADV/Documents/'))
clear all

X = xlsread('Table5_12.xlsx');
[n,p] = size(X);
xbar = mean(X)';
S = cov(X);
S_inv = inv(S);
mu0 = [190;275];

[vec,val] = eig(S);
T2 = n*(xbar-mu0)'*S_inv*(xbar-mu0);
c2 = (n-1)*p/(n-p)*3.2;

tail=X(:,1);
wing=X(:,2);

figure(1); clf;
subplot(121)
qqplot(tail);
ylabel('Tail length quantiles')
title('Tail length QQ plot')

subplot(122)
qqplot(wing)
ylabel('Wing length quantiles')
title('Wing length QQ plot')

figure(2); clf;
scatter(tail,wing);
xlabel('Tail length')
ylabel('Wing length')
title('Tail vs. Wing length scatter plot')
lsline
