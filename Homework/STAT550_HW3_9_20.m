%STAT550_HW3_9_20.m
%% Import data from text file.
clear all
%% Initialize variables.
filename = '/Users/ADV/Documents/MATLAB/T1-5.dat';
delimiter = ' ';

%% Format for each line of text:
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
VarName2 = dataArray{:, 1};
VarName3 = dataArray{:, 2};
VarName4 = dataArray{:, 3};
VarName5 = dataArray{:, 4};
VarName6 = dataArray{:, 5};
VarName7 = dataArray{:, 6};
VarName8 = dataArray{:, 7};


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Rename
X1 = VarName2; X2 = VarName3; X3 = VarName4; X4 = VarName5;
X5 = VarName6; X6 = VarName7; X7 = VarName8;

%% Combine
X = [X1 X2 X5 X6];

%% Covariance
S = cov(X);
R = corr(X);

%% Part a - PCA soln for m =1 and m =2
[V,D] = eig(S);
for i = 1:length(D)
    eigs(i) = D(i,i);
end
eigs_sort = fliplr(sort(eigs));
for i = 1:length(D)
    ind = find(eigs_sort(i) == eigs);
    V_sort(:,i) = V(:,ind);
end

lambda = eigs_sort;
ev = V_sort;


%% m = 1:
A1.L = sqrt(lambda(1))*V(:,1);
A1.Psi = zeros(length(ev));
for i = 1:4
    A1.Psi(i,i) = 1 - A1.L(i).^2;
end
%% m = 2:
A2.L = [sqrt(lambda(1))*V(:,1) sqrt(lambda(2))*V(:,2)];
A2.Psi = zeros(4);
for i = 1:4
    A2.Psi(i,i) = 1 - A2.L(i,1).^2 - A2.L(i,2).^2;
end

%% Part b

[B1.Lambda,B1.Psi,B1.T] = factoran(S,1,'xtype','cov');
B1Psi = zeros(4);
for i = 1:4
    B1Psi(i,i) = B1.Psi(i);
end
B1.Psi = B1Psi;
[B2.Lambda,B2.Psi,B2.T] = factoran(S,2,'xtype','cov');

%% Part c

