%STAT550 - Final Project
%Adam Vaccaro

clear all;
addpath(genpath('/Users/ADV/Documents/MATLAB/'))
addpath(genpath('/Users/ADV/Google Drive/CSULB/Spring 2017/STAT 550/FinalProject/data'))
%% Initialize variables.
filename = '/Users/ADV/Google Drive/CSULB/Spring 2017/STAT 550/FinalProject/data/AQI_SB_01_03.csv';
delimiter = ',';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{1} = datetime(dataArray{1}, 'Format', 'MM/dd/yyyy', 'InputFormat', 'MM/dd/yyyy');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
        dates{1} = datetime(dataArray{1}, 'Format', 'MM/dd/yyyy', 'InputFormat', 'MM/dd/yyyy');
    catch
        dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
    end
end

anyBlankDates = cellfun(@isempty, dataArray{1});
anyInvalidDates = isnan(dates{1}.Hour) - anyBlankDates;
dates = dates(:,1);

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
X.date1 = dates{:, 1};
X.oz1h = cell2mat(rawNumericColumns(:, 1));
X.day1 = cell2mat(rawNumericColumns(:, 2));
X.tempmax = cell2mat(rawNumericColumns(:, 3));
X.soilmax = cell2mat(rawNumericColumns(:, 4));
X.solar = cell2mat(rawNumericColumns(:, 5));
X.eto = cell2mat(rawNumericColumns(:, 6));
X.rhmax = cell2mat(rawNumericColumns(:, 7));
X.rhmin = cell2mat(rawNumericColumns(:, 8));
X.prec = cell2mat(rawNumericColumns(:, 9));
X.windu = cell2mat(rawNumericColumns(:, 10));
X.windv = cell2mat(rawNumericColumns(:, 11));
X.month1 = cell2mat(rawNumericColumns(:, 12));
X.weekend = cell2mat(rawNumericColumns(:, 13));
X.AQI = cell2mat(rawNumericColumns(:, 14));
X.AQIg = cell2mat(rawNumericColumns(:,15));

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% date1=datenum(date1);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me dates blankDates anyBlankDates invalidDates anyInvalidDates rawNumericColumns R;

%% Convert and correct data vector using datenum technology
date_vec = datevec(X.date1);
date_vec(:,1) = date_vec(:,1) + 2000;
X.datevec = date_vec;
X.datenum = datenum(X.datevec);
%% Check for normality
% figure(1); clf;
% subplot(221); hist(X.tempmax); title('Temp max');
% subplot(222); hist(X.soilmax); title('Soil max');
% subplot(223); hist(X.solar); title('Solar');
% subplot(224); hist(X.eto); title('Eto');
% figure(2); clf;
% subplot(221); hist(X.rhmax); title('RH max');
% subplot(222); hist(X.windu); title('Wind u');
% subplot(223); hist(X.windv); title('Wind v');
% subplot(224); hist(boxcox(X.rhmax)); title('RH max (T)');
% [T.rhmax,rhlam] = boxcox(X.rhmax);
% [T.windu,ulam] = boxcox(X.windu);
% [T.windu,ulam] = boxcox(X.windu - min(X.windu) +1);
% [T.windv,vlam] = boxcox(X.windv);
% [T.windv,vlam] = boxcox(X.windv - min(X.windv)+1);
% figure(3); clf;
% subplot(221); hist(T.rhmax); title('RH Max (T)');
% subplot(222); hist(T.windu); title('Wind u (T)');
% subplot(223); hist(T.windv); title('Wind v (T)');
% 
% figure(4); clf;
% subplot(221); hist(boxcox(X.tempmax)); title('Temp max (T)');
% subplot(222); hist(boxcox(X.soilmax)); title('Soil max (T)');
% subplot(223); hist(boxcox(X.solar)); title('Solar (T)');
% subplot(224); hist(boxcox(X.eto)); title('ETO (T)');

% figure(1); clf;
% subplot(221); qqplot(X.tempmax); title('Temp max');
% subplot(222); qqplot(boxcox(X.tempmax)); title('Temp max (T)');
% subplot(223); qqplot(X.soilmax); title('Soil max');
% subplot(224); qqplot(boxcox(X.soilmax)); title('Soil max (T)');
% 
% figure(2); clf;
% subplot(221); qqplot(X.solar); title('Solar');
% subplot(222); qqplot(boxcox(X.solar)); title('Solar (T)');
% subplot(223); qqplot(X.eto); title('ETO');
% subplot(224); qqplot(boxcox(X.eto)); title('ETO (T)');
% 
% figure(3); clf;
% subplot(221); qqplot(X.windu); title('Wind U');
% subplot(222); qqplot(boxcox(X.windu-min(X.windu)+1)); title('Wind U (T)');
% subplot(223); qqplot(X.windv); title('Wind V');
% subplot(224); qqplot(boxcox(X.windv-min(X.windu)+1)); title('Wind V (T)');
% 
% figure(4); clf;
% subplot(221); qqplot(X.rhmax); title('RH Max');
% subplot(222); qqplot(boxcox(X.rhmax)); title('RH Max (T)');

figure(1); clf;
subplot(221); qqplot(X.tempmax); title('Temp max');
subplot(222); qqplot(X.soilmax); title('Soil max');
subplot(223); qqplot(X.solar); title('Solar');
subplot(224); qqplot(X.eto); title('ETO');

figure(2); clf;
subplot(221); qqplot(X.rhmax); title('RH max');
subplot(222); qqplot(X.windu); title('Wind U');
subplot(223); qqplot(X.windv); title('Wind V');

%% Create season categorical variable
season = zeros(length(X.month1),1);
in_season = [5:9];
for i = 1:length(season)
    month = X.month1(i);
    if ismember(month,in_season) == 1
        season(i) = 1;
    end
end
X.season = season;    

%% Selects Vars
DATA = [X.tempmax X.soilmax X.solar X.eto X.rhmax X.windu X.windv];
Ni = length(DATA); %initial length
inds = ones(Ni,1); %pre-allocate storage space for index vector

%% Select rows w/o missing values
for i = 1:Ni
    nansum = sum(isnan(DATA(i,:)));
    if nansum > 0
        inds(i) = 0;
    end
end
inds = logical(inds);

%% Create dataset w/o missing values
DATAX = DATA(inds,:);

%% Select dates w/o missing values
DATEX = X.datenum(inds,:);

%% Select other vars w/o missing values
P.weekend = X.weekend(inds,:);
P.month1 = X.month1(inds,:);
P.date1 = datetime(datevec(DATEX));
P.day1 = X.day1(inds,:);
P.AQI = X.AQI(inds,:);
P.AQIg = X.AQIg(inds,:);
P.season =  X.season(inds,:);
%% Calculate Correlation and Covariance (do in SAS to account for NaNs)
S = cov(DATAX);
R = corr(DATAX);
%R = corrcoef(S);

%% Get eigenvalues and eigenvectors (PCs) from correlation matrix
[V,D] = eig(R);
V = -V;

for i = 1:length(V)
    eig_val(i) = D(i,i); %extract eigenvalues from diagonal matrix
end

[eig_val_s,ind_s] = sort(eig_val); %sorts low to high
eig_val_s = fliplr(eig_val_s); %flip to be high to low
ind_s = fliplr(ind_s);


eig_vec_s = zeros(size(V));
for i = 1:length(V)
    ind = ind_s(i); %extract index
    eig_vec_s(:,i) = V(:,ind);
end
    
PC1coef = eig_vec_s(:,1); 
PC2coef = eig_vec_s(:,2); 
PC3coef = eig_vec_s(:,3);

PC1 = DATAX*PC1coef;
PC2 = DATAX*PC2coef;
PC3 = DATAX*PC3coef;





%% Prepare data for export to SAS
DATA_RD = [PC1 PC2 PC3 P.AQIg P.day1 P.month1 P.season];
csvwrite('SB_01_03_PCs_and_categorical.csv',DATA_RD);

%% Discriminant Analysis
PCs = [PC1 PC2 PC3];
PCcat = [PC1 PC2 PC3 P.season P.weekend];
% Get quadratic discriminant model:
MdlQuadratic = fitcdiscr(PCs,P.AQIg,'DiscrimType','quadratic');
%MdlQuadratic2 = fitcdiscr(PCcat,P.AQIg,'DiscrimType','quadratic');
MdlQuadratic3 = fitcdiscr(PCs,P.AQIg,'DiscrimType','quadratic');

% Make predictions:
predictions = predict(MdlQuadratic,PCs);
correct = 0;
for i = 1:length(predictions)
    if predictions(i) == P.AQIg(i)
        correct = correct + 1;
    end
end
p_correct = correct/length(predictions);



%% Import 2004 Data
%% Initialize variables.
filename = '/Users/ADV/Google Drive/CSULB/Spring 2017/STAT 550/FinalProject/data/AQI_SB_04.csv';
delimiter = ',';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{1} = datetime(dataArray{1}, 'Format', 'MM/dd/yyyy', 'InputFormat', 'MM/dd/yyyy');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
        dates{1} = datetime(dataArray{1}, 'Format', 'MM/dd/yyyy', 'InputFormat', 'MM/dd/yyyy');
    catch
        dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
    end
end

anyBlankDates = cellfun(@isempty, dataArray{1});
anyInvalidDates = isnan(dates{1}.Hour) - anyBlankDates;
dates = dates(:,1);

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
Y.date1 = dates{:, 1};
Y.oz1h = cell2mat(rawNumericColumns(:, 1));
Y.day1 = cell2mat(rawNumericColumns(:, 2));
Y.tempmax = cell2mat(rawNumericColumns(:, 3));
Y.soilmax = cell2mat(rawNumericColumns(:, 4));
Y.solar = cell2mat(rawNumericColumns(:, 5));
Y.eto = cell2mat(rawNumericColumns(:, 6));
Y.rhmax = cell2mat(rawNumericColumns(:, 7));
Y.windu = cell2mat(rawNumericColumns(:, 8));
Y.windv = cell2mat(rawNumericColumns(:, 9));
Y.month1 = cell2mat(rawNumericColumns(:, 10));
Y.weekend = cell2mat(rawNumericColumns(:, 11));
Y.AQI = cell2mat(rawNumericColumns(:, 12));
Y.AQIg = cell2mat(rawNumericColumns(:, 13));

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me dates blankDates anyBlankDates invalidDates anyInvalidDates rawNumericColumns R;

%% Convert and correct data vector using datenum technology
date_vec = datevec(Y.date1);
date_vec(:,1) = date_vec(:,1) + 2000;
Y.datevec = date_vec;
Y.datenum = datenum(Y.datevec);
%% Create season categorical variable
season = zeros(length(X.month1),1);
in_season = [5:9];
for i = 1:length(season)
    month = X.month1(i);
    if ismember(month,in_season) == 1
        season(i) = 1;
    end
end
Y.season = season;    

%% Selects Vars
DATA04 = [Y.tempmax Y.soilmax Y.solar Y.eto Y.rhmax Y.windu Y.windv];
Ni04 = length(DATA04); %initial length
inds04 = ones(Ni04,1); %pre-allocate storage space for index vector

%% Select rows w/o missing values
for i = 1:Ni04
    nansum = sum(isnan(DATA04(i,:)));
    if nansum > 0
        inds04(i) = 0;
    end
end
inds04 = logical(inds04);

%% Create dataset w/o missing values
DATAX04 = DATA04(inds04,:);
Q.AQIg = Y.AQIg(inds04,:);
Q.season = Y.season(inds04,:);
Q.weekend = Y.weekend(inds04,:);
Q.month1 = Y.month1(inds04,:);
Q.day1 = Y.day1(inds04,:);
%% Select dates w/o missing values
DATEX04 = Y.datenum(inds04,:);

%% Get PCs from 2004 data
PC104 = DATAX04*PC1coef;
PC204 = DATAX04*PC2coef;
PC304 = DATAX04*PC3coef;

PCs04 = [PC104 PC204 PC304];



%% Cross-validation
% Make predictions:
predictions04 = predict(MdlQuadratic,PCs04);
correct04 = 0;
%correct2 = 0;
for i = 1:length(predictions04)
    if predictions04(i) == Q.AQIg(i)
        correct04 = correct04 + 1;
    end
 %   if classes(i) == Q.AQIg(i)
  %      correct2 = correct2 + 1;
   % end
end
p_correct04 = correct04/length(predictions04);
%p_correct2 = correct2/length(classes);
%% Histogram
figure(7); clf;
subplot(221); hist(P.AQIg); title('2001-2003 raw AQI group');
axis([1 4 0 800])
subplot(222); hist(predictions); title('2001-2003 AQI predictions');
axis([1 4 0 800])
subplot(223); hist(Q.AQIg); title('2004 raw AQI groups');
axis([1 4 0 300])
subplot(224); hist(predictions04); title('2004 AQI predictions');
axis([1 4 0 300])

%% Find type of misclassification
MC = zeros(4,4);
for i = 1:4
    for j = 1:4
        AQIr = i;
        AQIp = j;
        indr = find(P.AQIg == AQIr);
        indp = find(predictions == AQIp);
        ind_mc = ismember(indr,indp);
        MC(i,j) = sum(ind_mc);
        clear indr indp ind_mc
    end
end

MC_p = MC/sum(sum(MC)) * 100;


MC04 = zeros(4,4);
for i = 1:4
    for j = 1:4
        AQIr = i;
        AQIp = j;
        indr = find(Q.AQIg == AQIr);
        indp = find(predictions04 == AQIp);
        ind_mc = ismember(indr,indp);
        MC04(i,j) = sum(ind_mc);
        clear indr indp ind_mc
    end
end

MC04_p = MC04/sum(sum(MC04)) * 100;



% Low to high:
p_lh = sum([MC_p(1,2:4)])+sum([MC_p(2,3:4)])+ sum([MC_p(3,4)]);
p04_lh = sum([MC04_p(1,2:4)])+sum([MC04_p(2,3:4)])+sum([MC_p(3,4)]);

% High to low
p_hl = sum([MC_p(2,1)])+sum([MC_p(3,1:2)])+sum([MC_p(4,1:3)]);
p04_hl = sum([MC04_p(2,1)])+sum([MC04_p(3,1:2)])+sum([MC04_p(4,1:3)]);

% 2 or 3 to below
p_23l = sum([MC_p(2,1)])+sum([MC_p(3,1:2)]);
p04_23l = sum([MC04_p(2,1)])+sum([MC04_p(3,1:2)]);

% 4 to below
p_4l = sum([MC_p(4,1:3)]);
p04_4l = sum([MC04_p(4,1:3)]);

%% Scatterplots
figure(8); clf;
subplot(121);
indr1 = find(P.AQIg==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(P.AQIg==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(P.AQIg==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(P.AQIg==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of raw AQI groups 2001-2003');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

subplot(122);
indr1 = find(predictions==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(predictions==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(predictions==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(predictions==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of predicted AQI groups 2001-2003');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

figure(9); clf;
subplot(121);
indr1 = find(Q.AQIg==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(Q.AQIg==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(Q.AQIg==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(Q.AQIg==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of raw AQI groups 2004');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

subplot(122);
indr1 = find(predictions04==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(predictions04==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(predictions04==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(predictions04==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of predicted AQI groups 2004');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');


%% Logistic Regression
DATA_LR = [PC1 PC2 PC3 P.weekend P.month1 P.season];
DATA01F = [DATAX P.weekend P.month1 P.season];
DATA04_LR = [PC104 PC204 PC304 Q.weekend Q.month1 Q.season];
DATA04F = [DATAX04 Q.weekend Q.month1 Q.season];
AQIg = categorical(P.AQIg);
%get logistic model:
[B,dev,stats] = mnrfit(DATA_LR,AQIg);
[B2,dev2,stats2]=mnrfit(DATA01F,AQIg);
%predicted probabilities:
AQI_LR = mnrval(B,DATA_LR);
AQI01F = mnrval(B2,DATA01F);
AQI04_LR = mnrval(B,DATA04_LR);
AQI04F = mnrval(B2,DATA04F);
% convert to predicted AQI:
for i = 1:length(AQI_LR)
    AQI_LR_p(i) = find(AQI_LR(i,:)==max(AQI_LR(i,:)));
    AQI01F_P(i) = find(AQI01F(i,:)==max(AQI01F(i,:)));
end
for i = 1:length(AQI04_LR)
    AQI04_LR_p(i) = find(AQI04_LR(i,:)==max(AQI04_LR(i,:)));
    AQI04F_P(i) = find(AQI04F(i,:)==max(AQI04F(i,:)));
end
%2001:
correct_LR = 0;
correct_01F = 0;
for i = 1:length(AQI_LR_p)
    if AQI_LR_p(i) == P.AQIg(i)
        correct_LR = correct_LR + 1;
    end
    if AQI01F_P(i) == P.AQIg(i)
        correct_01F = correct_01F + 1;
    end
end
p_correct_LR = correct_LR/length(AQI_LR_p);
p_correct_01F = correct_01F/length(AQI_LR_p);
%2004:
correct04_LR = 0;
correct_04F = 0;
for i = 1:length(AQI04_LR_p)
    if AQI04_LR_p(i) == Q.AQIg(i)
        correct04_LR = correct04_LR + 1;
    end
    if AQI04F_P(i) == Q.AQIg(i)
        correct_04F = correct_04F + 1;
    end
end
p_correct04_LR = correct04_LR/length(AQI04_LR_p);
p_correct_04F = correct_04F/length(AQI04_LR_p);

%% Histogram
figure(10); clf;
subplot(221); hist(P.AQIg); title('2001-2003 raw AQI group');
axis([1 4 0 800])
subplot(222); hist(AQI_LR_p); title('2001-2003 AQI predictions');
axis([1 4 0 800])
subplot(223); hist(Q.AQIg); title('2004 raw AQI groups');
axis([1 4 0 300])
subplot(224); hist(AQI04_LR_p); title('2004 AQI predictions');
axis([1 4 0 300])

%% Find type of misclassification
MC_LR = zeros(4,4);
for i = 1:4
    for j = 1:4
        AQIr = i;
        AQIp = j;
        indr = find(P.AQIg == AQIr);
        indp = find(AQI_LR_p == AQIp);
        ind_mc = ismember(indr,indp);
        MC_LR(i,j) = sum(ind_mc);
        clear indr indp ind_mc
    end
end

MC_LR_p = MC_LR/sum(sum(MC_LR)) * 100;


MC04_LR = zeros(4,4);
for i = 1:4
    for j = 1:4
        AQIr = i;
        AQIp = j;
        indr = find(Q.AQIg == AQIr);
        indp = find(AQI04_LR_p == AQIp);
        ind_mc = ismember(indr,indp);
        MC04_LR(i,j) = sum(ind_mc);
        clear indr indp ind_mc
    end
end

MC04_LR_p = MC04_LR/sum(sum(MC04_LR)) * 100;



% Low to high:
p_lh_LR = sum([MC_LR_p(1,2:4)])+sum([MC_LR_p(2,3:4)])+ sum([MC_LR_p(3,4)]);
p04_lh_LR = sum([MC04_LR_p(1,2:4)])+sum([MC04_LR_p(2,3:4)])+sum([MC_LR_p(3,4)]);

% High to low
p_hl_LR = sum([MC_LR_p(2,1)])+sum([MC_LR_p(3,1:2)])+sum([MC_LR_p(4,1:3)]);
p04_hl_LR = sum([MC04_LR_p(2,1)])+sum([MC04_LR_p(3,1:2)])+sum([MC04_LR_p(4,1:3)]);

% 2 or 3 to below
p_23l_LR = sum([MC_LR_p(2,1)])+sum([MC_LR_p(3,1:2)]);
p04_23l_LR = sum([MC04_LR_p(2,1)])+sum([MC04_LR_p(3,1:2)]);

% 4 to below
p_4l_LR = sum([MC_LR_p(4,1:3)]);
p04_4l_LR = sum([MC04_LR_p(4,1:3)]);

%% Scatterplots
figure(11); clf;
subplot(121);
indr1 = find(P.AQIg==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(P.AQIg==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(P.AQIg==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(P.AQIg==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of raw AQI groups 2001-2003 (Logistic Regression)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

subplot(122);
indr1 = find(AQI_LR_p==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(AQI_LR_p==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(AQI_LR_p==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(AQI_LR_p==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of predicted AQI groups 2001-2003 (Logistic Regression)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

figure(12); clf;
subplot(121);
indr1 = find(Q.AQIg==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(Q.AQIg==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(Q.AQIg==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(Q.AQIg==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of raw AQI groups 2004 (Logistic Regression)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

subplot(122);
indr1 = find(AQI04_LR_p==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(AQI04_LR_p==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(AQI04_LR_p==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(AQI04_LR_p==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of predicted AQI groups 2004 (Logistic Regression)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

%% Plot PCs
figure(5); clf;
subplot(321); plot(DATEX,PC1); datetick('x','yyyy');
title('PC1 time series 2001-2003'); ylabel('PC1'); xlabel('Year');
subplot(323); plot(DATEX,PC2); datetick('x','yyyy');
title('PC2 time series 2001-2003'); ylabel('PC2'); xlabel('Year');
subplot(325); plot(DATEX,PC3); datetick('x','yyyy');
title('PC3 time series 2001-2003'); ylabel('PC3'); xlabel('Year');
%% Plot PCs from 2004
%figure(6); clf;
subplot(322); plot(DATEX04,PC104); datetick('x', 'yyyy');
title('PC1 times series 2004'); xlabel('Year'); ylabel('PC1');
subplot(324); plot(DATEX04,PC204); datetick('x','yyyy');
title('PC2 time series 2004'); xlabel('Year'); ylabel('PC2');
subplot(326); plot(DATEX04,PC304); datetick('x','yyyy');
title('PC3 time series 2004'); xlabel('Year'); ylabel('PC3');


%% Scatterplots
figure(13); clf;
subplot(131);
indr1 = find(P.AQIg==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(P.AQIg==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(P.AQIg==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(P.AQIg==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('Raw AQI groups 2001-2003');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

subplot(132);
indr1 = find(predictions==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(predictions==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(predictions==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(predictions==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('Predicted AQI groups 2001-2003 (DA)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

subplot(133);
indr1 = find(AQI_LR_p==1);
scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
indr2 = find(AQI_LR_p==2);
scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
indr3 = find(AQI_LR_p==3);
scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
indr4 = find(AQI_LR_p==4);
scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('Predicted AQI groups 2001-2003 (LR)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

% figure(13); clf;
% subplot(131);
% indr1 = find(Q.AQIg==1);
% scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
% indr2 = find(Q.AQIg==2);
% scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
% indr3 = find(Q.AQIg==3);
% scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
% indr4 = find(Q.AQIg==4);
% scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
% legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
% title('3D scatter plot of raw AQI groups 2004');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

% subplot(132);
% indr1 = find(AQI04_LR_p==1);
% scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
% indr2 = find(AQI04_LR_p==2);
% scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
% indr3 = find(AQI04_LR_p==3);
% scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
% indr4 = find(AQI04_LR_p==4);
% scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
% legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
% title('3D scatter plot of predicted AQI groups 2004 (LR)');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');




% %% Scatterplots
% figure(8); clf;
% subplot(121);
% indr1 = find(P.AQIg==1);
% scatter3(PC1(indr1),PC2(indr1),PC3(indr1),'k'); hold on;
% indr2 = find(P.AQIg==2);
% scatter3(PC1(indr2),PC2(indr2),PC3(indr2),'b^');
% indr3 = find(P.AQIg==3);
% scatter3(PC1(indr3),PC2(indr3),PC3(indr3),'g');
% indr4 = find(P.AQIg==4);
% scatter3(PC1(indr4),PC2(indr4),PC3(indr4),'r*');
% legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
% title('3D scatter plot of raw AQI groups 2001-2003');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');



% figure(14); clf;
% subplot(131);
figure(15); clf;
indr1 = find(Q.AQIg==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(Q.AQIg==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(Q.AQIg==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(Q.AQIg==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of raw AQI groups 2004');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

%subplot(132);
figure(16); clf;
indr1 = find(predictions04==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(predictions04==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(predictions04==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(predictions04==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of predicted AQI groups 2004 (DA)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

%subplot(133); clf;
figure(17); clf;
indr1 = find(AQI04_LR_p==1);
scatter3(PC104(indr1),PC204(indr1),PC304(indr1),'k'); hold on;
indr2 = find(AQI04_LR_p==2);
scatter3(PC104(indr2),PC204(indr2),PC304(indr2),'b^');
indr3 = find(AQI04_LR_p==3);
scatter3(PC104(indr3),PC204(indr3),PC304(indr3),'g');
indr4 = find(AQI04_LR_p==4);
scatter3(PC104(indr4),PC204(indr4),PC304(indr4),'r*');
legend('Group 1', 'Group 2', 'Group 3', 'Group 4');
title('3D scatter plot of predicted AQI groups 2004 (LR)');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
