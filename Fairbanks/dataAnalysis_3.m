% reading and analysing the data
% In this data analysis the linearly increasing trend in the mean after a
% certain date (chosen visually) is subtracted from the data.

% Analysis of data from FAIRBANKS, AK

clear all

% read the data file
% credit the authors of csvimport:  
% Ashish Sadanandan (2020). CSVIMPORT 
% (https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport), 
% MATLAB Central File Exchange. Retrieved June 12, 2020. 
[Yr, Tmp] = ...
    csvimport('mean-annual-temp_1906-2015_modified.csv', ...
    'columns', {'Year', 'Temp'}, 'noHeader', false);
h_rawdata = figure('Name', 'Annual Temperature', 'Position', [100 100 400 500]);
plot(Yr, Tmp)
xleft = floor(Yr(1)/10)*10;
xright = ceil(Yr(end)/10)*10;
axis([xleft xright 0.8*min(Tmp) 1.2*max(Tmp)])
xlabel('Year')
ylabel('Total Tmptation (mm)')
title('Total Annual Temperature')
imax = length(Yr);

h_10yrAvg = figure('Name', '10-year Average');
l = length(Tmp);
for i = 10:l
    Tmp_avg(i-9) = mean(Tmp(i-9:i));
end
plot(Yr(10:end), Tmp_avg)
xlabel('Year')
ylabel('Tmptation (mm)')
title('10-Year Annual Average Temperature')

% Find the average Temperature in 10-year intervals, not a sliding window
h_steps = figure('Name', 'Climate Steps', 'Position', [100 100 400 500]);
TmpSteps = zeros(size(Tmp));   % step-function for Temperature
steps_begin = (floor(Yr(1)/10))*10;  % line up the steps with decades
steps_end = (ceil(Yr(end)/10))*10;
stepsize = 10;                       % number of years per step
sl = steps_begin + stepsize - Yr(1) + 1;  % number of years in first step
TmpSteps(1:sl) = mean(Tmp(1:sl));   % mean Temperature at first step
shift = stepsize-sl;                 % shift in the index since first step short
i = 2;
while i*stepsize+Yr(1) < Yr(end)
    TmpSteps((i-1)*stepsize+1-shift:i*stepsize-shift) = ...
        mean(Tmp((i-1)*stepsize+1-shift:i*stepsize-shift));
    i = i+1;
end
TmpSteps((i-1)*stepsize+1-shift:end) = ...
    mean(Tmp((i-1)*stepsize+1-shift:end));
plot(Yr, Tmp, 'o', Yr, TmpSteps, '-');
axis([xleft xright 0.8*min(Tmp) 1.2*max(Tmp)])
legend('data', [num2str(stepsize), '-yr averages'])
xlabel('year')
ylabel('Temperature')
title(['Temperature with ', num2str(stepsize), '-yr averages'])

% subtract the 10-year means from the data
h_DT = figure('Name', 'Detrended Data', 'Position', [100 100 400 500]);
Tmp_DT = Tmp - TmpSteps;
plot(Yr, Tmp_DT);
axis([xleft xright 1.2*min(Tmp_DT) 1.2*max(Tmp_DT)])
xlabel('year')
ylabel('detrended Temperature')
title('Annual Temperature with 10-Yr Averages Removed')

% select the year at which the mean begins to increase linearly
yrchange = input('at what year does the mean begin to increase? ');
[mindiff, ichange] = min(abs(Yr-yrchange));
mean_before = mean(Tmp(1:ichange));

% need to transform into a binary series "good" and "bad" or 1 and 0
% definition 4: "Good" = within a certain distance of the first half mean
YrMean = mean(Tmp_DT);
YrMax_before = max(abs(Tmp_DT(1:ichange) - YrMean));
width = 0.2;
YrType = logical(abs(Tmp_DT-YrMean) <= width*YrMax_before);
delimiter = [YrMean+width*YrMax_before; ...
    YrMean-width*YrMax_before]*ones(size(Tmp'));
delim = 'delim4';
%----------------------------------------------------------------------

h_delim = figure('Name', 'Annual Temperature With Thresholds', ...
    'Position', [100 100 400 500]);
clf
hold on
plot(Yr, Tmp_DT)
plot(Yr, delimiter',':r')
axis([xleft xright 1.2*min(Tmp_DT) 1.2*max(Tmp_DT)])
hold off
xlabel('Year')
ylabel('Temperature')
title(['Annual Temperature, Means Removed'])

% now find length of each string of type 1 or 0
SwitchMarker = [0; diff(YrType)];
numswitches = sum(abs(SwitchMarker))+1;
SwInterval = zeros(size(YrType(1:numswitches)));  % length of each switch interval
iswitch = 0;  % index at which switch occurs
intervalcounter = 0; % length of current type interval
nswitch = 0;  % number of switch
for i = 1:imax
    intervalcounter = intervalcounter + 1;
    if (SwitchMarker(i) ~= 0)
        iswitch = i;
        nswitch = nswitch+1;
        SwInterval(nswitch) = intervalcounter;
        intervalcounter = 0;
    end
end
% add on last set of years with no switch
SwInterval(nswitch+1) = intervalcounter;
SortedSwIs = sort(SwInterval);
% figure(4)
% plot(SortedSwIs)
% xlabel('interval number')
% ylabel('interval length')

% now find the probability for each string length - all data
h_hist = figure('Name', 'Frequency of Switching Intervals');
histogram(SortedSwIs, 'Normalization', 'pdf')
xlabel('interval length')
ylabel('frequency')

% plot several geometric distributions on the histogram - all data
x = 0:2:30;
h_geom = figure('Name', 'Geometric Distributions');
hold on
y0 = geopdf(x,0.1);
plot(x,y0,'s:')
y1 = geopdf(x,0.2);
plot(x,y1,'*:')
y2 = geopdf(x,0.3);
plot(x,y2,'o:')
hold off
legend('data', 'p=0.1', 'p=0.2', 'p=0.3')
title(['PDF Year Type Interval, ', delim, ', Burlington'])

% repeat switch interval calculations, but now split into 1st part of data
% set and climate-changed part of data set
startClimateChangeYr = yrchange;
beforeClimateChange = startClimateChangeYr-Yr(1)+1;
sumSwInterval = cumsum(SwInterval);
[M,I] = min(abs(sumSwInterval-beforeClimateChange));
SwInterval_beforeCC = SwInterval(1:I);
SwInterval_afterCC = SwInterval(I+1:end);

% create the histograms
figure(6)
subplot(1,2,1)
sortedSwInterval_before = sort(SwInterval_beforeCC);
histogram(sortedSwInterval_before, 'Normalization', 'pdf')
hold on
x = 0:2:15;
p = [0.1, 0.2, 0.3];
y0 = geopdf(x,p(1));
plot(x,y0,'s:')
y1 = geopdf(x,p(2));
plot(x,y1,'*:')
y2 = geopdf(x,p(3));
plot(x,y2,'o:')
hold off
xlabel('interval length')
ylabel('frequency')
legend('data', ['p=', num2str(p(1))], ['p=', num2str(p(2))], ...
    ['p=', num2str(p(3))])
title(['before ', num2str(yrchange)])
subplot(1,2,2)
sortedSwInterval_after = sort(SwInterval_afterCC);
histogram(sortedSwInterval_after, 'Normalization', 'pdf')
hold on
x = 0:2:15;
p = [0.1, 0.2, 0.3];
y0 = geopdf(x,p(1));
plot(x,y0,'s:')
y1 = geopdf(x,p(2));
plot(x,y1,'*:')
y2 = geopdf(x,p(3));
plot(x,y2,'o:')
hold off
xlabel('interval length')
ylabel('frequency')
legend('data', ['p=', num2str(p(1))], ['p=', num2str(p(2))], ...
    ['p=', num2str(p(3))])
title(['after ', num2str(yrchange)])










