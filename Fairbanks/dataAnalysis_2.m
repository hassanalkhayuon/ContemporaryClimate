% reading and analysing the data
% In this data analysis the linearly increasing trend in the mean after a
% certain date (chosen visually) is subtracted from the data.

% Analysis of data from Fairbanks, AK

clear all

% read the data file
% credit the authors of csvimport:  
% Ashish Sadanandan (2020). CSVIMPORT 
% (https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport), 
% MATLAB Central File Exchange. Retrieved June 12, 2020. 
[Yr, Tmp] = ...
    csvimport('mean-annual-temp_1906-2015.csv', ...
    'columns', {'Year', 'Temp'}, 'noHeader', false);
h_rawdata = figure('Name', 'Annual Temperature', 'Position', [100 100 400 500]);
plot(Yr, Tmp)
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
legend('data', [num2str(stepsize), '-yr averages'])
xlabel('year')
ylabel('Temperature')
title(['Temperature with ', num2str(stepsize), '-yr averages'])


% select the year at which the mean begins to increase linearly
yrchange = input('at what year does the mean begin to increase? ');
[mindiff, ichange] = min(abs(Yr-yrchange));
mean_before = mean(Tmp(1:ichange));

% assume that the mean increases linearly between yrchange and the end of
% the data set, and subtract this linear increase from the data
yrsafter = Yr(ichange:end);
Tmpafter = Tmp(ichange:end);
h_fit = figure('Name', 'Linearly Increasing Mean');
p = polyfit(yrsafter, Tmpafter,1);
f = polyval(p,yrsafter);
hold on
plot(Yr, Tmp)
%plot(yrsafter, Tmpafter)
plot(yrsafter, f, '--r')
plot(Yr(1:ichange), mean_before*ones(size(Yr(1:ichange))), '--r')
hold off
xlabel('Time (years)')
ylabel('Tmptation')
title('Annual Temperature with Mean Trends')

% now subtract the linearly increasing mean from the data
Tmp_nochange = Tmp;
Tmp_nochange(ichange:end) = Tmp(ichange:end) - f + mean_before;
h_check = figure('Name', 'Adjusted Data')
plot(Yr, Tmp_nochange)
xlabel('Time (year)')
ylabel('Temperature')
title('Increasing Mean Trend Removed')


% need to transform into a binary series "good" and "bad" or 1 and 0
% SELECT ONE OF THE DEFINITIONS BELOW
%----------------------------------------------------------------------
% definition 1: "Good" = above the mean of the entire data set
% YrMean = mean(Tmp_nochange);
% YrType = logical(Tmp_nochange >= YrMean); 
% delimiter = YrMean*ones(size(Tmp_nochange));
% delim = 'delim1';
%----------------------------------------------------------------------
% definition 2: "Good" = above the mean of the data set before climate
% change increases are observed
% YrMean_firsthalf = mean(Tmp_nochange(1:ichange));
% YrType = logical(Tmp_nochange >= YrMean_firsthalf);
% delimiter = YrMean_firsthalf*ones(size(Tmp_nochange));
% delim = 'delim2';
%----------------------------------------------------------------------
% definition 3: "Good" = above the mean of the first quarter of the data set
% YrMean_firstquarter = mean(Tmp_nochange(1:ceil(end/4)));
% YrType = logical(Tmp_nochange >= YrMean_firstquarter);
% delimiter = YrMean_firstquarter*ones(size(Tmp_nochange));
% delim = 'delim3';
%----------------------------------------------------------------------
% definition 4: "Good" = within a certain distance of the first half mean
YrMean_firsthalf = mean(Tmp_nochange(1:ichange));
YrMax_firsthalf = max(abs(Tmp_nochange-YrMean_firsthalf));
width = 0.2;
YrType = logical(abs(Tmp_nochange-YrMean_firsthalf) <= width*YrMax_firsthalf);
delimiter = [YrMean_firsthalf+width*YrMax_firsthalf; ...
    YrMean_firsthalf-width*YrMax_firsthalf]*ones(size(Tmp'));
delim = 'delim4';
%----------------------------------------------------------------------

h_delim = figure('Name', 'Annual Temperature With Thresholds');
clf
hold on
plot(Yr, Tmp_nochange)
plot(Yr, delimiter',':r')
hold off
xlabel('Year')
ylabel('Temperature')
title(['Annual Tmp Adj,', delim, ', Fairbanks'])

h_var = figure('Name', 'Temperature Variability with Thresholds');
clf
hold on
plot(Yr, Tmp_nochange-mean_before);
plot(Yr, delimiter'-mean_before,':r')
hold off
xlabel('Year')
ylabel('Temperature')
title(['Temperature Variability & Normal Range, Fairbanks'])


% figure(3)
% plot(YrType)
% %axis([0 imax -0.5 1.5])
% xlabel('Year')
% ylabel('Year Type')
% title('Above mean = good year, Below mean = bad year')

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
title(['PDF Year Type Interval, ', delim, ', Fairbanks'])

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
x = 0:2:10;
p = [0.3, 0.4, 0.5];
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
title('before CC')
subplot(1,2,2)
sortedSwInterval_after = sort(SwInterval_afterCC);
histogram(sortedSwInterval_after, 'Normalization', 'pdf')
hold on
x = 0:2:10;
p = [0.1, 0.3, 0.5];
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
title('with CC')










