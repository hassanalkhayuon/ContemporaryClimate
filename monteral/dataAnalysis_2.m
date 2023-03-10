% reading and analysing the data
% In this data analysis the linearly increasing trend in the mean after a
% certain date (chosen visually) is subtracted from the data.

clear all

% read the data file
% credit the authors of csvimport:  
% Ashish Sadanandan (2020). CSVIMPORT 
% (https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport), 
% MATLAB Central File Exchange. Retrieved June 12, 2020. 
[year month MthPrecip] = ...
    csvimport('Montreal-precip-data_mthly_MatlabFormat.csv', ...
    'columns', {'Year', 'Month', 'TtlP'}, 'noHeader', false);
YrMth = year + month/12;
h_precip = figure('Name', 'Monthly Precipitation');
plot(YrMth, MthPrecip)
xlabel('Year')
ylabel('Precipitation (mm)')
title('Total Monthly Precipitation (missing data = 0)')

% convert to annual data
Yr = year(1):year(end);
imax = length(Yr);
for i = 1:imax
    Precip(i) = sum(MthPrecip((i-1)*12+1:i*12));
end
h_Precip = figure('Name', 'Annual Precipitation', ...
    'Position', [100 100 400 500]);
plot(Yr, Precip)
xleft = floor(Yr(1)/10)*10;
xright = ceil(Yr(end)/10)*10;
axis([xleft xright 0.8*min(Precip) 1.2*max(Precip)])
xlabel('Year')
ylabel('Precipitation (mm)')
title('Annual Precipitation')

h_10yrAvg = figure('Name', '10-year Average', ...
    'Position', [100 100 400 500]);
l = length(Precip);
for i = 10:l
    Precip_avg(i-9) = mean(Precip(i-9:i));
end
scatter(Yr(10:end), Precip_avg, 5, 'filled')
axis([xleft xright 0.8*min(Precip) 1.2*max(Precip)])
xlabel('Year')
ylabel('Preciptation (mm)')
title('10-Year Annual Average Precipitation')

% Find the average precipitation in 10-year intervals, not a sliding window
h_steps = figure('Name', 'Climate Steps', 'Position', [100 100 400 500]);
PrecipSteps = zeros(size(Precip));   % step-function for precipitation
steps_begin = (floor(Yr(1)/10))*10;  % line up the steps with decades
steps_end = (ceil(Yr(end)/10))*10;
stepsize = 10;                       % number of years per step
sl = steps_begin + stepsize - Yr(1) + 1;  % number of years in first step
PrecipSteps(1:sl) = mean(Precip(1:sl));   % mean precipitation at first step
shift = stepsize-sl;                 % shift in the index since first step short
i = 2;
while i*stepsize+Yr(1) < Yr(end)
    PrecipSteps((i-1)*stepsize+1-shift:i*stepsize-shift) = ...
        mean(Precip((i-1)*stepsize+1-shift:i*stepsize-shift));
    i = i+1;
end
PrecipSteps((i-1)*stepsize+1-shift:end) = ...
    mean(Precip((i-1)*stepsize+1-shift:end));
plot(Yr, Precip, 'o', Yr, PrecipSteps, '-');
axis([xleft xright 0.8*min(Precip) 1.2*max(Precip)])
legend('data', [num2str(stepsize), '-yr averages'])
xlabel('year')
ylabel('precipitation')
title(['Precipitation with ', num2str(stepsize), '-yr averages'])

% subtract this step-wise average from the data to de-trend it
h_DT = figure('Name', 'Detrended Data', 'Position', [100 100 400 500]);
Precip_DT = Precip - PrecipSteps;
plot(Yr, Precip_DT);
axis([xleft xright 1.2*min(Precip_DT) 1.2*max(Precip_DT)])
xlabel('year')
ylabel('detrended precipitation')
title('Annual Precipitation with Means Removed')

% select the year after which the amplitude begins to increase
yrchange = input('after what year does the amplitude begin to increase? ');
[mindiff, ichange] = min(abs(Yr-yrchange));

% need to transform into a binary series "good" and "bad" or 1 and 0
% definition 4: "Good" = within a certain distance of the mean
YrMean = mean(Precip_DT);
YrMax_before = max(abs(Precip_DT(1:ichange) - YrMean));
width = 0.2;
YrType = logical(abs(Precip_DT-YrMean) <= width*YrMax_before);
delimiter = [YrMean+width*YrMax_before; ...
    YrMean-width*YrMax_before]*ones(size(Precip));
delim = 'delim4';
%----------------------------------------------------------------------

h_delim = figure('Name', 'Annual Precipitation With Thresholds', ...
    'Position', [100 100 400 500]);
clf
hold on
plot(Yr, Precip_DT)
plot(Yr, delimiter',':r')
hold off
axis([xleft xright 1.2*min(Precip_DT) 1.2*max(Precip_DT)])
xlabel('Year')
ylabel('Precipitation')
title(['Annual Precipitation with Means Removed'])

% now find length of each string of type 1 or 0
SwitchMarker = [0; (diff(YrType))'];
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

% now find the probability for each string length - all data
h_hist = figure('Name', 'Frequency of Switching Intervals');
histogram(SortedSwIs, 'Normalization', 'pdf')
xlabel('interval length')
ylabel('frequency')

% plot several geometric distributions on the histogram - all data
x = 0:2:15;
h_geom = figure('Name', 'Geometric Distributions');
histogram(SortedSwIs, 'Normalization', 'pdf')
hold on
y0 = geopdf(x,0.1);
plot(x,y0,'s:')
y1 = geopdf(x,0.2);
plot(x,y1,'*:')
y2 = geopdf(x,0.3);
plot(x,y2,'o:')
hold off
legend('data', 'p=0.1', 'p=0.2', 'p=0.3')
title(['PDF Year Type Interval, ', delim, ', Montreal'])

% repeat switch interval calculations, but now split into 1st part of data
% set and climate-changed part of data set
startClimateChangeYr = yrchange;
beforeClimateChange = startClimateChangeYr-Yr(1)+1;
sumSwInterval = cumsum(SwInterval);
[M,I] = min(abs(sumSwInterval-beforeClimateChange));
SwInterval_beforeCC = SwInterval(1:I);
SwInterval_afterCC = SwInterval(I+1:end);

% create the histograms
h_geom_before_after = figure('Name', 'Before/After Geometric Distributions');
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
title(['before ', num2str(yrchange)])
subplot(1,2,2)
sortedSwInterval_after = sort(SwInterval_afterCC);
histogram(sortedSwInterval_after, 'Normalization', 'pdf')
hold on
x = 0:2:10;
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










