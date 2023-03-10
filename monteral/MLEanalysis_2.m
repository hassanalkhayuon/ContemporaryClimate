% reading and analysing the data
% In this data analysis the linearly increasing trend in the mean after a
% certain date (chosen visually) is subtracted from the data.

% Analysis of data from Montreal, PQ

clear all

% read the data file
% credit the authors of csvimport:  
% Ashish Sadanandan (2020). CSVIMPORT 
% (https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport), 
% MATLAB Central File Exchange. Retrieved June 12, 2020. 
[year month MthPrecip] = ...
    csvimport('Montreal-precip-data_mthly_MatlabFormat.csv', ...
    'columns', {'Year', 'Month', 'TtlP'}, 'noHeader', false);
YrMth = year + (month-1)/12;
figure(1)
plot(YrMth, MthPrecip)
xlabel('Year')
ylabel('Precipitation (mm)')
title('Total Monthly Precipitation')

% find the yearly precipitation
imax = floor(length(MthPrecip)/12);
i = 0;
Precip = zeros(imax,1);
Yr = zeros(size(Precip));
Yr(1) = floor(YrMth(1));
Precip(1) = sum(MthPrecip((1:12)));
for i = 2:imax
    Precip(i) = sum(MthPrecip((i-1)*12+1:i*12));
    Yr(i) = Yr(i-1)+1;
end
xleft = floor(Yr(1)/10)*10;
xright = ceil(Yr(end)/10)*10;
figure(2)
plot(Yr, Precip)
xlabel('month counter')
ylabel('precipitation (mm)')
title('yearly precipitation')

% h_10yrAvg = figure('Name', '10-year Average');
l = length(Precip);
for i = 10:l
    Precip_avg(i-9) = mean(Precip(i-9:i));
end
% plot(Yr(10:end), Precip_avg)
% xlabel('Year')
% ylabel('Preciptation (mm)')
% title('10-Year Annual Average Precipitation')

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
plot(Yr, Precip, Yr, PrecipSteps);
axis([xleft xright 0.8*min(Precip) 1.2*max(Precip)])
legend('data', [num2str(stepsize), '-yr averages'])
xlabel('year')
ylabel('precipitation')
title(['Precipitation with ', num2str(stepsize), '-yr averages'])

% subtract the 10-year means from the data
% h_DT = figure('Name', 'Detrended Data', 'Position', [100 100 400 500]);
Precip_DT = Precip - PrecipSteps;
% plot(Yr, Precip_DT);
% axis([xleft xright 1.2*min(Precip_DT) 1.2*max(Precip_DT)])
% xlabel('year')
% ylabel('detrended precipitation')
% title('Annual Precipitation with 10-Yr Averages Removed')

% select the year at which the mean begins to increase linearly
yrchange = input('at what year does the mean begin to increase? ');
[mindiff, ichange] = min(abs(Yr-yrchange));
mean_before = mean(Precip(1:ichange));

% need to transform into a binary series "good" and "bad" or 1 and 0
% definition 4: "Good" = within a certain distance of the first half mean
YrMean = mean(Precip_DT);
YrMax_before = max(abs(Precip_DT(1:ichange) - YrMean));
width = 0.2;
YrType = logical(abs(Precip_DT-YrMean) <= width*YrMax_before);
delimiter = [YrMean+width*YrMax_before; ...
    YrMean-width*YrMax_before]*ones(size(Precip'));
delim = 'delim4';
%----------------------------------------------------------------------

h_delim = figure('Name', 'Annual Precipitation With Thresholds', ...
    'Position', [100 100 400 500]);
clf
hold on
plot(Yr, Precip_DT)
plot(Yr, delimiter',':r')
axis([xleft xright 1.2*min(Precip_DT) 1.2*max(Precip_DT)])
hold off
xlabel('Year')
ylabel('Precipitation')
title(['Annual Precipitation, Means Removed'])

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
% h_hist = figure('Name', 'Frequency of Switching Intervals');
% histogram(SortedSwIs, 'Normalization', 'pdf')
% xlabel('interval length')
ylabel('frequency')

% repeat switch interval calculations, but now split into 1st part of data
% set and climate-changed part of data set
startClimateChangeYr = yrchange;
beforeClimateChange = startClimateChangeYr-Yr(1)+1;
sumSwInterval = cumsum(SwInterval);
[M,I] = min(abs(sumSwInterval-beforeClimateChange));
SwInterval_beforeCC = SwInterval(1:I);
SwInterval_afterCC = SwInterval(I+1:end);
sortedSwInterval_before = sort(SwInterval_beforeCC);
sortedSwInterval_after = sort(SwInterval_afterCC);

% use Maximum Likelihood Estimation to compute the most likely value of p
% analytically (see Hassan's notes)
% All Data
pvec = 0:0.01:1;
Lfn = (nswitch+1)*log(pvec) + log(1-pvec)*sumSwInterval(end);
hMLE = figure('Name', 'MLE', 'Position', [110 110 400 500]);
plot(pvec,Lfn, 'b', 'LineWidth', 2)
pMLE = (nswitch+1)/((nswitch+1) + sumSwInterval(end));
LMLE = (nswitch+1)*log(pMLE) + log(1-pMLE)*sumSwInterval(end);
hold on
scatter(pMLE, LMLE, 'k', 'filled')
plot([pMLE pMLE], [min(get(gca,'YLim')) LMLE], 'k--');
plot([0 pMLE], [LMLE LMLE], 'k--');
text(1.1*pMLE,0.9*LMLE, ['p=',num2str(pMLE)])
hold off
xlabel('geometric parameter, p')
ylabel('Likelihood')

% Before Climate Change Data
n = length(SwInterval_beforeCC);
sum_beforeCC = cumsum(SwInterval_beforeCC);
Lfn_before = n*log(pvec) + log(1-pvec)*sum_beforeCC(end);
hMLE_before = figure('Name', 'MLE_before', 'Position', [120 120 400 500]);
plot(pvec,Lfn_before, 'b', 'LineWidth', 2)
pMLE_before = n/(n + sum_beforeCC(end));
LMLE_before = n*log(pMLE_before) + log(1-pMLE_before)*sum_beforeCC(end);
hold on
scatter(pMLE_before, LMLE_before, 'k', 'filled')
plot([pMLE_before pMLE_before], [min(get(gca,'YLim')) LMLE_before], 'k--');
plot([0 pMLE_before], [LMLE_before LMLE_before], 'k--');
text(1.1*pMLE_before,0.9*LMLE_before, ['p=',num2str(pMLE_before)])
hold off
xlabel('geometric parameter, p')
ylabel('Likelihood')
title(['before ', num2str(yrchange)])

% During Climate Change Data
n = length(SwInterval_afterCC);
sum_afterCC = cumsum(SwInterval_afterCC);
Lfn_after = n*log(pvec) + log(1-pvec)*sum_afterCC(end);
hMLE_after = figure('Name', 'MLE_after', 'Position', [120 120 400 500]);
plot(pvec,Lfn_after, 'b', 'LineWidth', 2)
pMLE_after = n/(n + sum_afterCC(end));
LMLE_after = n*log(pMLE_after) + log(1-pMLE_after)*sum_afterCC(end);
hold on
scatter(pMLE_after, LMLE_after, 'k', 'filled')
plot([pMLE_after pMLE_after], [min(get(gca,'YLim')) LMLE_after], 'k--');
plot([0 pMLE_after], [LMLE_after LMLE_after], 'k--');
text(1.1*pMLE_after,0.9*LMLE_after, ['p=',num2str(pMLE_after)])
hold off
xlabel('geometric parameter, p')
ylabel('Likelihood')
title(['after ', num2str(yrchange)])

% plot the histogram - all data - and three geometric distributions: the
% pMLE, and two "boundary" values
x = 0:2:30;
% p = [pMLE, 0.2, 0.5];
% h_geom = figure('Name', 'Geometric Distributions', ...
%     'Position', [120 120 400 500]);
% hold on
% y0 = geopdf(x,p(1));
% plot(x,y0,'s:')
% y1 = geopdf(x,p(2));
% plot(x,y1,'*:')
% y2 = geopdf(x,p(3));
% plot(x,y2,'o:')
% % hold off
% legend('data', ['p (MLE) = ', num2str(p(1))], ['p = ', num2str(p(2))], ...
%     ['p = ', num2str(p(2))])
% title(['PDF Year Type Interval, ', delim, ', Burlington'])

% plot the histogram - all data - and the most likely geometric distribution
x = 0:2:30;
h_geom2 = figure('Name', 'Geometric Distributions', ...
    'Position', [120 120 400 500]);
hold on
y0 = geopdf(x,pMLE);
plot(x,y0,'s:')
hold off
legend('data', ['p (MLE) = ', num2str(pMLE)])
title(['PDF Year Type Interval, ', delim, ', Burlington'])


% create the histograms
% h_histplotsBA = figure('Name', 'BA Histograms with MLE geom distns');
% subplot(1,2,1)
% sortedSwInterval_before = sort(SwInterval_beforeCC);
% histogram(sortedSwInterval_before, 'Normalization', 'pdf')
% hold on
% x = 0:2:10;
% y0 = geopdf(x,pMLE_before);
% plot(x,y0,'s:')
% hold off
% xlabel('interval length')
% ylabel('frequency')
% legend('data', ['p (MLE) = ', num2str(pMLE_before)])
% title(['before ', num2str(yrchange)])
% subplot(1,2,2)
% sortedSwInterval_after = sort(SwInterval_afterCC);
% histogram(sortedSwInterval_after, 'Normalization', 'pdf')
% hold on
% x = 0:2:10;
% y0 = geopdf(x,pMLE_after);
% plot(x,y0,'s:')
% hold off
% xlabel('interval length')
% ylabel('frequency')
% legend('data', ['p (MLE) = ', num2str(pMLE_after)])
% title(['after ', num2str(yrchange)])

% create the histograms
% h_histplotsBA_multi = figure('Name', 'BA Histograms with multiple geom distns');
% subplot(1,2,1)
% sortedSwInterval_before = sort(SwInterval_beforeCC);
% histogram(sortedSwInterval_before, 'Normalization', 'pdf')
% hold on
% x = 0:2:10;
% p = [pMLE_before, 0.2, 0.5];
% y0 = geopdf(x,p(1));
% plot(x,y0,'s:')
% y1 = geopdf(x,p(2));
% plot(x,y1,'*:')
% y2 = geopdf(x,p(3));
% plot(x,y2,'o:')
% hold off
% xlabel('interval length')
% ylabel('frequency')
% legend('data', ['p (MLE) = ', num2str(p(1))], ['p = ', num2str(p(2))], ...
%     ['p = ', num2str(p(3))]) 
% title(['before ', num2str(yrchange)])
% subplot(1,2,2)
% sortedSwInterval_after = sort(SwInterval_afterCC);
% histogram(sortedSwInterval_after, 'Normalization', 'pdf')
% hold on
% x = 0:2:10;
% p = [pMLE_after, 0.2, 0.5];
% y0 = geopdf(x,p(1));
% plot(x,y0,'s:')
% y1 = geopdf(x,p(2));
% plot(x,y1,'*:')
% y2 = geopdf(x,p(3));
% plot(x,y2,'o:')
% hold off
% xlabel('interval length')
% ylabel('frequency')
% legend('data', ['p (MLE) = ', num2str(p(1))], ['p = ', num2str(p(2))], ...
%     ['p = ', num2str(p(3))])
% title(['after ', num2str(yrchange)])

%-------------------------------------------------------------------
% Composite plots for paper

h_allplots = figure('Name', 'Composite plot for paper', 'Position', ...
    [100 400 900 900]);
subplot(2,3,1)
plot(Yr, Precip, Yr, PrecipSteps);
axis([xleft xright 0.8*min(Precip) 1.2*max(Precip)])
legend('data', [num2str(stepsize), '-yr averages'])
xlabel('year')
ylabel('precipitation')
title(['Precip w/ ', num2str(stepsize), '-yr avgs'])
subplot(2,3,4)
hold on
plot(Yr, Precip_DT)
plot(Yr, delimiter',':r')
axis([xleft xright 1.2*min(Precip_DT) 1.2*max(Precip_DT)])
hold off
xlabel('Year')
ylabel('Precipitation')
title(['Ann Precip, Means Removed'])
subplot(2,3,[2 3])
histogram(sortedSwInterval_before, 'Normalization', 'pdf')
hold on
x = 0:2:10;
y0 = geopdf(x,pMLE_before);
plot(x,y0,'s:')
hold off
xlabel('interval length')
ylabel('frequency')
legend('data', ['p (MLE) = ', num2str(pMLE_before)])
title(['before ', num2str(yrchange)])
subplot(2,3,[5 6])
histogram(sortedSwInterval_after, 'Normalization', 'pdf')
hold on
x = 0:2:10;
y0 = geopdf(x,pMLE_after);
plot(x,y0,'s:')
hold off
xlabel('interval length')
ylabel('frequency')
legend('data', ['p (MLE) = ', num2str(pMLE_after)])
title(['after ', num2str(yrchange)])

%-------------------------------------------------------------------
% Composite plots for paper

h_allplots = figure('Name', 'Composite plot for paper', 'Position', ...
    [100 400 900 900]);
subplot(2,3,1)
plot(Yr, Precip, Yr, PrecipSteps);
axis([xleft xright 0.8*min(Precip) 1.2*max(Precip)])
legend('data', [num2str(stepsize), '-yr averages'])
xlabel('year')
ylabel('Precipitation')
title(['Precip w/ ', num2str(stepsize), '-yr avgs'])
subplot(2,3,4)
hold on
plot(Yr, Precip_DT)
plot(Yr, delimiter',':r')
axis([xleft xright 1.2*min(Precip_DT) 1.2*max(Precip_DT)])
hold off
xlabel('Year')
ylabel('Precipitation')
title(['Ann Precip, Means Removed'])
subplot(1,3,[2 3])
histogram(SortedSwIs, 'Normalization', 'pdf')
hold on
x = 0:2:10;
y0 = geopdf(x,pMLE);
plot(x,y0,'s:')
hold off
xlabel('interval length')
ylabel('frequency')
legend('data', ['p (MLE) = ', num2str(pMLE)])
title(['all years'])






