% reading and analysing the Fairbanks temperature data

clear all

% read the data file
% credit the authors of csvimport:  
% Ashish Sadanandan (2020). CSVIMPORT 
% (https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport), 
% MATLAB Central File Exchange. Retrieved June 12, 2020. 
[Yr, Tmp] = ...
    csvimport('mean-annual-temp_1904-2015.csv', ...
    'columns', {'Year', 'Temp'}, 'noHeader', false);
figure(1)
plot(Yr, Tmp)
xlabel('Year')
ylabel('Mean Temperature (F)')
title('Annual Mean Temperature')
imax = length(Yr);

% need to transform into a binary series "good" and "bad" or 1 and 0
% SELECT ONE OF THE DEFINITIONS BELOW
%----------------------------------------------------------------------
% definition 1: "Good" = above the mean of the entire data set
% YrMean = mean(Tmp);
% YrType = logical(Tmp >= YrMean); 
% delimiter = YrMean*ones(size(Tmp));
%----------------------------------------------------------------------
% definition 2: "Good" = above the mean of the first half of the data set
% YrMean_firsthalf = mean(Tmp(1:ceil(end/2)));
% YrType = logical(Tmp >= YrMean_firsthalf);
% delimiter = YrMean_firsthalf*ones(size(Tmp));
%----------------------------------------------------------------------
% definition 3: "Good" = above the mean of the first quarter of the data set
% YrMean_firstquarter = mean(Tmp(1:ceil(end/4)));
% YrType = logical(Tmp >= YrMean_firstquarter);
% delimiter = YrMean_firstquarter*ones(size(Tmp));
%----------------------------------------------------------------------
% definition 4: "Good" = within a certain distance of the first half mean
YrMean_firsthalf = mean(Tmp(1:ceil(end/2)));
YrMax_firsthalf = max(abs(Tmp-YrMean_firsthalf));
width = 0.2;
YrType = logical(abs(Tmp-YrMean_firsthalf) <= width*YrMax_firsthalf);
delimiter = [YrMean_firsthalf+width*YrMax_firsthalf; ...
    YrMean_firsthalf-width*YrMax_firsthalf]*ones(size(Tmp'));
%----------------------------------------------------------------------

figure(2)
hold on
plot(Yr, Tmp)
plot(Yr, delimiter',':r')
hold off
xlabel('Year')
ylabel('Temperature')

figure(3)
plot(YrType)
%axis([0 imax -0.5 1.5])
xlabel('Year')
ylabel('Year Type')
title('Above mean = good year, Below mean = bad year')

% now find length of each string of type 1 or 0
SwitchMarker = [0; diff(YrType)];
numswitches = sum(abs(SwitchMarker));
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
SortedSwIs = sort(SwInterval);
figure(4)
plot(SortedSwIs)
xlabel('interval number')
ylabel('interval length')

% now find the probability for each string length
figure(5)
histogram(SortedSwIs, 'Normalization', 'pdf')
xlabel('interval length')
ylabel('frequency')

% plot several geometric distributions on the histogram
x = 0:2:30;
figure(5)
hold on
y0 = geopdf(x,0.1);
plot(x,y0,'s:')
y1 = geopdf(x,0.2);
plot(x,y1,'*:')
y2 = geopdf(x,0.3);
plot(x,y2,'o:')
hold off
legend('data', 'p=0.1', 'p=0.2', 'p=0.3')
title('PDF Year Type Interval')







