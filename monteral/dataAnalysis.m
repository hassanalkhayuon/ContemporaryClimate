% reading and analysing the Pelly Ranch precipitation data

% read the data file
% credit the authors of csvimport:  
% Ashish Sadanandan (2020). CSVIMPORT 
% (https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport), 
% MATLAB Central File Exchange. Retrieved June 12, 2020. 
[year month MthPrecip] = ...
    csvimport('Montreal-precip-data_mthly_MatlabFormat.csv', ...
    'columns', {'Year', 'Month', 'TtlP'}, 'noHeader', false);
YrMth = year + month/12;
figure(1)
plot(YrMth, MthPrecip)
xlabel('Year')
ylabel('Precipitation (mm)')
title('Total Monthly Precipitation')

% find the yearly precipitation
i = 0;
YrWindowPrecip = zeros(size(MthPrecip(12:end)));
imax = length(YrWindowPrecip);
for i = 1:imax
    Yr(i) = sum(MthPrecip(i:i+11));
end
figure(2)
plot(Yr)
xlabel('month counter')
ylabel('precipitation (mm)')
title('yearly precipitation (mthly summed over 12 mth window)')

% need to transform into a binary series "good" and "bad" or 1 and 0
% SELECT ONE OF THE DEFINITIONS BELOW
%----------------------------------------------------------------------
% definition 1: "Good" = above the mean of the entire data set
% YrMean = mean(Yr);
% YrType = logical(Yr >= YrMean); 
% delimiter = YrMean*ones(size(Yr));
%----------------------------------------------------------------------
% definition 2: "Good" = above the mean of the first half of the data set
% YrMean_firsthalf = mean(Yr(1:ceil(end/2)));
% YrType = logical(Yr >= YrMean_firsthalf);
% delimiter = YrMean_firsthalf*ones(size(Yr));
%----------------------------------------------------------------------
% definition 3: "Good" = above the mean of the first quarter of the data set
YrMean_firstquarter = mean(Yr(1:ceil(end/4)));
YrType = logical(Yr >= YrMean_firstquarter);
delimiter = YrMean_firstquarter*ones(size(Yr));
%----------------------------------------------------------------------
% definition 4: "Good" = within a certain distance of the first half mean
% YrMean_firsthalf = mean(Yr(1:ceil(end/2)));
% YrMax_firsthalf = max(abs(Yr-YrMean_firsthalf));
% width = 0.2;
% YrType = logical(abs(Yr-YrMean_firsthalf) <= width*YrMax_firsthalf);
% delimiter = [YrMean_firsthalf+width*YrMax_firsthalf; ...
%     YrMean_firsthalf-width*YrMax_firsthalf]*ones(size(Yr));
%----------------------------------------------------------------------

figure(2)
hold on
plot(delimiter',':r')
hold off
figure(3)
plot(YrType)
axis([0 imax -0.5 1.5])
xlabel('month counter')
ylabel('year type')
title('yearly window above or below mean')

% now find length of each string of type 1 or 0
SwitchMarker = [0 diff(YrType)];
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
x = 0:2:50;
figure(5)
hold on
y0 = geopdf(x,0.05);
plot(x,y0,'s:')
y1 = geopdf(x,0.1);
plot(x,y1,'*:')
y2 = geopdf(x,0.2);
plot(x,y2,'o:')
hold off
legend('data', 'p=0.05', 'p=0.1', 'p=0.2')
title('PDF Year Type Interval, Sliding 12-month window')







