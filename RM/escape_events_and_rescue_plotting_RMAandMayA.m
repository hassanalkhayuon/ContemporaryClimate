% Creating plot for escape events and rescue events for both the RMA model
% and the MayA model
load('escape_events_and_rescue_RMAandMayA.mat')
set(0,'defaulttextInterpreter','latex');

%% interpolation 
rho = linspace(0.001,0.99,50);
PP  = linspace(0.001,0.99,20);
ee_may = interp1(PP,escape_number_MayA,rho,'spline');
re_may = interp1(PP,rescue_number_MayA,rho,'spline');
ee_rma = interp1(PP,escape_number_RMA,rho, 'spline');
re_rma = interp1(PP,rescue_number_RMA,rho, 'spline');
%% plotting 
subplot(1,2,1)
fill_between(rho, ee_rma, re_rma,[])
subplot(1,2,2)
fill_between(rho, ee_may, re_may,[])

% blue = 0.10,0.40,0.69
% 
% green = 0.31,0.70,0.40
axis([0 1 0 1001])

set(gcf, 'renderer', 'painters')
