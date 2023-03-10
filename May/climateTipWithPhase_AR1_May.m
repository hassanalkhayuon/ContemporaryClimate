% climateTipWithPhase_AR1_RM
% simulate one time series that tips to exctinction. 
% del is going to be fixes through out, and R is going to be varied btween
% Rstar and Rend acording to normlized autoregressve model.
% r(t) = phi*r(t-1) + eps_t where time in years. 

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\10-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmax   = 3.5;       % /yr.              Prey intrinsic growth rate
Rmin    = 2;
Rmean   = (Rmax+Rmin)/2;
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);

% rng(86) 
% rng(1150) %for Rstar = 2.5
%% forcing 
Tend  = 1000;
climate(1) = rand;
climate(2) = rand;
climate(3) = rand;
phi1 = 0.5;
phi2 = 0.5;
for ind_t = 3: Tend 
    climate(ind_t) = randn + ...
        phi1*climate(ind_t-1) + ...
        phi2*climate(ind_t-2);
end
% normlise climate
normlised_climate = (climate - min(climate))/(max(climate) - min(climate));
modifid_climate = (Rmax-Rmin)*normlised_climate + Rmin;
figure(1)
subplot(6,12,[1:7,13:19,25:31])
hold on
plot([0:1:Tend-1],modifid_climate)
box on
axis([0 Tend 1 3.1])
hA=gca;
set(gca,'XMinorTick','on','YMinorTick','on')
hA.XAxis.MinorTickValues = 0:5:100;
xticks([0 20 40 60 80 100])
xticklabels([])
hA.YAxis.MinorTickValues = 1:.1:3;
yticks([1.5 2 2.5 3 3.5 4])
yticklabels({1.5 2 2.5 3 3.5})
plot([0 200],[Rmin Rmin],':k','LineWidth',1)
plot([0 200],[Rmax Rmax],':k','LineWidth',1)
set(gca,'FontSize',20)
ylabel('$R(t)$','Rotation',0)

%% time sieries
mult     = 1000;
initcond = [11,0.004];
for ind_t = 1:Tend
    R = modifid_climate(ind_t);
    odefun   = @(t,var)Mayodefun_cli(var,R);
    [t,var]  = ode45(odefun,[ind_t-1:0.1:ind_t],initcond,opts);
    initcond = var(end,:);
    figure(10);
    hold on
    plot(t, var(:,2), '-k')
end

% subplot(6,12,[37:43,49:55,61:67])
% % hold on
% plot(...
%     t,var(:,1),'-k','LineWidth',2)
% plot(...
%     t,mult*var(:,2),'-','LineWidth',2,'Color',[.65 .65 .65])
% box on
% axis([0 Tend -2 30])
% hA=gca;
% set(gca,'XMinorTick','on','YMinorTick','on')
% hA.XAxis.MinorTickValues = 0:5:100;
% xticks([0 20 40 60 80 100])
% xticklabels([0 20 40 60 80])
% hA.YAxis.MinorTickValues = -2:1:45;
% yticks([0 10 20 30])
% yticklabels([0 10 20])
% set(gca,'FontSize',20)
% xlabel('$t$')
% ylabel('$\tilde{P},N$','Rotation',90)
% % var(end,1)
% % crit(I) = var(end,1);
% % end
% % 
%% functions

%the ode fuction of May model
function [dvar] = Mayodefun_cli(var,R)

%parameters
q       = 205;
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.


N = var(1);
P = var(2);



dN = R * N *(1-((C/R)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = s*P*( 1-((q*P)/(N+eps)) );

dvar = [dN;dP];

end

