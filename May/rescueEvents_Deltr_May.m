% to compute how many rescue events
warning off
% clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB\10_Phase_Sensitivity\P-tippingPaper\may_model')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmid = 2.65;
amp = [0.6 0.65 0.7];

q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(4)   % random seeding

%% Simulations

Tend      = 5000;
RR        = 1;  %avreage length of Type-L/H period
PnBin     = .3;
mult      = 1000;


K = 10000;

t_tip   = NaN(1,K);
R_aft   = NaN(1,K);
R_bef   = NaN(1,K);
var_bef = NaN(2,K);
resecueEvents = zeros(K,3);

for ind_p = 1:3
    
    Rstar  = Rmid + amp(ind_p);
    Rend   = Rmid - amp(ind_p);
    
    ind_sim = 0;
    while (ind_sim < K)
        
        vars_cl = {'T','ind_T','t','var','R'};
        clear(vars_cl{:})
        
        ind_sim = ind_sim + 1;
        
        T(1)           = nbininv(rand,RR,PnBin);
        R(1)           = NaN;
        initcond(:,1)  = [3;0.002];
        while (T(1) == 0)
            T(1) = nbininv(rand,RR,PnBin);
        end
        
        ind_T = 2;
        
        NORM = inf;
        Ttip = 0;
        while and(sum(T)<Tend , NORM>1e-3)
            T(ind_T) = nbininv(rand,RR,PnBin);
            while (T(ind_T) == 0)
                T(ind_T) = nbininv(rand,RR,PnBin);
            end
            R(ind_T) = Rend + rand*(Rstar - Rend);
            odefun   = @(t,var)Mayodefun(var,R(ind_T));
            tspan    = [sum(T(1:ind_T - 1)),sum(T(1:ind_T))];
            [t,var]  = ode45(odefun,tspan,initcond(:,ind_T-1),opts);
            initcond(:,ind_T) = var(end,:);
            NORM = norm(initcond(:,ind_T));
            
            [~,var_ersc]  = ode45(odefun,[tspan(1),tspan(1)+50],initcond(:,ind_T-1),opts);
            rescNorm = norm(var_ersc(end,:));
            
            if (NORM > 1e-3 && rescNorm<1e-3)
                resecueEvents(ind_sim,ind_p) = resecueEvents(ind_sim,ind_p) + 1;
            end
            
            Ttip = Ttip + T(ind_T);
            ind_T = ind_T + 1;
        end
        disp([ind_sim,K,ind_p])
    end
end
%%
figure
boxplot([resecueEvents(:,1),resecueEvents(:,2),resecueEvents(:,3)],'whisker',10,'Color','k')
xlim([0.4 3.8])
%% plotting
% % load('freq_climate_v2_rand.mat')
% 
% figure;
% subplot(2,2,1)
% hold on
% 
% plot(R_bef,'r.','MarkerSize',7)
% plot([0 K],mean(R_bef)*[1 1],'-r','LineWidth',2)
% plot(R_aft,'b.','MarkerSize',7)
% plot([0 K],mean(R_aft)*[1 1],'-b','LineWidth',2)
% xticks([0 200 300 500 700 900])
% xticklabels([0 100 300 500])
% yticks([2 2.5 3 3.5 4])
% yticklabels([2 2.5 3 3.5 4])
% axis([0 K 1.8 4])
% hA = gca;
% set(gca,'XMinorTick','on','YMinorTick','on');
% hA.XAxis.MinorTickValues = 0:20:1000;
% hA.YAxis.MinorTickValues = 1.5:0.05:5;
% set(gca,'FontSize',15)
% box on
% xlabel('Simulation index')
% 
% subplot(2,2,2)
% hold on
% plot(var_bef(1,:),mult*var_bef(2,:),'k.','Markersize',7)
% axis([0 15 0 35])
% 
% xticks([0 3 6 9 12 15])
% xticklabels([0 3 6 9 12])
% yticks([0 5 15 25 35])
% yticklabels([0 5 15 25 ])
% hA = gca;
% set(gca,'XMinorTick','on','YMinorTick','on');
% hA.XAxis.MinorTickValues = 0:0.5:15;
% hA.YAxis.MinorTickValues = 0:1:35;
% xlabel('$N$')
% ylabel('\tilde{P}','Rotation',0)
% set(gca,'FontSize',15)
% box on
% 
% subplot(2,2,3)
% hold on
% histogram(t_tip,50,'Normalization','probability','FaceColor',[.8 .8 .8],'LineWidth',1.5)
% 
% xticks([0 500 1000 1500 2000])
% xticklabels([0 0.5 1 1.5])
% yticks([0 0.05 .1])
% yticklabels([0 0.05 .1])
% axis([0 2000 0 .1])
% hA = gca;
% set(gca,'XMinorTick','on','YMinorTick','on');
% hA.XAxis.MinorTickValues = 0:100:2000;
% hA.YAxis.MinorTickValues = 0:0.005:0.11;
% xlabel('Time until tipping')
% set(gca,'FontSize',15)
% box on
% 
% subplot(2,2,4)
% hold on
% histogram(phi_bef,50,'Normalization','probability','FaceColor',[.8 .8 .8],'LineWidth',1.5)
% 
% xticks([-pi -pi/2 0 pi/2 pi])
% xticklabels({'-\pi' '-\pi/2' '0' '\pi/2' ' '})
% yticks([0 .1 .2 .3 .4])
% yticklabels([0 .1 .2 .3 .4])
% axis([-pi pi 0 .4])
% hA = gca;
% set(gca,'XMinorTick','on','YMinorTick','on');
% hA.XAxis.MinorTickValues = -pi:0.05*pi:pi;
% hA.YAxis.MinorTickValues = 0:0.01:0.5;
% xlabel('$\varphi$')
% set(gca,'FontSize',15)
% box on
%%

% ode of the May Model
function [dvar] = Mayodefun(var,R)

%parameters
% R       = RRR(t);
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

% myEvent to stop the integration when it tips
function [value, isterminal, direction] = myEvent(t,y)
TOL = 1e-3;
value      = norm(y)<TOL;
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions?
end