% frequency of tipping in terms of phases \varphi (picewize auotnomous integration)

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar  = 2.55;
Rend   = 1.6;
del    = 2.2;       %/yr.            Predator death rate in absent of prey
C      = 0.19;      %/yr.            Nonlinear death rate of prey
gamma  = 0.004;     % pred/prey.      Prey-predator conversion rate
beta   = 1.5;       % prey/ha        Predator half-saturating constant
alpha  = 800;       % prey/(pred.yr) Predator saturating kill rate
mu     = 0.03;      % NA             Allee parameter
nu     = 0.003;     % NA             Allee parameter

opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(4)
%% Simulations

Tend      = 600;
RR        = 1;  %avreage length of Type-L/H period
PP    = 0.001:0.05:0.999;


K = 1000;

t_tip   = NaN(K,length(PP));
observable = zeros(size(PP));
for ind_p = 1:length(PP)
    PnBin = PP(ind_p);
    ind_sim = 0;
    tippingFreq = 0;
    while (ind_sim < K)
        
        vars_cl = {'T','ind_T','t','var','R'};
        clear(vars_cl{:})
        
        ind_sim = ind_sim + 1;
        
        T(1) = nbininv(rand,RR,PnBin);
        R(1) = Rend + rand*(Rstar - Rend);
        initcond(:,1)  = [3;0.002];
        while (T(1) == 0)
            T(1) = nbininv(rand,RR,PnBin);
        end
        
        odefun   = @(t,var)RModefun(var,R(1));
        tspan    = [0,T(1)];
        [t,var]  = ode45(odefun,tspan,initcond(:,1),opts);
        initcond(:,1) = var(end,:);
        
        ind_T = 2;
        
        NORM = inf;
        Ttip = 0;
        
        while and(Ttip<Tend , NORM>1e-5)
            T(ind_T) = nbininv(rand,RR,PnBin);
            while (T(ind_T) == 0)
                T(ind_T) = nbininv(rand,RR,PnBin);
            end
            R(ind_T) = Rend + rand*(Rstar - Rend);
            
            odefun   = @(t,var)RModefun(var,R(ind_T));
            tspan    = [sum(T(1:ind_T - 1)),sum(T(1:ind_T))];
            [t,var]  = ode45(odefun,tspan,initcond(:,ind_T-1),opts);
            initcond(:,ind_T) = var(end,:);
            NORM = norm(initcond(:,ind_T));
            Ttip = Ttip + T(ind_T);
            ind_T = ind_T + 1;
        end
        
        if Ttip <= 0
            ind_sim = ind_sim -1;
        elseif and(Ttip<Tend-50,Ttip~=0)
            t_tip(ind_sim,ind_p) = Ttip;
            tippingFreq = tippingFreq + 1;
        end
 
    end
    observable(ind_p) = tippingFreq/K; %nanmean(t_tip(:,ind_p));
    disp([ind_p,PP(ind_p),observable(ind_p)])
end
%%
figure(1);hold on
plot(...
    PP,observable,'-','LineWidth',2,'MarkerSize',10)
box on
set(gca,'FontSize',20)
axis([0 1 0 0.7]);
xticks([0 0.2 0.4 0.6 0.8 1])
xticklabels([0 0.2 0.4 0.6 0.8])

yticks([0 0.2 0.4 0.6 0.8])
yticklabels([0 0.2 0.4])

xlabel('$p$','Position',[1 -0.001 0])
ylabel('$\textrm{Pr}(t_{ex} \leq t)$','Rotation',90,'Position',[-0.03 0.8 0])
%%
tTipFreq = zeros(3,Tend);
for ind_p = 1:3
    for ind_freq = 1:length(tTipFreq)
        for ind_sim = 1:K
            if t_tip(ind_sim,ind_p)== ind_freq
                tTipFreq(ind_p,ind_freq) = tTipFreq(ind_p,ind_freq) + 1;
            end
        end
        tTipFreq(ind_p,ind_freq) = tTipFreq(ind_p,ind_freq)/K;
    end
end


%%
%the ode fuction of RM model
function [dvar] = RModefun(var,R)
%parameters

delta      = 2.2;
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter


N = var(1);
P = var(2);

dN = R * N *(1-((C/R)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = gamma * ((alpha*P*N)/(beta + N)) - (delta*P);

dvar = [dN;dP];
end

% myEvent to stop the integration when it tips
function [value, isterminal, direction] = myEvent(t,y)
TOL = 1e-5;
value      = norm(y)<TOL;
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions?
end