% diterbution of tipping time on varying \Delta_r and rho linearly

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmid = 2.075;
ampMin = 0.4;
ampMax = 0.5;

Pmin  = 0.5;
Pmax  = 0.1;

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

Tend    = 5000;
RR      = 1;
PnBin   = @(tt)( ((Pmax   - Pmin  )/Tend)*tt + Pmin);
amp     = @(tt)( ((ampMax - ampMin)/Tend)*tt + ampMin );
Rstar   = @(tt) Rmid + amp(tt);
Rend    = @(tt) Rmid - amp(tt);

K = 10000; %10000

t_tip   = NaN(1,K);

ind_sim = 0;
while (ind_sim < K)
    vars_cl = {'T','ind_T','t','var','R'};
    clear(vars_cl{:})
    
    ind_sim = ind_sim + 1;
    
    T(1) = nbininv(rand,RR,PnBin(0));
    R(1) = Rend(0) + rand*(Rstar(0) - Rend(0));
    initcond(:,1)  = [3;0.002];
    while (T(1) == 0)
        T(1) = nbininv(rand,RR,PnBin(0));
    end
    
%     odefun   = @(t,var)RModefun(var,R(1));
%     tspan    = [0,T(1)];
%     [t,var]  = ode45(odefun,tspan,initcond(:,1),opts);
%     initcond(:,1) = var(end,:);
    
    ind_T = 2;
    
    NORM = inf;
    Ttip = 0;
    while and(Ttip<Tend , NORM>1e-3)
        T(ind_T) = nbininv(rand,RR,PnBin(T(ind_T - 1)));
        while (T(ind_T) == 0)
            T(ind_T) = nbininv(rand,RR,PnBin(T(ind_T - 1)));
        end
        R(ind_T) = Rend(T(ind_T-1)) + rand*(Rstar(T(ind_T-1)) - Rend(T(ind_T-1)));
        
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
    elseif and(Ttip<Tend ,Ttip~=0)
        t_tip(ind_sim) = Ttip;
    end
    
    disp([ind_sim,K])
end
%%
tTipLength = 50;
plotResilotion = round(Tend/(tTipLength));
tTipFreq = zeros(1,tTipLength+1);

for ind_freq = 2:length(tTipFreq)
    for ind_sim = 1:K
        if ( t_tip(ind_sim)> plotResilotion*(ind_freq-2) && ...
                t_tip(ind_sim)< plotResilotion*(ind_freq-1))
            tTipFreq(ind_freq) = tTipFreq(ind_freq) + 1;
        end
        
    end
    tTipFreq(ind_freq) = tTipFreq(ind_freq)/K;
end

timeSpan = 0:plotResilotion:Tend;

figure
plot(...
    timeSpan,tTipFreq,'b-',...
    'LineWidth',3)
%%
% tTip_0p1_sorted_tmp = sort(t_tip(:,1));
% tTip_0p3_sorted_tmp = sort(t_tip(:,2));
% tTip_0p5_sorted_tmp = sort(t_tip(:,3));
% 
% nanInd_0p1 = find(isnan(tTip_0p1_sorted_tmp),1);
% nanInd_0p3 = find(isnan(tTip_0p3_sorted_tmp),1);
% nanInd_0p5 = find(isnan(tTip_0p5_sorted_tmp),1);
% 
% MIN = min([nanInd_0p1,nanInd_0p3,nanInd_0p5]);
% 
% tTip_0p1_sorted = tTip_0p1_sorted_tmp(1:MIN-1);
% tTip_0p3_sorted = tTip_0p3_sorted_tmp(1:MIN-1);
% tTip_0p5_sorted = tTip_0p5_sorted_tmp(1:MIN-1);
% 
% 
% boxplot([tTip_0p1_sorted,tTip_0p3_sorted,tTip_0p5_sorted],'whisker',1000)

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