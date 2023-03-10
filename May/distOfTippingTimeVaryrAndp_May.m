% frequency of tipping in terms of phases \varphi
warning off
% clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmid   = 2.65;
ampMin = 0.5; 
ampMax = 0.65;

Pmin = 0.1;
Pmax = 0.5;

q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(4)  % random seeding

%% Simulations

Tend      = 5000;
RR        = 1;  %avreage length of Type-L/H period
PnBin   = @(tt)( ((Pmax   - Pmin  )/Tend)*tt + Pmin);
amp     = @(tt)( ((ampMax - ampMin)/Tend)*tt + ampMin );
Rstar   = @(tt) Rmid + amp(tt);
Rend    = @(tt) Rmid - amp(tt);

K = 10000;

t_tip   = NaN(1,K);

ind_sim = 0;
while (ind_sim < K)
    vars_cl = {'T','ind_T','t','var','R'};
    clear(vars_cl{:})
    
    ind_sim = ind_sim + 1;
    
    T(1)           = nbininv(rand,RR,PnBin(0));
    R(1)           = NaN;
    initcond(:,1)  = [3;0.002];
    while (T(1) == 0)
        T(1) = nbininv(rand,RR,PnBin(0));
    end
    
%     odefun   = @(t,var)Mayodefun(var,R(1));
%     tspan    = [0,T(1)];
%     [t,var]  = ode45(odefun,tspan,initcond(:,1),opts);
%     initcond(:,1) = var(end,:);
    
    ind_T = 2;
    
    NORM = inf;
    Ttip = 0;
    while and(sum(T)<Tend , NORM>5e-2)
        T(ind_T) = nbininv(rand,RR,PnBin(T(ind_T - 1)));
        while (T(ind_T) == 0)
            T(ind_T) = nbininv(rand,RR,PnBin(T(ind_T - 1)));
        end
        R(ind_T) = Rend(T(ind_T - 1)) + rand*(Rstar(T(ind_T - 1)) - Rend(T(ind_T - 1)));
        
        odefun   = @(t,var)Mayodefun(var,R(ind_T));
        tspan    = [sum(T(1:ind_T - 1)),sum(T(1:ind_T))];
        [t,var]  = ode45(odefun,tspan,initcond(:,ind_T-1),opts);
        initcond(:,ind_T) = var(end,:);
        NORM = norm(initcond(:,ind_T));
        Ttip = Ttip + T(ind_T);
        ind_T = ind_T + 1;
    end
    
    if Ttip <= 0
        ind_sim = ind_sim -1;
    elseif (Ttip<Tend-1&&Ttip~=0)
        t_tip(1,ind_sim) = Ttip;
    end
    
    disp([ind_sim,K])
end

%% 
tTipLength = 250;
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
xlim([0,4000])
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
TOL = 5e-2;
value      = norm(y)<TOL;
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions?
end