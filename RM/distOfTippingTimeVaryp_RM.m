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

Tend      = 5000;
RR        = 1;  %avreage length of Type-L/H period
PP    = [0.1,0.3,0.5];

K = 10000;

t_tip   = NaN(K,3);

for ind_p = 1:3
    PnBin = PP(ind_p);
    ind_sim = 0;
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
        
        while and(Ttip<Tend , NORM>1e-3)
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
        elseif and(Ttip<Tend,Ttip~=0)
            t_tip(ind_sim,ind_p) = Ttip;
        end
 
        disp([ind_p,ind_sim,K])
    end
    nontippingevents(ind_p) = sum(isnan(t_tip(:,ind_p)));
end

% 100.*nontippingevents./(10000 + nontippingevents)
% 11.5044    1.1174    0.9411
% 
%%
tTipLength = 50;
plotResilotion = round(Tend/(tTipLength));
tTipFreq = zeros(3,tTipLength+1);
for ind_p = 1:3
    for ind_freq = 2:length(tTipFreq)
        for ind_sim = 1:K
            if ( t_tip(ind_sim,ind_p)> plotResilotion*(ind_freq-2) && ...
                    t_tip(ind_sim,ind_p)< plotResilotion*(ind_freq-1))
                tTipFreq(ind_p,ind_freq) = tTipFreq(ind_p,ind_freq) + 1;
            end

        end
        tTipFreq(ind_p,ind_freq) = tTipFreq(ind_p,ind_freq)/K;
    end
end
timeSpan = 0:plotResilotion:Tend;
% tTipFreq_0p45 = interp1(linspace(1,400,400),tTipFreq(1,:),timeSpan);
% tTipFreq_0p5 = interp1(linspace(0,400,400),tTipFreq(2,:),timeSpan);
% tTipFreq_0p55 = interp1(linspace(0,400,400),tTipFreq(3,:),timeSpan);
figure
plot(...
    timeSpan,tTipFreq(1,:),'b-',...
    timeSpan,tTipFreq(2,:),'r-',...
    timeSpan,tTipFreq(3,:),'g-',...
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
% figure
% boxplot([tTip_0p1_sorted,tTip_0p3_sorted,tTip_0p5_sorted],'whisker',1000)

%%
% tTip_0p1_sorted_tmp = sort(t_tip(:,1));
% tTip_0p3_sorted_tmp = sort(t_tip(:,2));
% tTip_0p5_sorted_tmp = sort(t_tip(:,3));
% 
% nanInd_0p1 = find(isnan(tTip_0p1_sorted_tmp),1);
% nanInd_0p3 = find(isnan(tTip_0p3_sorted_tmp),1);
% nanInd_0p5 = find(isnan(tTip_0p5_sorted_tmp),1);
% 
% tTip_0p1_sorted = tTip_0p1_sorted_tmp(1:nanInd_0p1-1);
% tTip_0p3_sorted = tTip_0p3_sorted_tmp(1:nanInd_0p3-1);
% tTip_0p5_sorted = tTip_0p5_sorted_tmp(1:nanInd_0p5-1);
% 
% meanInd_0p1 = find(abs(mean(tTip_0p1_sorted)-tTip_0p1_sorted)<2,1);
% meanInd_0p3 = find(abs(mean(tTip_0p3_sorted)-tTip_0p3_sorted)<2,1);
% meanInd_0p5 = find(abs(mean(tTip_0p5_sorted)-tTip_0p5_sorted)<2,1);
% 
% uppInd_0p1 = meanInd_0p1 + 2000;
% uppInd_0p3 = meanInd_0p3 + 2000;
% uppInd_0p5 = meanInd_0p5 + 2000;
% 
% lowInd_0p1 = meanInd_0p1 - 2000;
% lowInd_0p3 = meanInd_0p3 - 2000;
% lowInd_0p5 = meanInd_0p5 - 2000;
% 
% eps = 0.03;
% figure; hold on
% % mean
% plot(...
%     [0.1-eps 0.1+eps],mean(tTip_0p1_sorted)*[1,1],'r',...
%     ...
%     [0.3-eps 0.3+eps],mean(tTip_0p3_sorted)*[1,1],'r',...
%     ...
%     [0.5-eps 0.5+eps],mean(tTip_0p5_sorted)*[1,1],'r',...
%     ...
%     'LineWidth',2)
% 
% % max
% plot(...
%     [0.1-eps 0.1+eps],max(tTip_0p1_sorted)*[1,1],'k',...
%     ...
%     [0.3-eps 0.3+eps],max(tTip_0p3_sorted)*[1,1],'k',...
%     ...
%     [0.5-eps 0.5+eps],max(tTip_0p5_sorted)*[1,1],'k',...
%     ...
%     'LineWidth',2)
% 
% %min
% plot(...
%     [0.1-eps 0.1+eps],min(tTip_0p1_sorted)*[1,1],'k',...
%     ...
%     [0.3-eps 0.3+eps],min(tTip_0p3_sorted)*[1,1],'k',...
%     ...
%     [0.5-eps 0.5+eps],min(tTip_0p5_sorted)*[1,1],'k',...
%     ...
%     'LineWidth',2)
% 
% % vertical
% plot(...
%     [0.1 0.1],[min(tTip_0p1_sorted),max(tTip_0p1_sorted)],'--k',...
%     ...
%     [0.3 0.3],[min(tTip_0p3_sorted), max(tTip_0p3_sorted)],'--k',...
%     ...
%     [0.5 0.5],[min(tTip_0p5_sorted),max(tTip_0p5_sorted)],'--k',...
%     ...
%     'LineWidth',1)
% 
% 
% % upper and lower bounds
% plot(...
%     [0.1-eps 0.1+eps],tTip_0p1_sorted(uppInd_0p1)*[1,1],'b',...
%     [0.1-eps 0.1+eps],tTip_0p1_sorted(lowInd_0p1)*[1,1],'b',...
%     ...
%     [0.3-eps 0.3+eps],tTip_0p3_sorted(uppInd_0p3)*[1,1],'b',...
%     [0.3-eps 0.3+eps],tTip_0p3_sorted(lowInd_0p3)*[1,1],'b',...
%     ...
%     [0.5-eps 0.5+eps],tTip_0p5_sorted(uppInd_0p5)*[1,1],'b',...
%     [0.5-eps 0.5+eps],tTip_0p5_sorted(lowInd_0p5)*[1,1],'b',...
%     ...
%     'LineWidth',1)
% 
% % sides
% plot(...
%     [0.1-eps 0.1-eps],[tTip_0p1_sorted(uppInd_0p1), tTip_0p1_sorted(lowInd_0p1)],'b',...
%     [0.1+eps 0.1+eps],[tTip_0p1_sorted(uppInd_0p1), tTip_0p1_sorted(lowInd_0p1)],'b',...
%     ...
%     [0.3-eps 0.3-eps],[tTip_0p3_sorted(uppInd_0p3), tTip_0p3_sorted(lowInd_0p3)],'b',...
%     [0.3+eps 0.3+eps],[tTip_0p3_sorted(uppInd_0p3), tTip_0p3_sorted(lowInd_0p3)],'b',...
%     ...
%     [0.5-eps 0.5-eps],[tTip_0p5_sorted(uppInd_0p5), tTip_0p5_sorted(lowInd_0p5)],'b',...
%     [0.5+eps 0.5+eps],[tTip_0p5_sorted(uppInd_0p5), tTip_0p5_sorted(lowInd_0p5)],'b',...
%     ...
%     'LineWidth',1)




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