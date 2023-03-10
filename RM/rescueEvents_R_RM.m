% to compute how many rescue events befor tipping in the RM model
% (take so much time)

warning off
% clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');


Rmid = 2.075;
amp = [0.4250 0.475 0.525];

del    = 2.2;       %/yr.            Predator death rate in absent of prey
C      = 0.19;      %/yr.            Nonlinear death rate of prey
gamma  = 0.004;     % pred/prey.      Prey-predator conversion rate
beta   = 1.5;       % prey/ha        Predator half-saturating constant
alpha  = 800;       % prey/(pred.yr) Predator saturating kill rate
mu     = 0.03;      % NA             Allee parameter
nu     = 0.003;     % NA             Allee parameter


opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(4)  % random seeding 

%% Simulations

Tend      = 10000;
RR        = 1;  %avreage length of Type-L/H period
PnBin     = 0.3;
mult      = 1000;

K = 10000;

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
        R(1)           = Rend + rand*(Rstar - Rend);
        initcond(:,1)  = [3;0.002];
        
        ind_T = 2;
        
        NORM = inf;
        
        while (NORM>1e-4)
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
            
            [~,var_ersc]  = ode45(odefun,[tspan(1),tspan(1)+50],initcond(:,ind_T-1),opts);
            rescNorm = norm(var_ersc(end,:));
            
            if (NORM > 1e-3 && rescNorm<1e-3)
                resecueEvents(ind_sim,ind_p) = resecueEvents(ind_sim,ind_p) + 1;
            end

            ind_T = ind_T + 1;
        end
        disp([ind_p,ind_sim,K])
    end
end
%%

figure
boxplot(...
    resecueEvents,'whisker',10)
xlim([0.5 3.8])
xlable('$p$')

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