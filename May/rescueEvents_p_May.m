% to compute how many rescue events befortipping in the May model
warning off
% clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB\10_Phase_Sensitivity\P-tippingPaper\may_model')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar   = 3.3;       % /yr.              Prey intrinsic growth rate
Rend    = 2;
Rmean   = (Rstar+Rend)/2;
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
pp     = [.1 .3 .5];
mult      = 1000;

K = 1000;

resecueEvents = zeros(K,3);

for ind_p = 1:3
    PnBin = pp(ind_p);
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

        while (sum(T)<Tend)
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
            
            ind_T = ind_T + 1;
        end
        disp([ind_p,ind_sim,K])
    end
end
%%
figure
boxplot([resecueEvents(:,1),resecueEvents(:,2),resecueEvents(:,3)],'whisker',10)
xlim([0.5 3.8])

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