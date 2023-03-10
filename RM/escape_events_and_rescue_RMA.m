% escape events and rescue events before t = 550

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

eps = 1e-5;
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

Tend      = 550;
RR        = 1;
PP = linspace(0.001,0.99,20);

filename = 'test.mat';
m = matfile(filename, 'Writable', true);

K = 1000;

% t_tip   = NaN(K,length(PP));
escape_number_sim = zeros(1,K);
rescue_number_sim = zeros(1,K);
escape_number = zeros(1,length(PP));
rescue_number = zeros(1,length(PP));
for ind_p = 1:length(PP)
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

        NORM_tipping = inf;
        NORM_rescue  = inf;
        while (NORM_tipping > eps && tspan(end)<Tend)
            T(ind_T) = nbininv(rand,RR,PnBin);
            while (T(ind_T) == 0)
                T(ind_T) = nbininv(rand,RR,PnBin);
            end
            R(ind_T) = Rend + rand*(Rstar - Rend);

            % Now I need to test whether the systems escapes
            odefun   = @(t,var)RModefun(var,R(ind_T));
            tspan    = [sum(T(1:ind_T - 1)),sum(T(1:ind_T))];
            [t,var]  = ode45(odefun,tspan,initcond(:,ind_T-1),opts);
            [~,var_esc]  = ode45(odefun,[0 100],initcond(:,ind_T-1),opts);
            initcond(:,ind_T) = var(end,:);
            NORM_tipping = norm(initcond(:,ind_T));
            NORM_rescue = norm(var_esc(end,:));
            ind_T = ind_T + 1;
            if (NORM_rescue < eps && NORM_tipping > eps)
                % there is escape  and resecue event
                escape_number_sim(ind_sim) =  1;
                rescue_number_sim(ind_sim) =  1;
            elseif (NORM_rescue < eps && NORM_tipping < eps)
                % there is escape only (P-tipping)
                escape_number_sim(ind_sim) = 1;
            end
        end
        escape_number(ind_p) = sum(escape_number_sim);
        rescue_number(ind_p) = sum(rescue_number_sim);
    end
    disp([ind_p,escape_number(ind_p),rescue_number(ind_p)])
    m.escape_number = escape_number;
    m.rescue_number = rescue_number;
    m.ind_p = ind_p;
end
%%
% figure
% histogram(escape_number(1:520))
% hold on
% histogram(rescue_number(1:520))

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