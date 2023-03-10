% to create a color map of tipping before Tend for the May model
warning off
% clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmid = 2.075;

del    = 2.2;       %/yr.            Predator death rate in absent of prey
C      = 0.19;      %/yr.            Nonlinear death rate of prey
gamma  = 0.004;     % pred/prey.      Prey-predator conversion rate
beta   = 1.5;       % prey/ha        Predator half-saturating constant
alpha  = 800;       % prey/(pred.yr) Predator saturating kill rate
mu     = 0.03;      % NA             Allee parameter
nu     = 0.003;     % NA             Allee parameter

opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(10)   % random seeding

%%
tol = 1e-3;

Tend      = 550;
RR        = 1;

pRes      = 50;
drRes      = 50;

pScan     = linspace(0.1,0.9,pRes);
drScan     = linspace(0.3,0.5,drRes);

numberOfRuns = 1000;

colorMat = NaN(pRes,drRes);
%%
for ind_p = 45:pRes
    for ind_dr = 1:drRes
        
        p = pScan(ind_p); % geometric dist parameter
        dr = drScan(ind_dr); % the diffrence from Rmid;
        
        R1 = Rmid + dr;
        R2 = Rmid - dr;
        
        tippingNumber = 0;
        parfor ind_runs = 1: numberOfRuns
            T = 0;
            tend = 0;
            initCond  = [3;0.002];
            while (tend < Tend)
                while (T == 0)
                    T = nbininv(rand,RR,p);
                end
                tend = tend + T;
                R        = R2 + rand*(R1 - R2);
                odefun   = @(t,var)RModefun(var,R);
                [~,var]  = ode45(odefun,[0,T],initCond,opts);
                initCond = var(end,:);
            end
            if abs(var(end,1)) < tol
                tippingNumber = tippingNumber + 1;
            end
        end
        colorMat (ind_p,ind_dr) = tippingNumber/numberOfRuns;
        disp([ind_p,ind_dr])
    end
end
%% plotting 
cMap = interp1([0;1],[0.4 0.7 0.1; 0.6 0.1 0.2],linspace(0,1,256));
figure;
hold on
PCOLOR = pcolor(drScan,pScan,colorMat);
PCOLOR.FaceAlpha = .7;
PCOLOR.LineStyle = 'none';
colormap(hot)
caxis([0 0.5])
%% functions
% ode of the RM Model
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

function [value, isterminal, direction] = myEvent(t,y)
TOL = 1e-3;
value      = norm(y(1))<TOL;
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions?
end