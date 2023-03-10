% to create a color map of tipping before Tend for the May model
warning off
% clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmid = 2.65;
amp = [0.6 0.65 0.7];



opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(10)   % random seeding
%%
tol = 1e-3;

Tend      = 550;
RR        = 1;

pRes      = 50;

pScan     = linspace(0.1,0.9,pRes);

numberOfRuns = 1000;

dr = 0.65; % the diffrence from Rmid;

colorMat = NaN(pRes,1);
%%
for ind_p = 1:pRes
       
        p = pScan(ind_p); % geometric dist parameter
        
        
        R1 = Rmid + dr;
        R2 = Rmid - dr;
        
        tippingNuber = 0;
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
                odefun   = @(t,var)Mayodefun(var,R);
                [~,var]  = ode45(odefun,[0,T],initCond,opts);
                initCond = var(end,:);
            end
            if abs(var(end,1)) < tol
                tippingNuber = tippingNuber + 1;
            end
        end
        colorMat (ind_p,1) = tippingNuber/numberOfRuns;
        disp([ind_p,50])
end
%% plotting 
% cMap = interp1([0;1],[0.4 0.7 0.1; 0.6 0.1 0.2],linspace(0,1,256));
% figure;
% hold on
% PCOLOR = pcolor(drScan,pScan,colorMat);
% PCOLOR.FaceAlpha = .7;
% PCOLOR.LineStyle = 'none';
% colormap(hot)
% caxis([0 1])
%% functions
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

function [value, isterminal, direction] = myEvent(t,y)
TOL = 1e-3;
value      = norm(y(1))<TOL;
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions?
end