% to illustrate resecue events by plotting two time series one tips and the
% other not

warning off
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

R1  = 1.6;
R2  = 2.4;
del    = 2.2;       %/yr.            Predator death rate in absent of prey
C      = 0.19;      %/yr.            Nonlinear death rate of prey
gamma  = 0.004;     % pred/prey.      Prey-predator conversion rate
beta   = 1.5;       % prey/ha        Predator half-saturating constant
alpha  = 800;       % prey/(pred.yr) Predator saturating kill rate
mu     = 0.03;      % NA             Allee parameter
nu     = 0.003;     % NA             Allee parameter

mult = 1000;

opts = odeset('RelTol',1e-16,'AbsTol',1e-16);
rng(4)  % random seeding 

%%
initCond = [4 21.5/1000];
TT = 1;
tspan1 = [0,TT];
tspan2 = [TT,9];


odefun1  = @(t,var)RModefun(var,R1);
odefun2  = @(t,var)RModefun(var,R2);

[~,var1] = ode45(odefun1,tspan1,initCond,opts);
[~,var2] = ode45(odefun2,tspan2,var1(end,:),opts);

figure; 

h1 = subplot(1,2,1);
hold on

plot(var1(:,1),mult*var1(:,2),'b-','LineWidth',2)
plot(var2(:,1),mult*var2(:,2),'r-','LineWidth',2)
set(h1,'XMinorTick','on','XScale','log','YMinorTick',...
    'on','YScale','log');

axis([1e-2 20 1e-2 40])

h2 = subplot(1,2,2);
hold on
%figure
hold on
plot(...
    tspan1, R1*ones(1,2),'b--',...
    tspan2, R2*ones(1,2),'r--',...
    'LineWidth',2);
plot(...
    [TT,TT], [R1,R2],':k',...
    'LineWidth',1);
axis([0 4 1.35 2.65])


%% theta and shaded area R2
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events', @myEvent);
odefun  = @(t,var)RModefun(var,R2);
e2      = [mu,0];
G       = @(var)RModefun(var,R2);
JJ      = MyJacobian(G,e2);
[eigvic,eigval] = eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-6*pert';
maninit = e2 + pert;
[~,varman]   = ode45(odefun,[10,0],maninit,opts);
figure;
hold on; 
X = [linspace(0,varman(1,1),10),varman(:,1)'];
Y = [1e-4*ones(10,1);mult*varman(:,2)];
area(X,[Y,100*ones(size(Y))],'LineStyle','none')
axis([0 15 0 100])
%% theta and shaded area R1

odefun  = @(t,var)RModefun(var,R1);
e2      = [mu,0];
G       = @(var)RModefun(var,R1);
JJ      = MyJacobian(G,e2);
[eigvic,eigval] = eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-6*pert';
maninit = e2 + pert;
[~,varman]   = ode45(odefun,[10,0],maninit,opts);
figure;
hold on; 
X = [linspace(0,varman(1,1),10),varman(:,1)'];
Y = [1e-4*ones(10,1);mult*varman(:,2)];
area(X,[Y,100*ones(size(Y))],'LineStyle','none')
axis([0 15 0 100])
%% ode of the May Model

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


