% climate_tip_May
% to illustrate resecue events by plotting two time series one tips and the
% other not

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
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
mult    = 1000;
opts = odeset('RelTol',1e-10,'AbsTol',1e-16);
%
% rng(34)% rng(316) %rng(136)
%%
initCond = [13.17 5.519/1000];
tspan1 = [0 3];
tspan2 = [3,4];
tspan3 = [4,9];

odefun1  = @(t,var)Mayodefun_cli(t,var,Rstar);
odefun2  = @(t,var)Mayodefun_cli(t,var,Rend);
odefun3  = @(t,var)Mayodefun_cli(t,var,2.8);

[~,var1] = ode45(odefun1,tspan1,initCond,opts);
[~,var2] = ode45(odefun2,tspan2,var1(end,:),opts);
[~,var3] = ode45(odefun3,tspan3,var2(end,:),opts);

figure; 
subplot(1,2,1)
hold on
plot(var1(:,1),mult*var1(:,2),'r-','LineWidth',2)
plot(var2(:,1),mult*var2(:,2),'b-','LineWidth',2)
plot(var3(:,1),mult*var3(:,2),'r-','LineWidth',2)

subplot(1,2,2)
hold on
plot(...
    tspan2, Rend*ones(1,2),'b-',...
    tspan3, 2.8*ones(1,2),'r-',...
    'LineWidth',2);
plot(...
    [tspan2(2),tspan2(2)], [Rend,2.8],':k',...
    'LineWidth',1);
axis([1e-5 20 0 40])

%% theta and shaded area! 
R2 = 2.8;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events', @myEvent);
odefun  = @(t,var)Mayodefun_cli(t,var,R2);
[N_eq] = May_eq(R2,q);
e5(1) =  N_eq(3);
e5(2) = (N_eq(3) + eps)/q;
G   = @(var)Mayodefun_cli(0,var,R2);
JJ = MyJacobian(G,e5);
[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-4*pert';
maninit1 = e5 + pert;
maninit2 = e5 - pert;
[~,varman1]   = ode45(odefun,[10,0],maninit1,opts);
[~,varman2]   = ode45(odefun,[10,0],maninit2,opts);
figure;
hold on; 
% linspace(0,varman(1,1),10)
% 1e-4*ones(10,1)
X = [linspace(0,varman2(1,1),10),varman2(:,1)',varman1(:,1)'];
Y = [1e-4*ones(10,1);mult*varman2(:,2);mult*varman(:,2)];
area(X,[Y,100*ones(size(Y))],'LineStyle','none')
axis([0 15 0 40])

%% ode of the May Model
function [dvar] = Mayodefun_cli(t,var,RRR)

%parameters
R       = RRR;
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


