% climate_tip_May
% simulate one time series that tips to exctinction. 
% q is going to be fixes through out, and R is going to be decreased up
% to vary between Rmin and Rmax and.

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rmax   = 3.3;       % /yr.              Prey intrinsic growth rate
Rmin    = 2;
Rmean   = (Rmax+Rmin)/2;
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
%
rng(34)% rng(316) %rng(136)
%%
Tend  = 100;
RR    = 1;
PnBin = .2;
% [swich,cli] = Climate(Tend,RR,PnBin);
[time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
climate = @(tt)interp1(time,cli1,tt);
R       = @(tt)Rmean + climate(tt)*(Rmax-Rmin)/2;

figure
subplot(6,12,[1:7,13:19,25:31])
hold on
plot(time,Rmean + cli1.*((Rmax-Rmin)/2),':k','LineWidth',1)
for ind = 1:length(swich)
    X = linspace(sum(swich(1:ind-1)),sum(swich(1:ind)),swich(ind));
    Y = Rmean + cli(sum(swich(1:ind-1))+1:sum(swich(1:ind)))*(Rmax-Rmin)/2;
    if Y(1) >= Rmean
        plot(X,Y,'r-','LineWidth',2)
    else
        plot(X,Y,'b-','LineWidth',2)
    end
    if swich(ind) == 1
        X1 = [X-1,X];
        Y1 = [Y Y];
        if Y(1) >= Rmean
            plot(X1,Y1,'r-','LineWidth',2)
        else
            plot(X1,Y1,'b-','LineWidth',2)
        end
    end
end
box on
axis([0 Tend 1.5 4])
ylabel('$R(t)$','Rotation',0)


%% timeseries

mult     = 1000;
initcond = [11,0.004];
odefun  = @(t,var)Mayodefun_cli(t,var,R);
[t,var] = ode45(odefun,[0 Tend],initcond,opts);
subplot(6,12,[37:43,49:55,61:67])
hold on
plot(...
    t,var(:,1),'-k','LineWidth',2)
plot(...
    t,mult*var(:,2),'-','LineWidth',2,'Color',[.65 .65 .65])
box on
axis([0 Tend -2 38])
hA=gca;
% set(gca,'XMinorTick','on','YMinorTick','on')
% hA.XAxis.MinorTickValues = 0:5:100;
% xticks([0 20 40 60 80 100])
% xticklabels([0 20 40 60 80])
% hA.YAxis.MinorTickValues = -2:1:45;
% yticks([0 10 20 30])
% set(gca,'FontSize',15)
% xlabel('$t$')
% ylabel('$P,N$','Rotation',0

%% phase

subplot(6,12,[8:12,20:24,32:36,44:48,56:60,68:72])
%%

hold on
% 1) the phase 
for ind = 1:length(swich)
    tspan = [sum(swich(1:ind-1)), sum(swich(1:ind))];
    rt =@(t)(Rmean + cli(sum(swich(1:ind-1))+1)*(Rmax-Rmin)/2);
    odefun  = @(t,var)Mayodefun_cli(t,var,rt);
    [t,var] = ode45(odefun,tspan,initcond,opts);
    if rt(1) >= Rmean
        plot(var(:,1),mult*var(:,2),'r-','LineWidth',2)
    else
        plot(var(:,1),mult*var(:,2),'b-','LineWidth',2)
    end
    initcond = var(end,:);
    vars_cl = {'t','var','rt'};
    clear(vars_cl{:})

end
%% plotting the limit cycles at r1 and r2 with thier threshold
mult     = 1000;
T        = 100;

% periodic orbit at r1
initcond = [8   0.01];
Tper1     = 9;%10;%8.1480;
R1       = 3.226;


ivpfun1   = @(t,var)Mayodefun_cli(t,var,@(t)R1);
perfun1   = @(t,var,Tper)Tper*Mayodefun_cli(t,var,@(t)R1);

[~,tempvar]  = ode45(ivpfun1,[0 T],initcond);
initper  = tempvar(end,:); %the top point is [6.9261, 0.0339];
ivpsol1 = ode45(@(t,var)perfun1(t,var,Tper1),[0 1],initper, opts);

ss = linspace(0,1,200);
tempinit1 = @(s)deval(ivpsol1,s);
solinit1 =bvpinit(ss,tempinit1,Tper1);

BC=@(WL,WR,Tper1)(...
    [WL(1)-WR(1);...
    WL(2)-WR(2);...
    WL(2)-initper(2);... %point phase condition
    ]);

per_r1 = bvp5c(perfun1,BC,solinit1);
var_r1 = deval(per_r1,ss);

% threshold for r1 
[N_eq] = May_eq(R1,q);
e5(1) =  N_eq(3);
e5(2) = (N_eq(3) + eps)/q;
G   = @(var)Mayodefun_cli(0,var,@(t)R1);
JJ = MyJacobian(G,e5);
[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-4*pert';
maninit = e5 + pert;
[t_man_1,varman1]   = ode45(ivpfun1,[10,0],maninit,opts);
Pthe_r1             = @(NN)interp1(varman1(:,1),mult*varman1(:,2),NN);

% periodic orbit at r2 
initcond = [8   0.01];
Tper2     = 11;
R2       = 2.068;

ivpfun2   = @(t,var)Mayodefun_cli(t,var,@(t)R2);
perfun2   = @(t,var,Tper)Tper*Mayodefun_cli(t,var,@(t)R2);

[~,tempvar2]  = ode45(ivpfun2,[0 T],initcond);
initper2  = tempvar2(end,:); %the top point is [6.9261, 0.0339];
ivpsol2 = ode45(@(t,var)perfun2(t,var,Tper2),[0 1],initper2, opts);

ss = linspace(0,1,200);
tempinit2 = @(s)deval(ivpsol2,s);
solinit2 =bvpinit(ss,tempinit2,Tper2);

BC=@(WL,WR,Tper2)(...
    [WL(1)-WR(1);...
    WL(2)-WR(2);...
    WL(2)-initper2(2);... %point phase condition
    ]);

per_r2 = bvp5c(perfun2,BC,solinit2);
var_r2 = deval(ivpsol2,ss);

% threshold for r2
[N_eq] = May_eq(R2,q);
e5(1) =  N_eq(3);
e5(2) = (N_eq(3) + eps)/q;
G   = @(var)Mayodefun_cli(0,var,@(t)R2);
JJ = MyJacobian(G,e5);
[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-4*pert';
maninit = e5 + pert;
[t_man_2,varman2]   = ode45(ivpfun2,[10,0],maninit,opts);
Pthe_r2             = @(NN)interp1(varman2(:,1),mult*varman2(:,2),NN);

% shaded area
resl = 500;
Nscan = linspace(0,15,resl);
Pscan = linspace(0,40,resl);
shading_mat = NaN(resl,resl);
for ind_N = 1:resl
    for ind_P = 1:resl
        if Pscan(ind_P) > Pthe_r2(Nscan(ind_N)) && Pscan(ind_P) < Pthe_r1(Nscan(ind_N))
           shading_mat(ind_P,ind_N) = 1;
        end
    end
end

% plotting 
figure(10); hold on
map = [1 1 1; 0.5 0.5 0.5 ;0 0 0];
plot(...
    var_r1(1,:),mult*var_r1(2,:),...
    'LineWidth',3,'Color','r')
plot(...
    var_r2(1,:),mult*var_r2(2,:),...
    'LineWidth',3,'Color','b')
plot(varman1(:,1),mult*varman1(:,2),'r--','LineWidth',1)
plot(varman2(1:end,1),mult*varman2(1:end,2),'b--','LineWidth',1)
PCOLOR = pcolor(Nscan,Pscan,shading_mat);
PCOLOR.FaceAlpha = 0.3;
PCOLOR.LineStyle = 'none';
colormap(map)
axis([0 15 0 40])
%%
% changing climate
function [time,swich,climate,climate1] = Climate(Tend,RR,PnBin)
%climatswich generats a sequence of negative binomial random variables that
%add up to Tend, has propablility of success in singil tril PnBin and
% corresponding number of successes RR
swichtemp = nbininv(rand(1,2*Tend),RR,PnBin); % how many swiches
indend = 1;
while (sum(swichtemp(1:indend))<Tend)
    indend = indend+1;
end
swichtemp1 = swichtemp(1:indend); % number of swiches befor Tend
ind_sw = 1;
for ind_sw1=1:length(swichtemp1)% clearing the swich list from zeros
    if swichtemp1(ind_sw1)~=0
        swich(ind_sw) = swichtemp1(ind_sw1);
        ind_sw = ind_sw+1;
    end
end
climate1 = NaN(1,sum(swich)); %climat vector sum(swich) long
ind_cl = 1;
for ind_sw = 1:length(swich)
    climatampl = rand;
    if rand<=0.5 %mod(ind_sw,2)==1 % good years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = climatampl;
            ind_cl = ind_cl + 1;
        end
    else  % bad years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = -climatampl;
            ind_cl = ind_cl + 1;
        end
    end
end
time = 0:0.001:sum(swich);
climate = NaN(size(time));
for ind_cl=1:length(time)-1
    climate(ind_cl) = climate1(floor(time(ind_cl)+1));
end
end


% ode of the May Model
function [dvar] = Mayodefun_cli(t,var,RRR)

%parameters
R       = RRR(t);
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



