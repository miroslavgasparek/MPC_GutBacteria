%% 27 May 2019 Miroslav Gasparek
% The script for running of the optimal nutrient dosing for the
% predator-prey model of the Loss of Microbial Diversity (LOMD)
% 
% The model equations and the parameters are adopted from:
% Astrom, K. J., & Murray, R. M. (2016). 
% Feedback Systems: An Introduction for Scientists and Engineers.
% https://doi.org/10.1086/596297
% 
% System equations are as follows:
% 
% dB/dt = (r+u)*B*(1-B/k) - a*B*P/(c+B), (prey dynamics)
% dP/dt = b*a*B*P/(c+B) - d*P,           (predator dynamics) 
% 
% Where:
% - B is the prey bacteria population which overgrew [g]
% - P is the predator bacteria population that we aim to increase [g] 
% 
% - u is the input representing the modulation of the food source to the
% bacteria [1/day] 
% 
% - r is the growth rate of prey bacteria [1/day]
% - k is carrying capacity of prey bacteria, what the population would be
% if the growth was not limited [g]
% - a is the interaction term describing how fast the predator bacteria can
% consume the prey bacteria [1/day]
% c - is the control of prey consumption rate for low prey bacterial
% population [g]
% - b is the growth coefficient of predator bacteria [dimensionless]
% - d is the mortality rate of predator bacteria [1/day]



clear; clc; close all;

% Load the model parameters, these are stored in genGutParameters.mat
% These parameters will be automatically used in the subsequent
% calculations and simulations
load gut_parameters;
%% Define the parameters of the system
a = 3.2; % [1/day]
b = 0.6; % [dimensionless]
c = 50; % [g]
d = 0.56; % [1/day]
k = 125;  % [g]
r = 1.6; % [1/day]

Ts=1/4; % Sampling time in days
Tf=1; % duration of prediction horizon in days
T = 20; % simulation duration in days

% Number of prediction steps 
N=ceil(Tf/Ts);

% Generate the linearized system matrices and compute the equilibrium
% points
[A,B,C,~, Beq, Peq] = genGutODE(a, b, c, d, k, r, Ts);

%% Declare initial conditions and target conditions
% Select the target mass of bacterial species
BTarget = Beq; % g, prey mass
PTarget = Peq; % g, predator mass

% Select the initial amounts of the bacterial species
% Prey bacterial mass
Binit = 100; % g
% Predator bacterial mass
Pinit = 5; % g
%% Declare penalty matrices:
% Q: the penalty on the state
% R: the penalty on the inputs
% P: the penalty on the terminal state
Q=1*eye(2);
P=10*eye(2);
R=0.1*eye(1);

%% Declare contraints
% The concentrations must be slightly above zero, constraints are 
% relative to the equilibrium points
cl=[-0.95*Beq; -0.95*Peq];
% Select the maximum masses of the species
ch=[3*Beq; 2*Peq];

% We cannot feed bacteria at "negative" rate, so minimum growth rate is 10%
% of the normal growth rate
ul= - 0.9*r;

% Input cannot be more than triple the normal nutrient rate
uh=2*r;

% Assume that the variables are independent, hence D is the identity matrix
D=eye(2);

%% Compute stage constraint matrices and vector
[Dt,Et,bt] = genStageConstraints(A,B,D,cl,ch,ul,uh);

%% Compute trajectory constraints matrices and vector
[DD,EE,bb] = genTrajectoryConstraints(Dt,Et,bt,N);

%% Compute QP constraint matrices
[Gamma,Phi] = genPrediction(A,B,N);
[F,J,L] = genConstraintMatrices(DD,EE,Gamma,Phi,N);

%% Compute QP cost matrices
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
% Prepare cost and constraint matrices for mpcqpsolver
% Calculating the inverse of the lower triangular H. see doc mpcqpsolver.
H = chol(H,'lower');
H = (H'\eye(size(H)))';

%% Run a linearized simulation of the MPC 
xTarget=[BTarget-Beq PTarget-Peq]';% target equilibrium state
x0 = [Binit-Beq Pinit-Peq]'; % starting weights of bacteria relative to the equilibrium
iA = false(size(bb));
t = 0:Ts:T;
x = [x0, zeros(2,length(t)-1)];
for t_step=1:length(t)-1
    % Solve for the optimal input at each step
    [u,~,iA] = genMPController(H,G,F,bb,J,L,x(:,t_step),xTarget,1,iA);
    % Perform the transition to the next state
    x(:,t_step+1)=A*x(:,t_step)+B*u;
    % Store the input to the vector
    u_vec(:,t_step) = u;
end

%% Run the nonlienar system simulation without the input 
sol = ode45(@(t,y) gut_bacteria_ode(t,y,a,b,c,d,k,r), [0 T], [Binit Pinit]);
sol_vals = deval(sol,t);

Bsim = sol_vals(1,:);
Psim = sol_vals(2,:);

%% Plot results
% Define some plotting variables
B = Beq + x(1,:); % Prey bacteria
P = Peq + x(2,:); % Predator bacteria 
U = [zeros(1), u_vec]; % Input

% Plot the input constraints
Umin = ul*ones(1,length(t));
Umax = uh*ones(1,length(t));

% Plot the target equilibrium concentrations
BTarget_plot = BTarget*ones(1,length(t));
PTarget_plot = PTarget*ones(1,length(t));

% Make the actual plot
figure(1);
sgtitle('Gut microbiome restoration using MPC for optimal nutrients dosing',...
        'fontsize',20,'interpreter','latex')

subplot(3,1,1)
hold on
plot(t,BTarget_plot,'k--','Linewidth',2)
plot(t, B,'Linewidth',2)
title('Time evolution of prey bacteria mass','fontsize',15,'interpreter','latex')
xlabel('Time [Days]','fontsize',12,'interpreter','latex')
ylabel('Bacteria Mass [g]','fontsize',12,'interpreter','latex')
ax1 = gca;
ax1.YLim  = [0 110];

subplot(3,1,2)
hold on
plot(t,PTarget_plot,'k--','Linewidth',2)
plot(t, P,'Linewidth',2)
title('Time evolution of predator bacteria mass','fontsize',15,'interpreter','latex')
xlabel('Time [Days]','fontsize',12,'interpreter','latex')
ylabel('Bacteria Mass [g]','fontsize',12,'interpreter','latex')
ax2 = gca;
ax2.YLim  = [0 110];

subplot(3,1,3)
hold on
pUmin = plot(t,Umin,'--','Linewidth',2,'Color','Red');
pUmax = plot(t,Umax,'--','Linewidth',2,'Color','Red');
pU = plot(t,U,'b.','markersize',15,'Linewidth',2);
title('Optimal nutrient feeding rate','fontsize',15,'interpreter','latex')
xlabel('Time [Days]','fontsize',12,'interpreter','latex')
ylabel('Rate [1/Day]','fontsize',12,'interpreter','latex')
ax3 = gca;
ax3.YLim  = [-2 5];
ax3.Children(1).Color = [0 0.4470 0.7410];

fig = gcf;
fig.Position = [488.0000  107.4000  861.0000  654.6000];


% Compare the plot with uncontrolled system
figure(2);
sgtitle('Gut microbiome evolution without nutrient dosing',...
        'fontsize',20,'interpreter','latex')

subplot(2,1,1)
hold on
plot(t,BTarget_plot,'k--','Linewidth',2)
plot(t, Bsim,'Linewidth',2)
title('Time evolution of prey bacteria mass','fontsize',15,'interpreter','latex')
xlabel('Time [Days]','fontsize',12,'interpreter','latex')
ylabel('Bacteria Mass [g]','fontsize',12,'interpreter','latex')
ax1 = gca;
ax1.YLim  = [-10 110];

subplot(2,1,2)
hold on
plot(t,PTarget_plot,'k--','Linewidth',2)
plot(t, Psim,'Linewidth',2)
title('Time evolution of predator bacteria mass','fontsize',15,'interpreter','latex')
xlabel('Time [Days]','fontsize',12,'interpreter','latex')
ylabel('Bacteria Mass [g]','fontsize',12,'interpreter','latex')
ax2 = gca;
ax2.YLim  = [0 110];

fig = gcf;
fig.Position = [488.0000 107.4000  861.0000  654.6000];