%% 27 May 2019 Miroslav Gasparek
% Running the simulation of the linearized receptor/ligand binding system
clear; clc; close all;
load gut_parameters;
Ts=1/5;
Tf=1; % duration of prediction horizon in days
T = 15; % simulation duration in days

N=ceil(Tf/Ts);
[A,B,C,~, Beq, Peq] = genGutODE(a, b, c, d, k, r, Ts);

%% Declare initial conditions and target conditions
BTarget = Beq;
PTarget = Peq;

Binit = 100;
Pinit = 5;
%% Declare penalty matrices:
Q=0.5*eye(2);
P=10*eye(2);
R=0.1*eye(1);

%% Declare contraints
% The concentrations must be non-negative
cl=[-Beq; -Peq];
% Select the maximum concentratation to be 
ch=[3*Beq; 2*Peq];

% We cannot feed bacteria at "negative" rate
ul=-r;

% Input cannot be more than triple the rate
uh=r;

% Assume that the variables are independent
D=eye(2);

%% Compute stage constraint matrices and vector
[Dt,Et,bt]=genStageConstraints(A,B,D,cl,ch,ul,uh);

%% Compute trajectory constraints matrices and vector
[DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);

%% Compute QP constraint matrices
[Gamma,Phi] = genPrediction(A,B,N);
[F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);

%% Compute QP cost matrices
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
% Prepare cost and constraint matrices for mpcqpsolver
% Calculating the inverse of the lower triangular H. see doc mpcqpsolver.
H = chol(H,'lower');
H=(H'\eye(size(H)))';

%% Run a linear simulation
xTarget=[BTarget-Beq PTarget-Peq]';% target equilibrium state
x0=[Binit-Beq Pinit-Peq]'; % starting weights of bacteria
iA = false(size(bb));
t=0:Ts:T;
x=[x0, zeros(2,length(t)-1)];
for t_step=1:length(t)-1
    [u,~,iA] = genMPController(H,G,F,bb,J,L,x(:,t_step),xTarget,1,iA);
    x(:,t_step+1)=A*x(:,t_step)+B*u;
    u_vec(:,t_step) = u;
end

%% Plot results
% Define some plotting variables
B = Beq + x(1,:);
P = Peq + x(2,:);
U = [zeros(1), u_vec];

Umin = ul*ones(1,length(t));
Umax = uh*ones(1,length(t));

BTarget_plot = BTarget*ones(1,length(t));
PTarget_plot = PTarget*ones(1,length(t));

% Make the actual plot
figure(1);
subplot(3,1,1)
hold on
plot(t,BTarget_plot,'k--','Linewidth',2)
plot(t, B,'Linewidth',2)
title('Prey Bacteria')
% xlabel('Time [Days]')
ylabel('Bacteria Mass [g]')
ax1 = gca;
ax1.YLim  = [0 110];

subplot(3,1,2)
hold on
plot(t,PTarget_plot,'k--','Linewidth',2)
plot(t, P,'Linewidth',2)
title('Predator Bacteria')
% xlabel('Time [Days]')
ylabel('Bacteria Mass [g]')
ax2 = gca;
ax2.YLim  = [0 110];

subplot(3,1,3)
hold on
pUmin = plot(t,Umin,'--','Linewidth',2,'Color','Red');
pUmax = plot(t,Umax,'--','Linewidth',2,'Color','Red');
pU = plot(t,U,'b.','markersize',10,'Linewidth',2);
title('Feeding rate')
xlabel('Time [Days]')
ylabel('Rate')
ax3 = gca;
ax3.YLim  = [-2 5];
ax3.Children(1).Color = [0 0.4470 0.7410];