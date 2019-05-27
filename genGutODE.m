function [A,B,C,D, Beq, Peq] = genGutODE(a, b, c, d, k, r, Ts)
% Inputs:

% Outputs:
% A,B,C,D = State Space matrices of a discrete-time or continuous-time state space model

% replace A,B,C,D with the correct values
A = zeros(2,2);
B = zeros(2,1);
C = zeros(2,2);

%% Compute equilibrium points
Beq = c*d/(a*b-d);
Peq = b*c*r*(a*b*k - c*d - d*k)/(k*(a*b-d)^2);

%% Compute the matrices of the linearized system
% State matrix
A(1,1) = r*(1 - 2*Beq/k) - a*Peq*(1/(c+Beq) - Beq/(c+Beq)^2);
A(1,2) = -a*Beq/(c+Beq);
A(2,1) = b*a*Peq*(1/(c+Beq) - Beq/(c+Beq)^2);
A(2,2) = b*a*Beq/(c+Beq) - d;

% Input matrix
B(1) = Beq*(1-Beq/k);
B(2) = 0;

% Observation matrix, consider that we can only observe the predator
% population
C=[0, 1;
   0, 0];

% Direct matrix
D=zeros(2,1);

% if Ts>0 then sample the model with a zero-order hold (piecewise constant) input, otherwise return a continuous-time model
if Ts>0
    system = ss(A,B,C,D);
    system_d = c2d(system, Ts);
    
    A = system_d.A;
    B = system_d.B;
    C = system_d.C;
    D = system_d.D;
end

end