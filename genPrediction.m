function [Gamma, Phi] = genPrediction(A,B,N)
% GENPREDICTION  [Gamma,Phi] = genPrediction(A,B,N). 
% A and B are discrete-time state space matrices for x[k+1]=Ax[k]+Bu[k]
% N is the horizon length. 
% Your code is suppose to work for any linear system, not just the gantry crane. 

% Write your code here
% Pre-allocate matrices

Phi_mat = [];
Gamma = [];

for i = 1:N
    Phi_mat = [Phi_mat; A^i];
end

Gamma_init = [];

for i = 1:N
    Gamma_init = [Gamma_init; (A^(i-1))*B];
end

Gamma_next = Gamma_init;
Gamma = Gamma_next;
for j = 1:N-1    
    Gamma_next = [zeros(size(B)); Gamma_next];
    Gamma_next = Gamma_next(1:end-(size(B,1)),:);
    Gamma = [Gamma, Gamma_next];
end


Phi = Phi_mat;
end