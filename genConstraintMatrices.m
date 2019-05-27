function [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N)

% Get the length of the state vector
len_xvec = size(Phi,1)/N;
% Get the length of the input vector
len_uvec = size(Gamma,2)/N;

% Get the "width" of the Gamma matrix
width_uvec = size(Gamma,2);

% Make the identity matrix that will be placed 
% at the top of the Phi_tilde vector
I_Phi = eye(len_xvec);

% Make the zero matrix that will be placed
% at the top of the Gamma_tilde vector
% Z_Gamma = zeros(len_uvec, width_uvec);
Z_Gamma = zeros(len_xvec, width_uvec);

% Compute Phi_tilde
Phi_cut = Phi(1:end-len_xvec,:); % Takes out the last matrix
Phi_tilde = [I_Phi; Phi_cut];

% Compute Gamma_tilde
% Gamma_cut = Gamma(1:end-len_uvec,:);
Gamma_cut = Gamma(1:end-len_xvec,:);
Gamma_tilde = [Z_Gamma; Gamma_cut];

% Compute the constraint matrices
F = DD*Gamma_tilde + EE;
J = -(DD*Phi_tilde);
L = DD*(Phi_tilde - (kron(ones(N,1), eye(len_xvec))));

end