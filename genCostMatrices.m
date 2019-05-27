function [H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N)
%% cost function matrices
% Gamma and Phi are the prediction matrices
% Q is the stage cost weight on the states, i.e. x'Qx
% R is the stage cost weight on the inputs, i.e. u'Ru
% P is the terminal weight on the final state

% Your code goes here
I_N_Rmat = eye(N);
I_N_Qmat = eye(N-1);
% Get R_mat
R_mat = kron(I_N_Rmat,R);


% Get Q_mat
Q_submat = kron(I_N_Qmat,Q);
Q_mat = [Q_submat, zeros(size(Q_submat,1), size(P,2));
         zeros(size(P,1), size(Q_submat,2) ), P];

% Get the H-matrix, weight on the input quadratic form     
H = (Gamma')*Q_mat*Gamma + R_mat;

G = (((Phi')*Q_mat*Gamma)');

end