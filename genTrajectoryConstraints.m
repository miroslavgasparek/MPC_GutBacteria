function [DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N)

% Generate the identity matrix of the appropriate size
I_N = eye(N);
% Generate the vector of ones of the appropriate size
one_vec = ones(N,1);


% Generate the matrix DD
DD = kron(I_N, Dt);

% Generate the matrix EE
EE = kron(I_N, Et);

% Generate the vector bb
bb = kron(one_vec,bt);

end