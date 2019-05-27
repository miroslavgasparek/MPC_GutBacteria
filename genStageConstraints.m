function [Dt,Et,bt]=genStageConstraints(A,B,D,cl,ch,ul,uh)

% Do not forget to augment the matrices by input constraints!
Dt = [D*A; -D*A; zeros(size(ul,1), size(A,1)); zeros(size(uh,1), size(A,1))];
Et = [D*B; -D*B; -eye(size(ul,1)); eye(size(uh,1))];
bt = [ch; -cl; -ul; uh];

end