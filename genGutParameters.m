%% 26 May 2019 Miroslav Gasparek
% Data for the MPC of the Ligand/Receptor interaction as described in
% "An Introduction to Feedback Control in Systems Biology" by 
% Carlo Cosentino and Declan Bates.
%
% Data in the book come from various sources

%% Parameters of the model
a = 3.2;
b = 0.6;
c = 50;
d = 0.56;
k = 125;
r = 1.6;

% Save the parameters 
save('gut_parameters.mat','a','b','c','d','k','r');