%% 28 May 2019 Miroslav Gasparek
% The ODEs of the predator-prey system 
function dydt = gut_bacteria_ode(t, y, a, b, c, d, k, r)
    dydt = zeros(2,1);
    dydt(1) = r*y(1)*(1-y(1)/k) - a*y(1)*y(2)/(c+y(1));
    dydt(2) = b*a*y(1)*y(2)/(c+y(1)) - d*y(2);
end