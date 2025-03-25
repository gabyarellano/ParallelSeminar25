function [t, y] = solveVdp(mu, nu)
%Solution to Van der Pol equation using ode23s
f = @(~,y) [nu*y(2); mu*(1-y(1)^2)*y(2)-y(1)];
% Modify the max bound, to run the algorithm longer.  Higher value takes
% longer to run. e.g. max = 1000;
max = mu; 
[t,y] = ode23s(f,[0 20*max],[2; 0]);
end
% Copyright 2023, The MathWorks, Inc.

