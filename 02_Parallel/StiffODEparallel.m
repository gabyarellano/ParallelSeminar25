function [X,Y,Z,CompTime] = StiffODEparallel(gSize)
% Function to compute parameter sweep study of the van der Pol equation.
% Input: 
%  gSize: Number of parameters nu and mu 
% Outpus:
   %CompTime: parfor-loop computation time
   %[X,Y,Z]: Grid of values from ODE solution
   %         Z is the mean period of the oscillator
% Copyright 2023 The MathWorks, Inc.

% Problem Setup
x = linspace(100, 150, gSize); % Parameter mu values
y = linspace(0.5, 2, gSize); % Parameter nu values
[X,Y] = meshgrid(x,y);
Z = nan(size(X));
% Solve Equation
t0 = tic;
parfor ii = 1:numel(X)
    [t, yy] = solveVdp(X(ii), Y(ii));
    l = islocalmax(yy(:, 2));
    Z(ii) = mean(diff(t(l))); 
end

CompTime = toc(t0);
end