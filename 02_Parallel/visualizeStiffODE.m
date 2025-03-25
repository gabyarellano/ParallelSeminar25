function visualizeStiffODE(varargin)
% Visualize parameter sweep output
%Inputs:
%  [X,Y,Z]: Grid of values from ODE solution
%           Z is the mean period of the oscillator
% Copyright 2020 The MathWorks, Inc.

if nargin == 3
    % Three inputs
    X = varargin{1};
    Y= varargin{2};
    Z = varargin{3};
elseif nargin == 1    
    % Job output as a cell array
    jobOutput = varargin{1}; % extract job
    assert(iscell(jobOutput),'jobOutput is expected to be a cell');
    % Extract pieces from job
    X = jobOutput{1};
    Y = jobOutput{2};
    Z = jobOutput{3};
else
    error('One or Three inputs expected');
end
figure
surf(X, Y, Z);
xlabel('\mu Values','Interpreter','Tex')
ylabel('\nu Values','Interpreter','Tex')
zlabel('Mean Period of  y_2')
view(137, 30)
axis([100 150 0.5 2 0 500]);
drawnow
end

