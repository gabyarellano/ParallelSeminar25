function updateSurface(s, d)
%
% Based on:
%    https://www.mathworks.com/help/distcomp/examples/...
%    plot-progress-during-parallel-computations-using-parfor-and-dataqueue.html

% Copyright 2021, The MathWorks, Inc.

s.ZData(d(1)) = d(2);
drawnow('limitrate');
end