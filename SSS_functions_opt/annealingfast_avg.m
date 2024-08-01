function newx = annealingfast_avg(optimValues,problem)
%ANNEALINGFAST Generates a point at distance equal to temperature
%   NEWX = ANNEALINGFAST(optimValues,problem) generates a
%   (multidimensional) point equal to the current point plus a random
%   step of size equal to the current temperature, modified to remain
%   within bounds.
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter
%
%   PROBLEM is a structure containing the following information:
%      objective: function handle to the objective function
%             x0: the start point
%           nvar: number of decision variables
%             lb: lower bound on decision variables
%             ub: upper bound on decision variables
%
%   Example:
%    Create an options structure using ANNEALINGFAST as the annealing
%    function
%    options = optimoptions('simulannealbnd','AnnealingFcn',@annealingfast);


%   Copyright 2006-2015 The MathWorks, Inc.

currentx = optimValues.x;
nvar = numel(currentx)/4;
newx = currentx;
y = randn(nvar,1);
y = y./norm(y);
y_144 = nvar*4;

for j=1:nvar
    if j==1
        y_144(j:4*j,:) = repmat(y(j),4,1);
    else
        y_144(4*j-3:4*j,:) = repmat(y(j),4,1);
    end
end
newx(:) = currentx(:) + optimValues.temperature.*y_144;
% newx(:) = currentx(:) + temperature.*y_144;

newx = bounds(newx,optimValues,problem);
