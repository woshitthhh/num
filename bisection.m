function [root, iter] = bisection(f, a, b, tol, maxiter)
% Bisection method for finding a root of f(x) = 0
% Inputs:
%   f: function handle for f(x)
%   a, b: interval [a, b] containing the root
%   tol: tolerance for stopping criterion |f(x)| < tol
%   maxiter: maximum number of iterations
% Outputs:
%   root: approximate root of f(x) = 0
%   iter: number of iterations used

% Check if f(a) and f(b) have opposite signs
if sign(f(a)) == sign(f(b))
    error('f(a) and f(b) must have opposite signs')
end

% Initialize variables
iter = 0;
fa = f(a);
fb = f(b);

% Main loop
while iter < maxiter
    % Compute midpoint and function value at midpoint
    c = (a + b) / 2;
    fc = f(c);
    
    % Check if |f(c)| < tol
    if abs(fc) < tol
        root = c;
        return
    end
    
    % Update interval
    if sign(fc) == sign(fa)
        a = c;
        fa = fc;
    else
        b = c;
        fb = fc;
    end
    a 
    fa
    b
    fb
    a-b
    % Increment iteration counter
    iter = iter + 1;
end

% Maximum number of iterations reached without convergence
warning('Maximum number of iterations reached without convergence')
root = c;
end
