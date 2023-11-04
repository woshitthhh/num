function [root, iter] = newton(f, df, x0, tol, maxiter)
% Newton's method for finding a root of f(x) = 0
% Inputs:
%   f: function handle for f(x)
%   df: function handle for f'(x)
%   x0: initial guess for root
%   tol: tolerance for stopping criterion |f(x)| < tol
%   maxiter: maximum number of iterations
% Outputs:
%   root: approximate root of f(x) = 0
%   iter: number of iterations used

% Initialize variables
iter = 0;
x = x0;

% Main loop
while iter < maxiter
    % Compute function value and derivative at x
    fx = f(x)
    dfx = df(x);
    
    % Check if |f(x)| < tol
    if abs(fx) < tol
        root = x;
        return
    end
    
    % Update x using Newton's method
    x = x - fx / dfx
    
    % Increment iteration counter
    iter = iter + 1;
end

% Maximum number of iterations reached without convergence
warning('Maximum number of iterations reached without convergence')
root = x;
end
