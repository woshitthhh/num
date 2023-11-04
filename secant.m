function [root, iter] = secant(f, x0, x1, tol, maxiter)
% Secant method for finding a root of f(x) = 0
% Inputs:
%   f: function handle for f(x)
%   x0, x1: initial guesses for root
%   tol: tolerance for stopping criterion |f(x)| < tol
%   maxiter: maximum number of iterations
% Outputs:
%   root: approximate root of f(x) = 0
%   iter: number of iterations used

% Initialize variables
iter = 0;
fx0 = f(x0);
fx1 = f(x1);

% Main loop
while iter < maxiter
    % Compute next approximation using secant method
    x = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    fx = f(x);
    
    % Check if |f(x)| < tol
    if abs(fx) < tol
        root = x;
        return
    end
    
    % Update variables for next iteration
    x0 = x1
    fx0 = fx1
    x1 = x;
    fx1 = fx;
    
    % Increment iteration counter
    iter = iter + 1;
end

% Maximum number of iterations reached without convergence
warning('Maximum number of iterations reached without convergence')
root = x;
end
