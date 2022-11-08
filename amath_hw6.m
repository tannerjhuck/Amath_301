% Tanner Huck    
% Math 301 B
% Homework 6

clear; clc; close all;

%% Coding Question 1
% Part a
% Finding the max of the function by taking fzero of the derivative with
% initial guess 2
equation = @(t) (10/3)*(exp(-t/24)-exp(-t/2));
derivative = @(t) (-5/36)*exp((-t)/2)*( exp((-11*t)/24)-12);
root = abs(fzero(derivative, 2));
A1 = [root, equation(root)];

% Part b
% Finding the max of the function again, but using fminbnd
negative_equation= @(t) -(10/3)*(exp(-t/24)-exp(-t/2));
root_bnd= fminbnd(negative_equation, 0, 12);
A2 = [root_bnd; equation(root_bnd)];

%% Question 2
% Part a
% finding the col vector that mimimizes f using fminsearch
Him_fun = @(x,y) (x.^2+y-11).^2+(x+y.^2-7).^2;
Him_fun = @(p) Him_fun(p(1), p(2));
initial_guess = [-3;-2];
A3 = fminsearch(Him_fun, initial_guess);

% Part b
% creating a function that calculates the gradient
fgrad = @(x,y) [4.*x.^3-42.*x+4.*x.*y+2.*y.^2-14; 
               4.*y.^3-26.*y+4.*x.*y+2.*x.^2-22]; 
fgrad = @(p) fgrad(p(1), p(2));
A4 = norm(fgrad(A3));

% Part c
% crating the tolerance
tol = 10e-8;
% for loop that will run 5000 iterations of gradient descent
p = [-3; -2]; % Choose an initial guess
for i = 1:5000
% one step in grad decent
grad = fgrad(p); % Find which direction to go
phi = @ (t) p - t*grad; % Define the path
f_of_phi = @ (t) Him_fun(phi(t)); % Create a function of "heights along path"
tmin = fminbnd(f_of_phi,0,1); % Find time it takes to reach min height
p = phi(tmin); % Find the point on the path and update your guess

% calculate the infinity norm of vector grad
infinity_norm = norm(grad, inf);
% stop if the infinity vector is less than tolerance
if (infinity_norm < tol) 
    break
end
end

% Part d
A5 = p;
A6 = i;

















