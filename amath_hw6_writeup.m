% Tanner Huck    
% Math 301 B
% 5/19/2022;
% Homework 6 Writeup

clear; clc; close all;

%% Coding Question 1
% Part a
% deining Himmelblau's function, same as coding section
Him_fun = @(x,y) (x.^2+y-11).^2+(x+y.^2-7).^2;
Him_fun = @(p) Him_fun(p(1), p(2));
% creating a mesh to be used in the gaph
points = linspace(-5,5,40);
[x,y] = meshgrid(points, points);
z = (x.^2+y-11).^2+(x+y.^2-7).^2;
% graphing the Himmelblau's function
surf(x, y, z)
% adding things to make the grah look better
colormap('hot') % make it hot colors
colorbar % adding colorbar
view(30,30) % changing the view
shading interp % getting rid of mesh lines
% adding labels
title('Himmelblau’s function')
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');

%% Part b
% new figure
figure
% making a new mesh for a contour plot
points =linspace(-5,5,100);
[x,y] = meshgrid(points, points);
z = (x.^2+y-11).^2+(x+y.^2-7).^2;
% makeing the contour plot
contourf(x, y, z, logspace(-1, 3, 22))
% using a logarithmic spacing of the colormap to focus more on larger changes
set(gca, 'colorscale', 'log')
% changing color scale
colormap('jet')
% adding colorbar
colorbar
% adding labels
title('Himmelblau’s function and Minimum Points')
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');

% Part c
% using fminsearch to find the 4 minimun points and adding them to the
% contour plot as yellow stars
x01 = [-4,-3,0];
min1 = fminsearch(Him_fun, x01);
hold on
plot(min1(1),min1(2), 'yp', 'linewidth', 5, 'MarkerSize',15)
x02 = [-3,3,0];
min2 = fminsearch(Him_fun, x02);
hold on
plot(min2(1),min2(2), 'yp', 'linewidth', 5, 'MarkerSize',15)
x03 = [3,2,0];
min3 = fminsearch(Him_fun, x03);
hold on
plot(min3(1),min3(2), 'yp', 'linewidth', 5, 'MarkerSize',15)
x04 = [4,-2,0];
min4 = fminsearch(Him_fun, x04);
hold on
plot(min4(1),min4(2), 'yp', 'linewidth', 5, 'MarkerSize',15)

% Part d
% code from problem 2a of coding section
initial_guess = [-3;-2];
A3 = fminsearch(Him_fun, initial_guess);
% adding this to the contour plot as a red circle
hold on
plot(A3(1),A3(2), 'ro', 'linewidth', 2, 'MarkerSize', 7)

% code from 2 of coding section
fgrad = @(x,y) [4.*x.^3-42.*x+4.*x.*y+2.*y.^2-14; 
               4.*y.^3-26.*y+4.*x.*y+2.*x.^2-22]; 
fgrad = @(p) fgrad(p(1), p(2));
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
A5 = p;

% adding this to the contour plot as a green square
hold on
plot(A5(1),A5(2), 'gs', 'linewidth', 2, 'MarkerSize', 15)

% adding a legend for these last two points
legend('','','','','','Coding 2a minimum', 'Grad descent')

%% Problem 2
% Part a
tic
% Timing to see how long it takes to do gradient descent
% crating the tolerance
tol = 10e-9;
% for loop that will run 12000 iterations of gradient descent
p = [2; 3]; % Choose an initial guess
for i = 1:12000
% one step in grad descent
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
iterations_grad_descent = i;
time_grad_descent = toc;

% Part b
tic
% Timing to see how long it takes to do gradient descent with fixed t
tol = 10e-9; % crating the tolerance
% for loop that will run 12000 iterations of gradient descent
p = [2; 3]; % Choose an initial guess
tstep = 0.01; % deining tstep, the fixed value of t
for i = 1:12000
% one step in grad decent
grad = fgrad(p); % Find which direction to go
% replace old code with this line instead, uses a fixed value of t instead
% of calculating tmin with fminbnd for every step 
p = p - tstep*grad;
% calculate the infinity norm of vector grad
infinity_norm = norm(grad, inf);
% stop if the infinity vector is less than tolerance
if (infinity_norm < tol) 
    break
end
end
iterations_tstep_01 = i;
time_tstep_01 = toc;

% Part c
% same as part b, but with tstep 0.02
tic
% Timing to see how long it takes to do gradient descent with fixed t
tol = 10e-9; % crating the tolerance
% for loop that will run 12000 iterations of gradient descent
p = [2; 3]; % Choose an initial guess
tstep = 0.02; % deining tstep, the fixed value of t
for i = 1:12000
% one step in grad decent
grad = fgrad(p); % Find which direction to go
% replace old code with this line instead, uses a fixed value of t instead
% of calculating tmin with fminbnd for every step 
p = p - tstep*grad;
% calculate the infinity norm of vector grad
infinity_norm = norm(grad, inf);
% stop if the infinity vector is less than tolerance
if (infinity_norm < tol) 
    break
end
end
iterations_tstep_02 = i;
time_tstep_02 = toc;

% Part d
% repeating part b, but with t step 0.025
tic
% Timing to see how long it takes to do gradient descent with fixed t
tol = 10e-9; % crating the tolerance
% for loop that will run 12000 iterations of gradient descent
p = [2; 3]; % Choose an initial guess
tstep = 0.025; % deining tstep, the fixed value of t
for i = 1:12000
% one step in grad decent
grad = fgrad(p); % Find which direction to go
% replace old code with this line instead, uses a fixed value of t instead
% of calculating tmin with fminbnd for every step 
p = p - tstep*grad;
% calculate the infinity norm of vector grad
infinity_norm = norm(grad, inf);
% stop if the infinity vector is less than tolerance
if (infinity_norm < tol) 
    break
end
end
iterations_tstep_025 = i;
time_tstep_025 = toc;

% Part e
% creating a table to store the results of the times and iterations of
% differenct methods of gradient discent
results_table = table(["tstep=0.01"; "tstep=0.02"; "tstep=0.025"; ...
    "fminbound"], ...
    [iterations_tstep_01; iterations_tstep_02; iterations_tstep_025; ...
    iterations_grad_descent], ...
    [time_tstep_01; time_tstep_02; time_tstep_025; time_grad_descent], ...
    ["Yes"; "Yes"; "No"; "Yes"])
results_table.Properties.VariableNames = {'Method', ...
    'Number Iterations','Time', 'Converged'}



