% Tanner Huck    
% Math 301 B
% Homework 7

clear; clc; close all;

%% Coding Question 1
% Part a
% loading in the data
load("Plutonium.csv");
% making vectors for the plutonium and time from the observations
time = Plutonium(1,:);
plutonium = Plutonium(2,:);

% usuing a second difference formula to approximate the derivative at t=0
A1 = (plutonium(3) - 4*plutonium(2) + 3*plutonium(1)) / -2;

% Part b
% usuing a second difference formula to approximate the derivative at t=40
A2 = (3*plutonium(41) - 4*plutonium(40) + plutonium(39)) / 2;

% Part c
% creating a vector to hold all of the deriatives
d_dp = 1:41;
% assigning the derivatives at the endpoints
d_dp(1) = A1;
d_dp(41) = A2;
% for loop using central differance to solve for the derivative at all the
% other points and putting it in the derivative vector
for i = 2:40
    d_dp(i) = plutonium(i+1)/2 - plutonium(i-1)/2;
end
% creating a vector to hold all the decay rates
decay_rate = 1:41;
% for loop to calculate all the decay rates
for i = 1:41
    decay_rate(i) = d_dp(i) * -1/plutonium(i);
end
% mean decay rate
A3 = mean(decay_rate);

% Part d
% calculating the half life
A4 = log(2)/A3;

% Part e
% calculating the second derivative at t=27
A5 = (plutonium(30) + plutonium(26) - 2*plutonium(28)) / 4;

%% Question 2
% Part a
% creating the integral P with given mu and sigma values
mu = 3.5;
sigma = 0.73;
P = @(x) (exp((-((x-mu).^2))/(2.*(sigma.^2))))./(sqrt(2.*pi.*(sigma.^2)));
% calculating the integral from 4 to 5
A6 = integral(P, 4, 5);

% Part b
% vector to hold left-sided rectangle rule to approximate P
left_rect = (1:16);
% Left-rectangle rule
for i = 1:16
    left_rect(i) = (2.^-i)*sum( P(4:2^-i:5-(2^-i)) ); 
end
% transposing the vecotor into a col vector
left_rect = transpose(left_rect);
A7 = left_rect;

% Part c
% vector to hold right-sided rectangle rule to approximate P
right_rect = (1:16);
% Right-rectangle rule
for i = 1:16
    right_rect(i) = (2.^-i)*sum( P(4+(2^-i):2^-i:5) ); 
end
% transposing the vecotor into a col vector
right_rect = transpose(right_rect);
A8 = right_rect;

% Part d
% vector to hold midpoint rule to approximate P
midpoint = (1:16);
% midpoint rule
for i = 1:16
    midpoint(i) = (2.^-i)*sum( P(4+((2^-i)/2):2^-i:5-((2^-i)/2)) ); 
end
% transposing the vecotor into a col vector
midpoint  = transpose(midpoint );
A9 = midpoint;

% Part e
% vector to hold rapezoidal rule to approximate P
trapezoidal = (1:16);
% trapezoidal rule
for i = 1:16
    trapezoidal(i) = (2.^-i)/2*(P(4)+P(5)+2*...
        sum( P(4+(2^-i):2^-i:5-(2^-i)) )); 
end
% transposing the vecotor into a col vector
trapezoidal = transpose(trapezoidal);
A10 = trapezoidal;

% Part f
% vector to hold simpson rule to approximate P
simpson = (1:16);
% simpson rule
for i = 1:16
    simpson(i) = (2.^-i)/3*(P(4)+P(5)+ 2*...
        sum( P(4+2*(2^-i):2*(2^-i):5-2*(2^-i)) ) + 4*...
        sum( P(4+(2^-i):2*(2^-i):5-(2^-i)) )); 
end
% transposing the vecotor into a col vector
simpson = transpose(simpson);
A11 = simpson;

%% Writeup Question 1
% Part a
% calculating the error of each method for each step side
error_left = abs(left_rect - A6);
error_right = abs(right_rect - A6);
error_mid = abs(midpoint - A6);
error_trap = abs(trapezoidal - A6);
error_simp = abs(simpson - A6);

% Part b
% creating a vector with all the step sizes
step_size = (1:16);
for i = 1:16
    step_size(i) = (2^-i);
end
% creating a log log plot for the error of each method versus the step size 
loglog(step_size, error_left, 'ko','MarkerSize',12) 
hold on
loglog(step_size, error_right, 'm+','MarkerSize',12,'linewidth', 0.75)
hold on
loglog(step_size, error_mid, 'g*','MarkerSize',12)
hold on
loglog(step_size, error_trap, 'bx','MarkerSize',12)
hold on
loglog(step_size, error_simp, 'r+','MarkerSize',12)
% adding a trend line that represents O(h)
hold on 
loglog(step_size, 0.2*step_size, 'k', 'linewidth', 2, 'Linestyle', '-.')
% adding a trend line that represents O(h^2)
hold on 
plot(step_size, 0.01*(step_size.^2), 'k', 'linewidth', 2, 'Linestyle', ':')
% adding a trend line that represents O(h^3)
hold on 
plot(step_size, 0.008*(step_size.^4), 'k', 'linewidth', 2, 'Linestyle', '--')
% adding a line where machine precision is
hold on
yline(10^-16, 'k', 'linewidth', 1)
% adding a legend
legend('Left Rectangle', 'Right Rectangle', 'Midpoint', 'Trapazoidal', ...
    'Simpson', 'O(h)', 'O(h^2)', 'O(h^4)', 'Machine precision', ...
    'location', 'best')
% adding a title and axis labels
title('Convergence of Numerical Integration Schemes')
xlabel('Step Size');
ylabel('Error');

%% Question 2. Extra credit
% Part a
% creating the function for the circumference with given a and b values
% (semi-major and semi-minor axis')
a = 1.5;
b = 0.3;
Cir = @(x) 4.*sqrt(((a.^2).*((cos(x)).^2))+((b.^2).*((sin(x)).^2)));
% calculating the integral of the function from 0 to pi/2
true = integral(Cir, 0, pi/2);

% Part b
% vector to hold left-sided rectangle rule to approximate Cir
left_rect_cir = (1:14);
% Left-rectangle rule
for i = 1:14
    step = 2^i;
    step_size = (pi/2)/step;
    left_rect_cir(i) = step_size*sum( Cir(linspace(0,pi/2,step)) ); 
end
% transposing the vecotor into a col vector
left_rect_cir = transpose(left_rect_cir);

% Part c
% vector to hold trapezoidal rule to approximate Cir
trapezoidal_cir = (1:14);
% trapezoidal rule
for i = 1:14
    step = 2^i;
    step_size = (pi/2)/step;
    points = linspace(0,pi/2,step);
    trapezoidal_cir(i) = step_size/2*(Cir(0)+Cir(pi/2)+2* ...
        sum( Cir(points(2:end-1)) )); 
end
% transposing the vecotor into a col vector
trapezoidal_cir = transpose(trapezoidal_cir);

% Part d
% calculating the error of each method for each step side
error_left_cir = abs(left_rect_cir - true);
error_trap_cir = abs(trapezoidal_cir - true);
% creating a vector with all the step sizes
step_size = (1:14);
for i = 1:14
    step_size(i) = (pi/2)/(2^i);
end
figure
% creating a log log plot for the error of each method versus the step size
loglog(step_size, error_left_cir, 'b*','MarkerSize',12) 
hold on
loglog(step_size, error_trap_cir, 'g+','MarkerSize',12,'linewidth', 0.75)
% adding a trend line that represents O(h)
hold on 
loglog(step_size, 0.4*step_size, 'k', 'linewidth', 2, 'Linestyle', '-.')
% adding a trend line that represents O(h^2)
hold on 
plot(step_size, 4*(step_size), 'k', 'linewidth', 2, 'Linestyle', ':')
% adding a legend
legend('Left Rectangle', 'Trapazoidal', 'O(h) for left-rectangle', ...
    'O(h) for trapazoid', 'location', 'best')
% adding a title and axis labels
title('Convergence of Numerical Integration Schemes')
xlabel('Step Size');
ylabel('Error');



