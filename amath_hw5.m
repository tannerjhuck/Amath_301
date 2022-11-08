% Tanner Huck    
% Math 301 B
% Homework 5 part 1

clear; clc; close all;

%% Coding Question 1
% Part a
% loading in the data
load("CO2_data.mat");

% Part b
% Finding the optimal parameters to get the least amount of error usuing
% the sum of squared error method
adapter = @(vars) sumSquaredError(vars(1), vars(2), vars(3), year);
min_vals_squared = fminsearch(adapter, [30,0.03,300]);
A1 = min_vals_squared;

% Part c
% Saving the total Sum of Squared Errors
A2 = sumSquaredError(A1(1),A1(2),A1(3),year);

% Part d
% Finding optimal parameters to get least amount of error, but using the
% average error method
adapter_avg = @(vars) avgError(vars(1), vars(2), vars(3), year);
min_vals_avg = fminsearch(adapter_avg, [30,0.03,300]);
A3 = min_vals_avg;

% Part e
% Finding optimal parameters to get least amount of error, but using the
% max error method
adapter_max = @(vars) maxError(vars(1), vars(2), vars(3), year);
min_vals_max = fminsearch(adapter_max, [30,0.03,300]);
A4 = min_vals_max;

% Part f
% Finding the optimal parameters to get the least amount of error usuing
% the sum of squared error method with a new equation
adapter_newSquared = @(vars) sumSquaredErrorNew(vars(1), vars(2), ...
    vars(3), vars(4), vars(5), vars(6), year);
options = optimset('MaxFunEvals', 2000);
min_vals_squared_new = fminsearch(adapter_newSquared, ...
    [A1(1),A1(2),A1(3),-5,4,0], options);
A5 = min_vals_squared_new;

% Part g
% Saving the total Sum of Squared Errors
A6 = sumSquaredErrorNew(A5(1),A5(2),A5(3),A5(4),A5(5),A5(6),year);

%% Write Up Question 1
% Part a
% Ploting the exponential + sinusoidal fit
plot(year, (A5(1).*exp(A5(2).*year))+A5(3)+(A5(4).*sin(A5(5).* ...
    (year-A5(6)))), 'b', 'LineWidth', 2);
hold on
% Plotting the Atmospheric CO_2
plot(year, CO2, '-k.')
hold on
% Plotting the exponential fit
plot(year, (A1(1).*exp(A1(2).*year))+A1(3), 'r', 'LineWidth', 2);
xlim([0, 63])
% Adding labels, title, and legend to the plot
title('Atmospheric CO_2 and two lines of best fit')
xlabel('Years since January 1958');
ylabel('Atmospheric CO_2');
legend('exponential + sinusoidal fit','Atmospheric CO_2', ...
    'exponential fit','location', 'NW');

% Part d
% Finding with fit is better at soloving atmospheric Co2 in Jan, 2022
ex_sin = (A5(1).*exp(A5(2).*(64 +(1/12))))+A5(3)+(A5(4).* ...
    sin(A5(5).*(((64 +(1/12)))-A5(6))));
ex = (A1(1).*exp(A1(2).*(64 +(1/12))))+A1(3);

% Part e
% Calculating the average predicted CO2 for all of 2022
avg_ex = 0; 
for i = 1:12
    avg_ex = avg_ex + (A1(1).*exp(A1(2).*(63 +(i/12))))+A1(3)
end
avg_ex = avg_ex/12;

avg_ex_sin = 0; 
for i = 1:12
    avg_ex_sin = avg_ex_sin + (A5(1).*exp(A5(2).*(63 +(i/12))))+A5(3) ...
    +(A5(4).*sin(A5(5).*(((63 +(i/12)))-A5(6))))
end
avg_ex_sin = avg_ex_sin/12;

%% Functions
% Sum of squared error
function error = sumSquaredError(a, r, b, t)
    load('CO2_data.mat');
    y = @(x) (a.*exp(r.*x))+b;
    error = sum( (CO2 - y(t)).^2 );
end

% Average error
function error_avg = avgError(a, r, b, t)
    load('CO2_data.mat');
    y = @(x) (a.*exp(r.*x))+b;
    error_avg = sum( abs(CO2 - y(t)) );
end

% Max error
function error_max = maxError(a, r, b, t)
    load('CO2_data.mat');
    y = @(x) (a.*exp(r.*x))+b;
    error_max = max( abs(CO2 - y(t)) );
end

% Sum of squared error new equation
function error_newSumSquared = sumSquaredErrorNew(a, r, b, c, d, e, t)
    load('CO2_data.mat');
    y = @(x) (a.*exp(r.*x))+b+(c.*sin(d.*(x-e)));
    error_newSumSquared = sum( (CO2 - y(t)).^2 );
end