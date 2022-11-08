% Tanner Huck    
% Math 301 B
% Homework 5 part 2

clear; clc; close all;

%% Question 2
% Part a
load("salmon_data.mat");

% Part b
% Finding the line of best fit for degree 1
lin_fit_degree1 = polyfit(year, salmon, 1);
A7 = lin_fit_degree1;

% Part c
% Finding the line of best fit for degree 3
lin_fit_degree3 = polyfit(year, salmon, 3);
A8 = lin_fit_degree3;

% Part d
% Finding the line of best fit for degree 5
lin_fit_degree5 = polyfit(year, salmon, 5);
A9 = lin_fit_degree5;

% Part e
% Finding the error for each of the different lines
err1 = abs(A7(1).*2021+A7(2)-489523)/489523; 
err2 = abs(A8(1).*2021.^3+A8(2).*2021.^2+A8(3).*2021+A8(4)-489523)/489523; 
err3 = abs(A9(1).*2021.^5+A9(2).*2021.^4+A9(3).* ...
    2021.^3+A9(4).*2021.^2+A9(5).*2021+A9(6)-489523)/489523; 
A10 = [err1,err2,err3];

%% WriteUp Question 2
% Part a
% Ploting the salmon data
plot(year, salmon, '-k.') 
hold on
% Plotting degree 1 best fit line
plot(year, A7(1).*year+A7(2), 'b', 'linewidth', 2)
hold on
% Plotting degree 3 best fit line
plot(year, A8(1).*year.^3+A8(2).*year.^2+A8(3).*year+A8(4), ...
    'r', 'linewidth', 2)
hold on
% Plotting degree 5 best fit line
plot(year, A9(1).*year.^5+A9(2).*year.^4+A9(3).*year.^3+A9(4).*year.^2 ...
    +A9(5).*year+A9(6), 'm', 'linewidth', 2)
% Setting axis bounds
xlim([1930,2020])
ylim([100000,1500000])
% Adding labels, title, and legend to the plot
title('Salmon Population and best fit lines')
xlabel('Year');
ylabel('Number of Salmon');
legend('Number of Salmon','Best-fit polynomial of degree 1', ...
    'Best-fit polynomial of degree 3','Best-fit polynomial of degree 5',...
    'NW');

% Part e
% Finding the predicted populatijn in 2050
degree1_prediction2050 = A7(1).*2050+A7(2); 
degree3_prediction2050 = A8(1).*2050.^3+A8(2).*2050.^2+A8(3)...
    .*2050+A8(4);
degree5_prediction2050 = A9(1).*2050.^5+A9(2).*2050.^4+A9(3).* ...
    2050.^3+A9(4).*2050.^2+A9(5).*2050+A9(6); 










