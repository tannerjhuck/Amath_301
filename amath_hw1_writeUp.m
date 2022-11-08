% Tanner Huck    
% Math 301 A
% Homework 1 Writeup

clear; clc;

%% Question 1
% A Calculations on the computer are exact if typed in correctly.
%       This is flase, this is becasue the computer will only save to 16
%       decimal places. SO for example if you type in pi, you are not
%       getting the exact value of pi, you are getting pi to 16 decimal
%       places.
% B Calculations done on the computer have an error proportional to the
% machine precision, which is the precision of stored digits on 
% the computer
%       This is true. The computer only saves to a certain amount of
%       decimal places. When you do calculations, the error is due to not
%       containing the other decimal places of the exact value. Meaning
%       that the error will be because of the precision of stored digits.
% C In most cases, when we see a number that is approximately 1e−16, 
% we shouldassume that it is actually 0.
%       This is true. Because when you get a number that is very close to
%       zero. It is most likly that it got somthing other than zero due to
%       computer error.
% D It is impossible to do calculations on any computer with numbers 
% that aresmaller than 1e-16.
%       This is false. For example you can do 0.1e-18 + 0.1e-20 and matlab
%       will give you an answer.

%% Question 2
% a. listed the values of x1-x4 that we found in problem 6 in smallest to
% largest order
x3 = 0;
x4 = 0;
x1 = 1.884836819954217e-08;
x2 = 0.018870549276471;

 % b. Discussing the differences of x values
 % We know that all four of these x vlaues should be 0, but due to computer
 % error, some of them are not exactly zero. I believe that some of
 % these are smaller/larger than others becasue of this error. x3 and x4
 % are the same and did calculate correctly as 0, whereas x2 is a little
 % less at about -0.019 and x1 being the smallest at -1.885e-8. ALthough x1
 % and x2 were not exactly zero, they are still fairly close to zero, so
 % they are not outliers. They may have had some computer error involved
 % when they were calculated. 


% c. why are some of the x values 0 and some are not
% I believe that some of these x values were not correct calculated to 0
% because of computer error. x3 and x4 were summations of 0.25 and 0.5,
% which are stored in matlab as 4-1 and 2-1 respectivly. x1 and x2 were
% summations of 0.1, which are stored as 10-1. My guess is that x3 and x4
% calculated correctly to zero because they have larger values, they have 
% a smaller chance of computer error. Whereas x1 and x2 did not calculate 
% correctly because they have smaller values, which may have a higher
% chance of error.


%% Question 3
% Ploting the n-th order Taylor-series approximation to cosine
xplot = -pi:0.01:pi; % Create an array from -pi to pi in step sizes of 0.01
yplot = cos(xplot);

% Creating three variables which will later be our vector of taylor series
% expansion
n_values_1 = 1; 
n_values_2 = 1; 
n_values_3 = 1;

% For loop to caclulate each of the taylor series expansions
% Don't need the first one because it is only has one calculation
n_values_1 = n_values_1 + (((-1).^1)*((xplot).^(2*1)))/factorial(2*1);
for i=1:3
    n_values_2 = n_values_2 + (((-1).^i)*((xplot).^(2*i)))/factorial(2*i);
end
for i=1:14
    n_values_3 = n_values_3 + (((-1).^i)*((xplot).^(2*i)))/factorial(2*i);
end

% Creating a plot for cos(x) and its Taylor approximations
plot(xplot,yplot,'LineWidth',2,'Color','k');
hold on 
plot(xplot, n_values_1,'LineWidth',2,'Color','b','LineStyle','--')
hold on
plot(xplot, n_values_2,'LineWidth',2,'Color','r','LineStyle','-.')
hold on 
plot(xplot, n_values_3,'LineWidth',2,'Color','m','LineStyle',':')

% Giving the plot axis labels
xlabel('x-values') 
ylabel('Approximations')

% Giving the plot a title
title('cos(x) and its Taylor approximations')

% Giving the plot a legend for all four of the lines
% Also moving the legend to the south location
legend('cos(x)','n = 1 Taylor approximation', ...
    'n = 3 Taylor approximation', ...
    'n = 14 Taylor approximation','Location', 'south')

% Setting the font of all labels to 10
set(gca, 'FontSize', 10)












