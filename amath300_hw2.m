% Tanner Huck    
% Math 301 A
% Homework 3

clear; clc;

%% Question 1
% Using the bisectin method to find the root of a function 
% 0 = f/(1-f) * sqrt(5/(2+f)) - 0.06from -0.15 to 0.7 with a 
% tolerance of 10^-10

% By using the bijection function I created later at the end,
% A1 is the root of the function and A2 is the number of steps it took to 
% get to that root
[A1, A2] = bisection(@(f) f/(1-f) * sqrt(5/(2+f)) - 0.06, -1.5, 0.7, 10.^-10);

%% Question 2
% part a
% Calculating the eact root of the function f(x)=2x^3-3x^2-9 using fzero
% with initial quess of 2.5
f = @(x) 2*x.^3 - 3*x.^2 - 9;
A3 = fzero(f, 2.5);

%% part b
% Using newtons method to find the root with 100 steps
finalVal = 2.5;
df = @(x) 6*x.^2 - 6*x;
for i = 1:100
    finalVal = finalVal - (f(finalVal)/(df(finalVal)));
end       
% Calulating the error of newtons method compared to the exact value
A4 = abs(A3 -finalVal);

%% part c
% Calculating how many steps of newtons method it takes to find the root
% within a certain tolerance
A5 = 0;
tol = 10.^-12;
newFinalVal = 2.5;
for i = 1:100
    if (abs(A3 - newFinalVal) < tol) 
        break
    end
    newFinalVal = newFinalVal - (f(newFinalVal)/(df(newFinalVal)));
    A5 = A5 + 1; 
end


%% Functions
% Function for bisections
function [root_value, iterations] = bisection(fun, left, right, tol)
    for k = 1:1000
        mid = (left + right)/2; % Finds the midpoint between two points        
        if abs(fun(mid)) < tol % If the midpoint is whithin our tol, we 
            % say this is our root value
            root_value = mid; 
            iterations = k; % Setting k to how many steps it took to find 
            % the root
            break % Stop, because we found our root
        elseif fun(left)*fun(mid)<0 
            right = mid;
            % Otherwise, if the point to the left of the mid point we set 
            % the mid to the new right point
        elseif fun(mid)*fun(right)<0 
            left = mid;
            % Otherwise, if the point to the right of the mid point we set 
            % the mid to the new left point
        else
            break
            % If it does none of these things it will stop and display no
            % root found
        end
    end
end









