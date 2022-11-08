% Tanner Huck    
% Math 301 B
% 4/21/2022;
% Homework 3 Writeup

clear; clc;

%% Question 1
% part a
% Creating the special 1500x1500 matrix
N = 1500;
A = -1*diag(ones(N,1)) + 4*diag(ones(N-1,1),1) + 4*diag(ones(N-1,1),-1);

%% part b 
% Using a for loop to solve Ax=b 120 times using the backslash method
% and timing it
tic 
residual_backslash = 0;
for i = 1:120
    b = rand(1500,1);
    x = A\b;
    residual_temp1 = A * x - b;
    residual_backslash = residual_backslash + residual_temp1;
end
norm_backslash = norm(residual_backslash);
total_backslash = toc;

%% part c 
% Using a loop to solve Ax=b 120 times using the LU decomposition method
% and timing it
tic 
residual_lu = 0;
[L, U, P] = lu(A);
for i = 1:120
    b = rand(1500,1);
    y_LU = L\(P*b);
    x_LU = U\y_LU;
    residual_temp2 = A * x_LU - b;
    residual_lu = residual_lu + residual_temp2;
end
norm_lu = norm(residual_lu);
total_lu = toc;


%% part d
% Using a loop to solve Ax=b 120 times using the matrix inverse method
% and timing it. 
tic 
residual_inverse = 0;
inverse = inv(A);
for i = 1:120
    b = rand(1500,1);
    y = inverse*b;
    residual_temp3 = A * y - b;
    residual_inverse = residual_inverse + residual_temp3;
end
norm_inverse = norm(residual_inverse);
total_inverse = toc;

%% part e
% Edited prvious code to calculate error of each method using the norm
% function. The total error for each method is assigned to the variables
% error_gaussian, error_lu, and error_inverse.

% part f
% Making a table of how long each method took to calculate Ax=b and how
% much error it had.
T = table(["Backslash time"; "LU time"; "inv time"], ...
    [total_backslash; total_lu; total_inverse], ...
    [norm_backslash; norm_lu; norm_inverse])
T.Properties.VariableNames = {'Method','Time for 120 vectors','||r||'}
