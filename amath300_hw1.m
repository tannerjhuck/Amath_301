% Tanner Huck    
% Math 301 A
% Homework 1

clear; clc;

%% Question 1
% Creating a variable A1 with givin elements:
%   cos(0),cos(π/4),cos(π/2),cos(3π/4),cos(π)
A1 = [cos(0), cos(pi/4), cos(pi/2), cos(3*pi/4), cos(pi)];

%% Question 2
% creatig more vectors
% u, 6 evenly spaced elements from 3 to 4
u = linspace(3, 4, 6);

% v, 5 given elements: 0 0.75 1.5 2.25 3.0 3.75 
v = 0:0.75:4;

%% 2a
% vector w with elements of u^3
w = u.^3;
A2= w;

%% 2b
% vector x with elements of tangent(u) plus e^(v)
x = tan(u) + exp(v);
A3 = x;

%% Question 3
% Creating matrix A4 with 2 rows and 4 columns
A4 = [ 1 2 -1 1
       3 -1 6 9];

%% Question 4
% vector z stars with -6, ends at 3, spacing of 1/100
z = (-6:1/100:3);

% vector z1 that consists of every third element in z
z1 = z(3:3:end);
A5 = z1;

%% Question 5
% Creating matrix A, a 21x21 matrix with rows and columns defined by
%       1/(i+j-1)

i_size = 21;
j_size = 21;
A = zeros(i_size, j_size); %creates an empty 21x21 matrix
% nested for loop that will go through every column and then every row
% and update each element to 1/(i+j-1) where i and j determines the element
for i_element=1:i_size 
    for j_element=1:j_size
        A(i_element, j_element) = 1/(i_element+j_element-1);
    end
end

%% 5a
% setting matrix A=A6
A6 = A;

%% 5b
% Creating matrix B which is the same as A, but changing the 9th row to
% the 8th row * 4 and saving it to matrix A7

A(9,:) = A(8,:) * 4;
A7 = A;

%% 5c
% Matrix A8 which is the same as A, but only the last 9 rows and first 7
% columns

A8 = A(end-8:end,1:7);

%% Question 6
% Using for loops to calculate y1 through y4
y_values = zeros(1,4); %creating an empty 1 by 4 vector
% changing element 1 of the vector to the sumation of y1
for i=1:100000
    y_values(:,1) = y_values(:,1) + 0.1;
end
% changing element 2,3,4 of the vector to the sumation of y2,y3,y4
for i=1:100000000
     y_values(:,2) = y_values(:,2) + 0.1;
     y_values(:,3) = y_values(:,3) + 0.25;
     y_values(:,4) = y_values(:,4) + 0.5;
end
%Assigning y1 through y4 the corrisponding element of the vector
y1 = y_values(:,1); y2 = y_values(:,2); 
y3 = y_values(:,3); y4 = y_values(:,4);

% now to calculate x1 throuugh x4
x_values = [10000-y_values(:,1), y_values(:,2)-10000000, 25000000-y_values(:,3), 50000000-y_values(:,4)];
x1 = x_values(:,1); x2 = x_values(:,2); 
x3 = x_values(:,3); x4 = x_values(:,4);
A9 = abs(x_values(:,1)); A10 = abs(x_values(:,2));
A11 = abs(x_values(:,3)); A12 = abs(x_values(:,4));



















