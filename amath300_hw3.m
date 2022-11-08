% Tanner Huck    
% Math 301 B
% Homework 3

clear; clc;

%% Question 1
% part a
% Using the rotationMatrix function to create a rotation matrix that will
% rotate 2pi/5 degrees around the z axis
A1 = rotationMatrix(2*pi/5);

% part b
% Finding the vector x if x multiplied by the rotaion matrix is eual to the
% vector [pi; 3; 4]. We can do this by using the "\" in matlab
vectorY = [pi; 3; 4];
A2 = A1\vectorY;

% part c
% Finding the vector x again, but using lu decomposition instead
[L1, U1, P1] = lu(A1);
y_LU = L1\(P1*vectorY);
x_LU = U1\y_LU;
A3 = x_LU;

%% Question 2
% part a&b
% Writing a series of linear equations Ax=b
% The matrix A in our equation
A = [-1/sqrt(17),1,0,0,0,0,0,0,0,0,0,1/sqrt(17),0,0,0,0,0,0,0;
    -4/sqrt(17),0,0,0,0,0,0,0,0,0,0,-4/sqrt(17),0,0,0,0,0,0,0;
    0,-1,1,0,0,0,0,0,0,0,0,0,-1/sqrt(17),1/sqrt(17),0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,-4/sqrt(17),-4/sqrt(17),0,0,0,0,0;
    0,0,-1,1,0,0,0,0,0,0,0,0,0,0,-1/sqrt(17),1/sqrt(17),0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4/sqrt(17),-4/sqrt(17),0,0,0;
    0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,-1/sqrt(17),1/sqrt(17),0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4/sqrt(17),-4/sqrt(17),0;
    0,0,0,0,-1,1/sqrt(17),0,0,0,0,0,0,0,0,0,0,0,0,-1/sqrt(17);
    0,0,0,0,0,-4/sqrt(17),0,0,0,0,0,0,0,0,0,0,0,0,-4/sqrt(17);
    0,0,0,0,0,-1/sqrt(17),-1,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,-1/sqrt(17),1/sqrt(17);
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4/sqrt(17),4/sqrt(17);
    0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,-1/sqrt(17),1/sqrt(17),0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4/sqrt(17),4/sqrt(17),0,0;
    0,0,0,0,0,0,0,0,1,-1,0,0,0,-1/sqrt(17),1/sqrt(17),0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,4/sqrt(17),4/sqrt(17),0,0,0,0;
    0,0,0,0,0,0,0,0,0,1,-1,-1/sqrt(17),1/sqrt(17),0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,4/sqrt(17),4/sqrt(17),0,0,0,0,0,0];
% The vecotor b in our equation
b = zeros(19,1);
% Setting the 8th,9th,10th, and 11th values to 10000,9700,12000,21500
b(13,1)=10000; b(15,1)=9700; b(17,1)=12000; b(19,1)=21500;

% part c
% Setting our matrix to A4
A4 = A;

% part d
% Solving our system for the vector x using backlash
x = A4\b;
A5 = x;

% part e
% Creating L,U, and P from matrix A
[L2, U2, P2] = lu(A);
A6 = L2;

% part f
% Finding the largest value in vector x
A7 = max(abs(x));

% part g
maxval = -69;
while maxval <= 42000
    b(17,1) = b(17,1) + 5;
    y_LU = L2\(P2*b);
    x_LU = U2\y_LU;
    [maxval, index] = max(abs(x_LU));
end
A8 = b(17,1);
A9 = index;





%% Functions
% rotaitonMatrix function will take an input of an angle and create a
% rotation matrix that will rotate around the z axis by that angle
function [rotationMatrix] = rotationMatrix(theta)
    rotationMatrix = [cos(theta), -sin(theta), 0; ...
            sin(theta),cos(theta), 0; ...
            0, 0, 1];
end









