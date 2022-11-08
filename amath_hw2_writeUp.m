% Tanner Huck    
% Math 301 B
% 4/14/2002;
% Homework 2 Writeup

clear; clc; close all;

%% Question 1
% part a
% The root finding problem that we need to solve to find the two values 
% that the population oscillates between is 
% 0 = x.*(exp(6-((3.*(x.*(1+exp(3.*(1-x))))))))-x.

% part b
% Ploting the function and seeing approximately where the roots are to make
% guesses with the fzero funciton
f = @(x) x.*(exp(6-((3.*(x.*(1+exp(3.*(1-x))))))))-x;
x = linspace(-5, 5, 10000);
plot(x, f(x), "LineWidth", 2)
ylim([-5,5])
title('Caterpillar population function')
xlabel('Population of Catepillars (thousands)') 
ylabel('Difference of equations')
% From the plot we see that the roots are around 0.5 and 2

% part c 
% first root x_1 using fzero
x_1_exact = fzero(f, 0.5);

% part d
% first finding  root x_2 using fzero
x_2_exact = fzero(f, 2);

% part e
% Looking at the graph to determine good intervals for a bisection method
% to find the roots x_1 and x_2. An interval for x_1 would be [0.1,0.5] and
% an interval for x_2 would be [1,2]

% part f
% Calculating the roots x_1 and x_2 using the bisection method
% using bisection method to find the root x_1, how many steps it took to
% get there, and creating a vector of each calulated midpoint
[x_1_biValues, x_1_biIterations, x_1_biMidpointsVector] = bisection(...
    @(x) x.*(exp(6-((3.*(x.*(1+exp(3.*(1-x))))))))-x, 0.1, 0.5, 10.^-12);
% Fixing the midpoint vector to the correct number of elements, how many
% steps we acually took
x_1_biMidpointsVector = x_1_biMidpointsVector(1:x_1_biIterations);
% Creating a vector with the error of each midpoint compared to the exact
% value
x_1_biErrorVector = abs(x_1_exact - x_1_biMidpointsVector);

% Now doing all the same things for root x_2
[x_2_biValues, x_2_biIterations, x_2_biMidpointsVector] = bisection(...
    @(x) x.*(exp(6-((3.*(x.*(1+exp(3.*(1-x))))))))-x, 1, 2, 10.^-12);
x_2_biMidpointsVector = x_2_biMidpointsVector(1:x_2_biIterations);
x_2_biErrorVector = abs(x_2_exact - x_2_biMidpointsVector);

% part g
% Now if we want to find the roots again using newtons method instead of
% bisection, looking at our grah, we can use initial guesses of 0.1 and 2
% for roots x_1 and x_2 respectivily. 

% part h 
% Using newtons method to find the root x_1
% The derivative of the function
df = @(x) -1+exp(6-3.*x.*(1+exp(3-3.*x)))-(3.*exp(6-3.*x.*(1+exp(3-3.*x))) ...
    .*x.*(-3.*x.*exp(-3.*x+3)+exp(-3.*x+3)+1));
% Tolerance used for the newtons method function
tol = 10.^-12; 
% Initial guess
x_1_nwFinalVal = 0.1; 
% Newtons method in a for loop
for i = 1:100
    if (abs(f(x_1_nwFinalVal)) >= tol) 
        x_1_nwErrorVector(i) = abs(x_1_exact - x_1_nwFinalVal); % Finds the
        % error of each guess
        % Calculates each new guess
        x_1_nwFinalVal = x_1_nwFinalVal - (f(x_1_nwFinalVal)/(df(x_1_nwFinalVal))); 
        x_1_nwGuess(i) = x_1_nwFinalVal; % Saves each guess
    else
        % Saves the number of steps
        step = i;
        % Includes the last guess in the error vector
        x_1_nwErrorVector(i) = abs(x_1_exact - x_1_nwFinalVal);
        break
    end
end

% Save newtons method but for the x_2 root
x_2_nwFinalVal = 2; 
for j = 1:100
    if (abs(f(x_2_nwFinalVal)) >= tol) 
        x_2_nwErrorVector(j) = abs(x_2_exact - x_2_nwFinalVal); 
        x_2_nwFinalVal = x_2_nwFinalVal - (f(x_2_nwFinalVal)/(df(x_2_nwFinalVal))); 
        x_2_nwGuess(j) = x_2_nwFinalVal; 
    else
        step = j;
        x_2_nwErrorVector(j) = abs(x_2_exact - x_2_nwFinalVal);
        break
    end
end



%% Ploting
figure; % Figure 2
% Bisection method for x_1
plot((1:x_1_biIterations), x_1_biErrorVector,'k+') 
hold on
% Newtons method for x_1
plot(0:length(x_1_nwGuess), x_1_nwErrorVector, 'g*')
% Adding lables, title, changing font, and creating a legend
xlabel('number of iterations”') 
ylabel('error')
title('Error for finding x_1')
set(gca, 'FontSize', 15)
legend('Bisection Method','Newtons Method','Location', 'NorthEast')

figure; % Figure 3
% Same code, but for x_2 instead of x_1
plot((1:x_2_biIterations), x_2_biErrorVector,'k+')
hold on
plot(0:length(x_2_nwGuess), x_2_nwErrorVector, 'g*')
xlabel('number of iterations”') 
ylabel('error')
title('Error for finding x_2')
set(gca, 'FontSize', 15)
legend('Bisection Method','Newtons Method','Location', 'NorthEast')


%% Question 2
figure; % Figure 4
% Same code as figures 2 and 3, but with semilogy instead of plot
semilogy((1:x_1_biIterations), x_1_biErrorVector,'k+') 
hold on
semilogy(0:length(x_1_nwGuess), x_1_nwErrorVector, 'g*')
xlabel('number of iterations”') 
ylabel('error')
title('Error for finding x_1 (logarithmic scale)')
set(gca, 'FontSize', 15)
legend('Bisection Method','Newtons Method','Location', 'NorthEast')

figure; % Figure 6
% Same code as figure 5 but with x_2
semilogy((1:x_2_biIterations), x_2_biErrorVector,'k+')
hold on
semilogy(0:length(x_2_nwGuess), x_2_nwErrorVector, 'g*')
xlabel('number of iterations”') 
ylabel('error')
title('Error for finding x_2 (logarithmic scale)')
set(gca, 'FontSize', 15)
legend('Bisection Method','Newtons Method','Location', 'NorthEast')

% part a
% I believe that newtons method is better preforming for this equation. I
% believe this because it reaches the root more efficiently than
% the bisection method. As we can see from the graph, it only takes newtons
% method around 5 or 6 steps to rech the root, whereas it takes the
% bisection 10+ steps to reach the root. I have also found that the newtons
% method will run much faster on matlab. When I compare how quickly matlab
% can run bisection method vs newtons method, newtons method will complete
% running a noticable amount faster. Overall, newtons method is better
% preforming than the bisection method beause it is more efficient and
% faster at finding roots.





%% Functions
% Function for bisections
function [root_value, iterations, midpointsVector] = bisection(fun, left, ...
    right, tol)
    midpointsVector = (1:1000);
    for k = 1:1000
        mid = (left + right)/2; % Finds the midpoint between two points 
        midpointsVector(k) = mid; 
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











