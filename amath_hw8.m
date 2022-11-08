% Tanner Huck    
% Math 301 B
% Homework 8

clear; clc; close all;

%% Coding Question 1
% Setting up the ODE
% dp/dt
dpdt = @(t, p) p*(1-p)*(p-(1/2));

% Forward Euler Method
[tans, pans] = forward_euler(dpdt, 0:1:15, 0.75);
A1 = pans;

% Backward Euler Method
[tans, pans] = backward_euler(dpdt, 0:1:15, 0.75);
A2 = pans;

% Midpoint method
[tans, pans] = midpoint(dpdt, 0:1:15, 0.75);
A3 = pans;

% RK4 method
[tans, pans] = RK4(dpdt, 0:1:15, 0.75);
A4 = pans;

%% Question 2
% part a
a = 1/sqrt(2);

ddx = [a, 1; -1, -a];
x = [1, 0; 0, 1];

A5 = x\ddx;

% part b
[tans, yans] = forward_euler_system(0:0.01:20, [2;1], A5);
A6 = yans(:, end);

% part c 
[tans, yans] = backward_euler_systems(A5, 0:0.01:20, [2;1]);
A7 = yans(:, end);



%% Question 3
% part a
g = 9.8;
L = 21;
sig = 0.06;
dtheta_dt = @(theta, v) v;
dv_dt = @(theta, v) (-g/L)*sin(theta) - sig*v;

% part b,c
odefun = @(theta, v) [dtheta_dt(v(1), v(2)); dv_dt(v(1), v(2))];
[tsol, ysol] = ode45(odefun, [0, 50], [pi/8, -0.1]);
A8 = ysol(end, :);
A8 = transpose(A8);

% part d
[tsol, ysol] = forward_euler_pendulum(0:0.01:50, [pi/8, -0.1]);
A9 = ysol(:, end);
A9 = abs(A8-A9);

%% Write up questions
% part a
% using meshgrid to generate a grid of points to use in a plot
thetha_bounds = linspace(-3*pi,3*pi,25); 
v_bounds = linspace(-3,3,25);
[x,y] = meshgrid(thetha_bounds,v_bounds);

% part b
% parameters for the pendulum problem
g = 9.8;
L = 21;
sig = 0.15;
% theta' and v' 
dtheta_dt = @(theta, v) v;
dv_dt = @(theta, v) (-g/L)*sin(theta) - sig*v;
% creating a grid of arrows with components (theta', ̇v')
quiver(x,y,dtheta_dt(x,y),dv_dt(x,y)); 

% part c
% adding labels and a title to the plot
title('Nonlinear Damped Pendulum', 'Fontsize', 20)
xlabel('\theta','Fontsize',15)
ylabel('velocity','Fontsize',15)

% part d
% solving the system with different initial conditions
% defining the ODE
odefun = @(theta, v) [dtheta_dt(v(1), v(2)); dv_dt(v(1), v(2))];
hold on
% initial condition theta = pie, v = 0.1
% solving the system and adding it to the plot
[tsol, ysol] = ode45(odefun, (0: 0.01: 50), [pi, 0.1]); 
plot(ysol(:, 1), ysol(:, 2), 'r', 'linewidth', 2)
hold on
% initial condition theta = pie, v = -0.1
% solving the system and adding it to the plot
[tsol, ysol] = ode45(odefun, (0: 0.01: 50), [pi, -0.1]); 
plot(ysol(:, 1), ysol(:, 2), 'b', 'linewidth', 2)
hold on
% initial condition theta = pie, v = 0.1
% solving the system and adding it to the plot
[tsol, ysol] = ode45(odefun, (0: 0.01: 50), [2*pi, -3]); 
plot(ysol(:, 1), ysol(:, 2), 'y', 'linewidth', 2)
hold on
% initial condition theta = pie, v = 0.1
% solving the system and adding it to the plot
[tsol, ysol] = ode45(odefun, (0: 0.01: 50), [-2*pi, 3]); 
plot(ysol(:, 1), ysol(:, 2), 'g', 'linewidth', 2)

% part e
% seting the axes to display only part of the plot
xlim([-3*pi,3*pi]);
ylim([-3,3]);



%% Functions below
% Forward Euler
function [t, y] = forward_euler(odefun,tspan,y0)
    % Forward Euler method
    % Solves the differential equation y' = f(t,y) at the times
    % specified by the vector tspan and with initial condition y0.
    %  - odefun is an anonymous function of the form odefun = @(t, v) ...
    %  - tspan is a row or column vector
    %  - y0 is a number

    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(length(tspan), 1); % Setup our solution column vector
    y(1) = y0; % Define the initial condition
    for k = 1:length(y)-1
        y(k+1) = y(k) + dt*odefun(tspan(k), y(k)); % Forward Euler step
    end
    t = tspan;
end

% Backward Euler
function [t, y] = backward_euler(odefun,tspan,y0)
    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(length(tspan), 1); % Setup our solution column vector
    y(1) = y0; % Define the initial condition
    for k = 1:length(y)-1
        g = @(z) z - y(k) - dt*odefun(tspan(k+1), z); % z = y(k+1)
        y(k+1) = fzero(g, y(k));
    end
    t = tspan;
end

% Midpoint 
function [t, y] = midpoint(odefun,tspan,y0)
    % Midpoint method - explicit method
    % Solves the differential equation y' = f(t,y) at the times
    % specified by the vector tspan and with initial condition y0.
    %  - odefun is an anonymous function of the form odefun = @(t, v) ...
    %  - tspan is a row or column vector
    %  - y0 is a number

    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(length(tspan), 1); % Setup our solution column vector
    y(1) = y0; % Define the initial condition
    for k = 1:length(y)-1
        k1 = odefun(tspan(k), y(k));
        k2 = odefun(tspan(k) + dt/2, y(k) + dt/2*k1);
        y(k+1) = y(k) + dt*k2; % Midpoint step
    end
    t = tspan;
end

% RK4 by hand
function [t, y] = RK4(odefun,tspan,y0)
    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(length(tspan), 1); % Setup our solution column vector
    y(1) = y0; % Define the initial condition
    for k = 1:length(y)-1
        k1 = odefun(tspan(k), y(k)); % formula for rk4
        k2 = odefun(tspan(k) + dt/2, y(k) + dt/2*k1);
        k3 = odefun(tspan(k) + dt/2, y(k) + dt/2*k2);
        k4 = odefun(tspan(k) + dt/2, y(k) + dt*k3);
        y(k+1) = y(k) + (dt/6)*(k1+2*k2+2*k3+k4);
    end
    t = tspan;
end

% Forward Euler for systems
function [t, y] = forward_euler_system(tspan,y0, matrixA)
    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(2, length(tspan)); % Setup our solution column vector
    y(:, 1) = y0; % Define the initial condition
    for k = 1:length(y)-1
         y(:, k+1) = y(:, k) + dt*matrixA*(y(:, k)); % Forward Euler step
    end
    t = tspan;
end

% Backward Euler for systems
function [t, y] = backward_euler_systems(A5,tspan,y0)
    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(2, length(tspan)); % Setup our solution column vector
    y(:, 1) = y0; % Define the initial condition
    I2 = [1,0;0,1];
    A = (I2 - dt.*A5);
    [l, u, p] = lu(A);
    for k = 1:length(y)-1
        b = y(:, k);
        y = l\(p*b);
        y(:, k+1) = u\y;
    end
    t = tspan;
end

% Forward Euler for pendulum
function [t, y] = forward_euler_pendulum(tspan,y0)
    dt = tspan(2)-tspan(1); % Calculate dt from the t values
    y = zeros(2, length(tspan)); % Setup our solution column vector
    y(1,1) = y0(1); % Define the initial condition
    y(2,1) = y0(2);
    for k = 1:length(y)-1
         y(1, k+1) = y(1, k) + dt*(y(2, k)); % Forward Euler step
         y(2, k+1) = y(2, k) + dt*(-9.8/21*sin(y(1, k))-0.06*y(2,k));
    end
    t = tspan;
end

