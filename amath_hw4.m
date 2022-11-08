% Tanner Huck    
% Math 301 B
% Homework 4

clear; clc; close all;

%% Question 1
% Loading in the data
load("particle_position.mat");

% Calculating the mean of everyrow
M = mean(A,2);

% For loop that will go through each row of matrix A and subtract the mean
% value for the row
for i = 1:3 
    A(i,:) = A(i,:) - M(i,:); 
end

% Part a
% Calculating the svd of our matrix A
[U, S, V] = svd(A, 'econ');

% Column vector with all the singular values
A1 = [S(1,1); S(2,2); S(3,3)];

% Part b
% Calculating the error of the rank 1 approximation of A
Arank1 = U(:,1)*S(1,1)*V(:,1)';
A2 = norm(A-Arank1);

% Part c
% Calculating the error of the rank 2 approximation of A
Arank2 = Arank1 + U(:,2)*S(2,2)*V(:,2)';
A3 = norm(A-Arank2);

%% Write up Question 1
% Part a.
% Plotting the particles position in 3d
plot3(A(1,:),A(2,:),A(3,:), 'k.');
hold on
% Plotting the rank 1 approximation of the particles position
plot3(Arank1(1,:),Arank1(2,:),Arank1(3,:), 'r.');
% Adding labels, title, and legend to the plot
title('Paticle Position and rank 1 aproximation')
xlabel('x');
ylabel('y');
zlabel('z');
legend('Particle Position','Rank 1 Approximation','location', 'best');
% Changing the font size for all labels
set(gca, 'FontSize', 12);

% Part b.
% Plotting the particles position in 3d
figure;
plot3(A(1,:),A(2,:),A(3,:), 'k.');
hold on
% Plotting the rank 2 approximation of the particles position
plot3(Arank2(1,:),Arank2(2,:),Arank2(3,:), 'r.');
% Adding labels, title, and legend to the plot
title('Paticle Position and rank 2 aproximation')
xlabel('x');
ylabel('y');
zlabel('z');
legend('Particle Position','Rank 1 Approximation','location', 'best');
% Changing the font size for all labels
set(gca, 'FontSize', 12);

%% Question 2
% Taking the olive image and turing it into gray scale
B = imread('olive.jpg');
B = double(rgb2gray(B));

% Part a
% Calculating the fifteen largest singular values of the image
S_B = svd(B, 'econ');
A4 = maxk(S_B,15);

% Part b
% Calculating the energy of the reank 1 approximation, the largest singular
% value divided by the sum of the singular values
sum = 0;
for i = 1:length(S_B)
    sum = sum + S_B(i,:);
end
A5 = A4(1,1) / sum; 

% Part c
% Calculating the energy of the rank 15 approximation, sum of the top 15
% singular values divided by the sum of all singular values
sum_top15 = 0;
for i = 1:length(A4)
    sum_top15 = sum_top15 + A4(i,:);
end
A6 = sum_top15 / sum;

% Part d
% calculating smallest r value such that energy of rank approximation is
% greater than or equal to 0.92
for r = 1:999
    r_amount = maxk(S_B,r);
    r_sum = 0;
    for r2 = 1:length(r_amount)
        r_sum = r_sum + r_amount(r2,:);
    end
    r_aprox = r_sum / sum;
    if r_aprox >= 0.92
        break
    end
end
A7 = r;

%% Write up Question 2
% Part a
% Ploting the original pictur in grayscale
% creating a subplot for 4 plots
subplot(2,2,1)
imagesc(B)
title('Original')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])
% Plotting the rank-1 aprox in grayscale
% Calculating the svd and rank 1 of the image
[Ub, Sb, Vb] = svd(B, 'econ');
Brank1 = Ub(:,1)*Sb(1,1)*Vb(:,1)';
% Ploting in the subplot
subplot(2,2,2)
imagesc(Brank1)
title('Rank-1 Approximation')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])
% Plotting the rank-10 aprox in grayscale
% Calculating the rank 10 aprox of the image
Brank10 = Ub(:,1:10)*Sb(1:10,1:10)*Vb(:,1:10)';
% plotting it in the subplot
subplot(2,2,3)
imagesc(Brank10)
title('Rank-10 Approximation')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])
% Plotting the rank-629 aprox in grayscale
% Calculating the rank 629 rank aprox of the image
Brank629 = Ub(:,1:629)*Sb(1:629,1:629)*Vb(:,1:629)';
% Plotting it in the subplot
subplot(2,2,4)
imagesc(Brank629)
title('Rank-629 Approximation')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])

% Part b
% Calculating the number of pixels in the image
totalPixels = 4032 * 3024;
% Calculating how many values required to store rank 629 approximation
B_approx_pixel_count = 4032*629 + 629 + 3024*629;

%% Question 3
% loading in the noisy image
load("NoisyImage.mat");

% Part a
% calculating the energy of the rank 2 approximation of the noisy image
singleVal_noise = svd(A_noise, 'econ');
% Total sum of all singular values
sumSingleVal_noise = 0;
for p = 1:length(singleVal_noise)
    sumSingleVal_noise = sumSingleVal_noise + singleVal_noise(p,:);
end
% Calculating the sum of the top 2 single values
topTwoSingleVal_noise = maxk(singleVal_noise,2);
sum_topTwoSingleVal_noise = topTwoSingleVal_noise(1,:) + topTwoSingleVal_noise(2,:);
% Calculating energy of rank 2 aprox
A8 = sum_topTwoSingleVal_noise / sumSingleVal_noise; 

% Part b
% calculating the rank two aprox of noise
[U2, S2, V2] = svd(A_noise, 'econ');
noise_rank1 = U2(:,1)*S2(1,1)*V2(:,1)';
noise_rank2 = noise_rank1 + U2(:,2)*S2(2,2)*V2(:,2)';
% error between A_noise and its rank 2 aproximation
A9 = norm(A_noise - noise_rank2);

% error between A and the A_noise rank 2 aproximation
A10 = norm(A - noise_rank2);

%% Write Up Question 3
% Part a
% Plotting the singular values of A_noise
figure
semilogy(singleVal_noise,'bo')
% setting y to be a log scale
title('Singular values of A-noise')

% Part b
figure
% creating a subplot for 3 plots
subplot(1,3,1)
% Ploting the original image
imagesc(A)
title('True Image')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])
% Ploting the noisy image
subplot(1,3,2)
imagesc(A_noise)
title('Noisy Image')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])
% Ploting the rank 2 aprroximation of the noisy image
subplot(1,3,3)
imagesc(noise_rank2)
title('Rank-2 Approx')
colormap gray
% Remoiving axis number
set(gca,'XTick',[], 'YTick', [])







