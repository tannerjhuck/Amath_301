clear all; close all; clc

% This script reads in an image, shows it, but then converts it to a
% grayscale image so that we can use the SVD on it.

%% First read in the image
A = imread('lighthouse.png'); % Read in the image "lighthouse.png"
imshow(A) % plots the image

size(A) % See that this is a 640 x 480 x 3 array, so we can't do SVD on it.

%% Convert to gray and then redo it.
A = double(rgb2gray(A)); % Converts A to a grayscale image, 
                         % which is now a matrix.
                         
size(A) % See, it is a matrix!

% Now to show it:
imagesc(A)
colormap gray

%% Now that we have a matrix, let's take the SVD
[U, S, V] = svd(A, 'econ');

% I want to first look at the singular values of A, because they tell me
% the directions that carry the most variance or information.
sing_vals = diag(S); % S is a diagonal matrix, diag(S) gives an array from that diagonal

figure
plot(sing_vals, 'o')

semilogy(sing_vals, 'o')

%% Let's calculate the rank-1 approximation of this image
Arank1 = U(:,1)*S(1,1)*V(:,1)';

% Let's calculate how much "energy" this has
% We will get back to this
figure
subplot(1,2,1) % This creates a 1x2 dimensional array of plots (1 row, 2 columns)
                % Then the "1" in the last argument says "I am working with
                % the first thing in that 1x2 dimensional array"
imagesc(A) % True image
subplot(1,2,2) % Work with the 2nd element of that 1x2 array
imagesc(Arank1)
colormap gray

%% How much "energy" or "information" is in this rank-1 approx
rank1_energy = sing_vals(1) / sum(sing_vals)

%% What if we did a rank-10 approximation?
r = 10;
A_approx = U(:, 1:r)*S(1:r, 1:r)*V(:,1:r)'; %This has the sum 
                                            % built in
                                            % This is the rank-10
                                            % approximation
energy = sum(sing_vals(1:r))/sum(sing_vals)

figure
subplot(1,2,1) % 1 row, 2 columns, first image
imagesc(A)
subplot(1,2,2) % 1 row, 2 cols, second image
imagesc(A_approx)
colormap gray

%% Let's say we want to calculate 90% of the energy
% We could do this in a for loop, or we can use a built-in
% MATLAB function called cumsum
% cumsum calculates the cumulative sum of a vector
energies = cumsum(sing_vals)/sum(sing_vals)
figure
plot(energies)

%% We are going to use a rank-150 approximation
r = 150;
A_approx = U(:, 1:r)*S(1:r, 1:r)*V(:,1:r)'; %This has the sum 
                                            % built in
                                            % This is the rank-10
                                            % approximation
energy = sum(sing_vals(1:r))/sum(sing_vals)

figure
subplot(1,2,1) % 1 row, 2 columns, first image
imagesc(A)
subplot(1,2,2) % 1 row, 2 cols, second image
imagesc(A_approx)
colormap gray

%% I claim that we are capturing 90.34% of the "energy" but
% storing less data, let's check that.
size(A)
% How many pixels are stored in A
A_pixel_count = 640*480

% To recreate the rank-r approximation you need the formula:
% A_approx = U(:, 1:r)*S(1:r, 1:r)*V(:,1:r)';
size(U(:, 1:r))
size(S(1:r, 1:r))
size(V(:, 1:r))

A_approx_pixel_count = 640*150 + 150 + 480*150
% 150*(640 + 1 + 480)

% So we store
A_approx_pixel_count/A_pixel_count % The percentage of pixels
                                   % stored for the approx.
% So we are capturing 90% of the information/energy in the 
% picture using 54% of the data. 

