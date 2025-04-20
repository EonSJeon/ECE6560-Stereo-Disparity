% -------------------------------------------------------------------------
% Main Script
% -------------------------------------------------------------------------
clc; clear; close all;

% Load two stereo images (left, right)
IL = im2double(imread('./tsukuba/scene1.row3.col5.ppm'));
IR = im2double(imread('./tsukuba/scene1.row3.col4.ppm'));

% Convert to grayscale if needed
if size(IL,3) == 3, IL = rgb2gray(IL); end
if size(IR,3) == 3, IR = rgb2gray(IR); end

% Parameters
lambda       = 0.6;      % data vs. smoothness weight
numIters     = 100000;   % max iterations
fixedStep    = 1e-3;     % gradient descent step size
maskRatio    = 0.03;     % fraction of columns to mask on right side
patchRadius  = 1;        % e.g. 3x3 patch => radius=1

% Run disparity estimation with patch-based data term
[d_est, energyHistory] = depthMapPatch(...
    IL, IR, lambda, numIters, fixedStep, maskRatio, patchRadius);

% Show the estimated disparity
figure; 
imshow(d_est, []);
title('Disparity Map (Patch-Based Data Term)');
colormap jet; colorbar;

% Plot the loss (energy) history
figure;
plot(energyHistory, 'LineWidth', 2);
xlabel('Iteration');  ylabel('Energy');
title('Energy Loss History (Patch-Based)');
grid on;

% Approximate depth as 1 / (disparity + eps)
depthMap_approx = 1./(d_est + eps);
figure;
imshow(depthMap_approx, []);
title('Approx. Depth = 1/(Disparity) -- Patch-Based');
colormap jet; colorbar;
