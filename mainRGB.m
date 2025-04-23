% -------------------------------------------------------------------------
% Main Script
% -------------------------------------------------------------------------
clc; clear; close all;

IL = im2double(imread('./tsukuba/scene1.row3.col4.ppm'));
IR = im2double(imread('./tsukuba/scene1.row3.col5.ppm'));
% IL = im2double(imread('./test1.png'));
% IR = im2double(imread('./test2.png'));

% Parameters
lambda       = 0.3;    
numIters     = 150000;         

[d_est, energyHistory] = depthMapColor(IL, IR, lambda, numIters);

% Estimated disparity
figure; 
imshow(d_est, []);
title('Disparity Map');
colormap jet; colorbar;

% Plot the loss history
figure;
plot(energyHistory, 'LineWidth', 2);
xlabel('Iteration');  ylabel('Energy');
title('Energy Loss History');
grid on;

% Approximate depth as 1 / (disparity + eps)
depthMap_approx = 1./(d_est + 0.01);
figure;
imshow(depthMap_approx, []);
title('Approx. Depth = 1/(Disparity)');
colormap jet; colorbar;



function [d, energyHistory] = depthMapColor(IL, IR, lambda, numIters)
    dt_cfl=0.25;
    dx=1;
    
    % IL, IR: H x W x 3 color images
    [h, w, ~] = size(IL);

    % Initialize disparity
    d = zeros(h, w);

    % Precompute derivative of IR in x for each channel
    IxR = zeros(h, w, 3);
    for c = 1:3
        IxR(:,:,c) = backwardDiffX(IR(:,:,c));  % or central difference
    end

    [X, Y] = meshgrid(1:w, 1:h);
    energyHistory = zeros(numIters,1);

    % PDE loop
    for iter = 1:numIters
    
        % 1) Warp IR each channel
        IR_warp = zeros(h, w, 3);
        IxR_warp = zeros(h, w, 3);

        mask = (X - d)>=1;    
        for c = 1:3
            IR_warp(:,:,c)  = interp2(IR(:,:,c), X - d, Y, 'linear', 0);
            IxR_warp(:,:,c) = interp2(IxR(:,:,c), X - d, Y, 'linear', 0);
        end

        % 2) Mismatch for each channel
        mismatchR = IL(:,:,1) - IR_warp(:,:,1).*mask;
        mismatchG = IL(:,:,2) - IR_warp(:,:,2).*mask;
        mismatchB = IL(:,:,3) - IR_warp(:,:,3).*mask;

        % sum partial derivatives => grad_data
        grad_dataR = mismatchR .* IxR_warp(:,:,1).*mask;
        grad_dataG = mismatchG .* IxR_warp(:,:,2).*mask;
        grad_dataB = mismatchB .* IxR_warp(:,:,3).*mask;

        grad_data = - lambda * ( grad_dataR + grad_dataG + grad_dataB ).*mask;

        % 3) Smoothness term (e.g. Neumann Laplacian)
        lap_d      = laplacianNeumann(d);
        grad_smooth= (1 - lambda) * lap_d;

        grad_total = grad_data + grad_smooth;

        % time-step
        grad_max = max(abs(grad_total(:)));
        stepSize = min(dt_cfl, 0.25*dx/grad_max);
        d = d + grad_total * stepSize;
        
        % --- Compute new energy ---
        energyHistory(iter) = computeColorEnergy(d, mismatchR,mismatchG, mismatchB, lambda);

        % --- Optional intermediate plots every 10000 iters ---
        if mod(iter,5000)==0
            h_iter = figure;
            subplot(4,1,1);
            imshow(d, []);
            title(sprintf('Disparity at iter %d', iter));
            colormap jet; colorbar;
            subplot(4,1,2);
            imshow(IL, []);
            title('IL reference');
            subplot(4,1,3);
            imshow(IR, []);
            title('IR reference');
            subplot(4,1,4);
            imshow(IR_warp);
            title('IR warp')
            saveas(h_iter, sprintf('iter_%05d.png', iter));
            close(h_iter);

            % figure;
            % imshow(mismatch);
            % title('Mismatch');
            % 
            % figure;
            % imshow(grad_data);
            % title('grad data')
            % ;
            
        end

        % Print progress
        if mod(iter, 1000) == 0
            fprintf('Iteration %d / %d, E = %.6f\n', iter, numIters, energyHistory(iter));
        end
    end
end

function Bx = backwardDiffX(I)
%BACKWARDDIFFX  Computes backward difference in x-direction
%   Bx(i,j) = I(i,j) - I(i,j-1), with the first column replicated.
    [rows, cols] = size(I);
    Bx = zeros(rows, cols);
    % For each row, backward difference in x
    Bx(:,2:cols) = I(:,2:cols) - I(:,1:cols-1);
    % For the first column, replicate the next value (or set to zero)
    Bx(:,1) = Bx(:,2);
end

function L = laplacianNeumann(U)
%LAPLACIANNEUMANN Computes the 2D Laplacian of U using a 5-point stencil
% with replicated boundaries to approximate Neumann boundary conditions.

    [h, w] = size(U);

    % Replicate boundaries
    U_up    = [U(1,:);    U(1:h-1,:)];   
    U_down  = [U(2:h,:); U(h,:)];     
    U_left  = [U(:,1),   U(:,1:w-1)];    
    U_right = [U(:,2:w), U(:,w)];

    % size(U_up)

    % Compute the Laplacian
    L = (U_up + U_down + U_left + U_right) - 4 * U;
    % size(L)
end

function E_val = computeColorEnergy(d, mismatchR,mismatchG, mismatchB, lambda)

    dataTerm = 0.5 * lambda * ...
        ( sum(mismatchR(:).^2) + sum(mismatchG(:).^2) + sum(mismatchB(:).^2) );

    [dx, dy] = gradient(d);
    smoothTerm = 0.5 * (1 - lambda) * sum(dx(:).^2 + dy(:).^2);

    E_val = dataTerm + smoothTerm;
end
