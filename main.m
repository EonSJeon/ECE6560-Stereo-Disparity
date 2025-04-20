% -------------------------------------------------------------------------
% Main Script
% -------------------------------------------------------------------------
clc; clear; close all;

IL = im2double(imread('./tsukuba/scene1.row3.col1.ppm'));
IR = im2double(imread('./tsukuba/scene1.row3.col5.ppm'));
% IL = im2double(imread('./test1.png'));
% IR = im2double(imread('./test2.png'));


if size(IL,3) == 3, IL = rgb2gray(IL); end
if size(IR,3) == 3, IR = rgb2gray(IR); end

% Parameters
lambda       = 0.8;    
numIters     = 150000;         

[d_est, energyHistory] = depthMap(IL, IR, lambda, numIters);

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

function [d, energyHistory] = depthMap(IL, IR, lambda, numIters)
    dt_cfl = 0.25;
    dx=1;

    IL = double(IL);
    IR = double(IR);

    [h, w] = size(IL);

    % Initialize disparity
    d = ones(h, w)*0;
    % d =15;
    % size(d)

    % Central difference in x
    IxR = backwardDiffX(IR);

    % Make coordinate grid
    [X, Y] = meshgrid(1:w, 1:h);

    % Preallocate energy history
    energyHistory = zeros(numIters,1);

    % Compute initial energy
    init_mismatch = IL-IR;
    energyHistory(1) = computeEnergy(d, init_mismatch, lambda);

    % ------------- Gradient Descent Loop -------------
    for iter = 2:numIters
        % --- Data Term Gradient (masked) ---
        mask = (X - d)>=1;       
        IR_warp  = interp2(IR, X - d, Y, 'linear', 0);
        IxR_warp = interp2(IxR, X - d, Y, 'linear', 0);

        mismatch = (IL - IR_warp).*mask;
        
        % figure;
        % subplot(3,1,1);
        % imshow(IL);
        % subplot(3,1,2);
        % imshow(IR);
        % subplot(3,1,3);
        % imshow(IR_warp);
        % 
        % figure;
        % imshow(mask);

        grad_data = -lambda * mismatch .* IxR_warp .* mask;

        % figure;
        % imshow(grad_data~=0);

        % --- Smoothness Gradient  ---
        % size(d)
        lap_d = laplacianNeumann(d);

        % figure;
        % imshow(lap_d);
        grad_smooth = (1 - lambda) .* lap_d;

        % --- Combine & Update (unconstrained) ---
        grad_total = grad_data + grad_smooth;
        grad_max = max(abs(grad_total(:)));
        stepSize = min(dt_cfl, 0.5*dx/grad_max);
        % update = min(max(dt_cfl *grad_total,-dx),dx);
        d = d + grad_total * stepSize;
        % d = max(d,0);

        % --- Compute new energy ---
        energyHistory(iter) = computeEnergy(d, mismatch, lambda);

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


function Fx = forwardDiffX(I)
%FORWARDDIFFX Computes forward difference in the x-direction
%   Fx(i,j) = I(i,j+1) - I(i,j), except for the last column which is replicated.

    [rows, cols] = size(I);
    Fx = zeros(rows, cols);
    % For each row, difference in x
    Fx(:,1:end-1) = I(:,2:end) - I(:,1:end-1);
    % For the last column, replicate or set to 0
    Fx(:,end) = Fx(:,end-1);  
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: laplacianNeumann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function E_val = computeEnergy(d, mismatch, lambda)
%COMPUTEENERGYDATAMASKED
% E(d) =  ∑ [ Sel(x)*(λ/2*(IL - IR(x+d))^2 ] + ∑ [ (1-λ)/2*||∇d||^2 ].
    
    dataTerm = 0.5 * lambda * sum( mismatch(:).^2 );

    [dx, dy] = gradient(d);
    smoothTerm = 0.5 * (1 - lambda) * sum(dx(:).^2 + dy(:).^2);

    E_val = dataTerm + smoothTerm;
end