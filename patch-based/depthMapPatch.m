function [d, energyHistory] = depthMapPatch( ...
    IL, IR, lambda, numIters, stepSize, maskRatio, patchRadius)
%DEPTHMAPPATCH Variational stereo with a patch-based data term:
%
%   E(d) = ∑_x [ Sel(x) * (λ/2) * ∑_{x' in N(x)}(IL(x') - IR(x'+d(x)))^2 ]
%           + ∑_x [ (1-λ)/2 * ||∇d(x)||^2 ],
%
% solved by explicit gradient descent + bounding 0 <= d <= 0.1*width.

    [rows, cols] = size(IL);

    % 1) Create mask: 0 in rightmost fraction, 1 elsewhere
    cutoffCol = round((1 - maskRatio) * cols);
    Sel = ones(rows, cols);
    if cutoffCol < cols
        Sel(:, cutoffCol+1:end) = 0;
    end

    % 2) Initialize disparity
    d = zeros(rows, cols);
    globalMax = 0.1 * cols;  % clamp upper bound

    % 3) Precompute partial derivative of IR wrt x (for data gradient)
    %    We'll do a simple forward difference. 
    IxR = forwardDiffX(IR);

    % 4) Energy history + initial energy
    energyHistory = zeros(numIters,1);
    energyHistory(1) = computeEnergyPatch(d, IL, IR, lambda, Sel, patchRadius);

    % Stop criteria
    increaseTolerance  = 1e-10;
    increasePatience   = 5;
    convergeTolerance  = 1e-4;
    convergePatience   = 10;
    incCount  = 0;
    convCount = 0;

    % ----------- Main Gradient Descent -----------
    for iter = 2:numIters

        % (A) Data Term Gradient (patch-based)
        grad_data = computeGradPatch(d, IL, IR, IxR, lambda, Sel, patchRadius);

        % (B) Smoothness Term Gradient (Neumann BC)
        lap_d = laplacianNeumann(d);
        grad_smooth = (1 - lambda) .* lap_d;

        % (C) Combine & Update
        grad_total = grad_data + grad_smooth;
        d = d - stepSize * grad_total;

        % (D) Enforce bounds: 0 <= d <= globalMax
        d = max(d, 0);
        d = min(d, globalMax);

        % (E) Compute new energy + check stop conditions
        energyHistory(iter) = computeEnergyPatch(d, IL, IR, lambda, Sel, patchRadius);

        if energyHistory(iter) > energyHistory(iter-1) + increaseTolerance
            incCount = incCount + 1;
        else
            incCount = 0;
        end

        if abs(energyHistory(iter) - energyHistory(iter-1)) < convergeTolerance
            convCount = convCount + 1;
        else
            convCount = 0;
        end

        if incCount >= increasePatience
            fprintf('Stopping early: energy increased %d times in a row.\n', incCount);
            energyHistory = energyHistory(1:iter);
            break;
        end
        if convCount >= convergePatience
            fprintf('Stopping early: energy converged for %d consecutive steps.\n', convCount);
            energyHistory = energyHistory(1:iter);
            break;
        end

        if mod(iter, 500) == 0
            fprintf('Iter %d/%d, E = %.6f\n', iter, numIters, energyHistory(iter));
        end
    end
end
