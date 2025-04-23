function [d, energyHistory] = depthMap(IL, IR, lambda, mode, kernel)

    % default
    if nargin < 4 || isempty(mode),  mode   = 'color';   end
    if nargin < 5 || isempty(kernel), kernel = ones(1,1); end  
    maxIters = 150000;

    IL_color = im2double(IL);
    IR_color = im2double(IR);

    % grayscale for loss
    IL_gray = rgb2gray(IL_color);

    [h,w,C] = size(IL_color);
    area    = h*w;
    kernel_weight = sum(kernel(:));
    

    % initialize disparity and log
    [X,Y]   = meshgrid(1:w,1:h);
    d = zeros(h,w);
    energyHistory = zeros(maxIters,1);

    

    %% Precompute Gradients
    if strcmpi(mode,'color')
        IxR_color = zeros(h,w,C);
        for c=1:C
            IxR_color(:,:,c) = backwardDiffX(IR_color(:,:,c));
        end
    else
        IxR_gray = backwardDiffX(rgb2gray(IR_color));
    end

    % initial loss
    energyHistory(1) = computeGrayDataEnergy(IL_gray, rgb2gray(IR_color), lambda);
    oldE = energyHistory(1);

    % PDE parameters
    dx       = 1;
    dt_cfl   = 0.25 * dx^2 / (1 - lambda);
    alpha    = 1;
    alpha_th = 1e-3;
    
    for iter=2:maxIters
        %% Warping & Dataterm
        mask = ((X-d)>=1);
        if strcmpi(mode,'color')
            % Color Mode
            IR_color_warp = zeros(h,w,C);
            for c=1:C
                IR_color_warp(:,:,c) = interp2(IR_color(:,:,c), X-d, Y, 'linear', 0);
            end
            IxR_color_warp = zeros(h,w,C);
            for c=1:C
                IxR_color_warp(:,:,c) = interp2(IxR_color(:,:,c), X-d, Y, 'linear', 0);
            end

            mismatch = (IL_color - IR_color_warp) .* mask;
            grad_chan = mismatch .* IxR_color_warp;

            p = zeros(h,w);
            for c=1:C
                p = p + conv2(grad_chan(:,:,c),rnel, 'same');
            end

            gradData = -lambda * p / kernel_weight * (1/3); 
            % Multiplying 1/3 makes this comparable.
        else
            % GrayScale Mode
            IR_gray_warp = interp2(rgb2gray(IR_color), X-d, Y, 'linear', 0);
            IxR_gray_warp = interp2(IxR_gray,    X-d, Y, 'linear', 0);

            mismatch_gray = (IL_gray - IR_gray_warp) .* mask;
            grad_chan    = mismatch_gray .* IxR_gray_warp;
            p = conv2(grad_chan, kernel, 'same');

            gradData = -lambda * p / kernel_weight;
        end

        %% Smoothness term
        lap        = laplacianNeumannBoundary(d);
        smoothTerm = 0.5*(1-lambda) * sum((gradient(d)).^2,'all');
        gradSmooth = (1-lambda) * lap;

        %% Compute loss
        if strcmpi(mode,'color')
            IRw_forLoss = rgb2gray(IR_color_warp);
        else
            IRw_forLoss = IR_gray_warp;
        end

        grayDataLoss    = computeGrayDataEnergy(IL_gray, IRw_forLoss, lambda, mask);
        mask_weight = sum(mask,'all');
        energyHistory(iter) = area/mask_weight * grayDataLoss + smoothTerm;
        % area/mask_weight -> should normalize the data loss to be applied
        % to whole pixels.

        %% PDE update
        gradTotal = gradData + gradSmooth;
        gmax      = max(abs(gradTotal(:)));
        dt        = min(dt_cfl, alpha*dx / max(gmax,eps));
        d         = max(d + dt * gradTotal, 0);

        %% Logging
        % Adaptive alpha
        if mod(iter,1000)==0
            fprintf('Iter %d, E=%.6f\n', iter, energyHistory(iter));
        end
        if mod(iter,50)==0
            if energyHistory(iter) > oldE
                alpha = alpha/4;
                fprintf('alpha has been changed to %f at iter %d.\n', alpha, iter);
            end
            if alpha < alpha_th
                energyHistory = energyHistory(1:iter);
                return;
            end
            oldE = energyHistory(iter);
        end
    end
end

function E = computeGrayDataEnergy(IL_gray, IR_gray, lambda, mask)
    if nargin < 4 || isempty(mask),  mask   = ones(size(IL_gray));   end
    m = (IL_gray - IR_gray).*mask;
    s = m.^2;
    E = 0.5 * lambda * sum(s(:));
end

function Bx = backwardDiffX(I)
    [r,c] = size(I);
    Bx = zeros(r,c);
    Bx(:,2:c) = I(:,2:c) - I(:,1:c-1);
    Bx(:,1)   = Bx(:,2);
end

function L = laplacianNeumannBoundary(U)
    [r,c] = size(U);
    Uup    = [U(1,:);    U(1:r-1,:)];
    Udown  = [U(2:r,:);  U(r,:)];
    Uleft  = [U(:,1),    U(:,1:c-1)];
    Uright = [U(:,2:c),  U(:,c)];
    L      = (Uup + Udown + Uleft + Uright) - 4*U;
end
