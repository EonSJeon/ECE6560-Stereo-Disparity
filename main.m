clc; clear; close all;

I_paths = { ...
    './tsukuba/scene1.row3.col1.ppm', ...
    './tsukuba/scene1.row3.col2.ppm', ...
    './tsukuba/scene1.row3.col3.ppm', ...
    './tsukuba/scene1.row3.col4.ppm', ...
    './tsukuba/scene1.row3.col5.ppm'  ...
};

% Experiment settings
lambdas     = [0.95, 0.9, 0.6, 0.3];   
kernelSizes = [5, 3, 1];               
baseDir     = fullfile(pwd, 'results');
if ~exist(baseDir, 'dir'), mkdir(baseDir); end

% Define pairs to run 
pairs = [3 5; 1 5; 4 5];  

for p = 1:size(pairs,1)
    leftIdx  = pairs(p,1);
    rightIdx = pairs(p,2);
    IL_path  = I_paths{leftIdx};
    IR_path  = I_paths{rightIdx};
    
    % create subfolder for this pair
    subfolder = fullfile(baseDir, sprintf('L%02d_R%02d', leftIdx, rightIdx));
    if ~exist(subfolder, 'dir'), mkdir(subfolder); end

    % loop over kernels and lambdas
    for k = kernelSizes
        kernel = ones(1, k);
        for lam = lambdas
            % run disparity experiment
            exp = runExp(IL_path, IR_path, lam, kernel);

            % save experiment struct
            fname = sprintf('exp_lambda%.2f_k%dx1.mat', exp.lambda, k);
            save(fullfile(subfolder, fname), 'exp');

            % plot and save results into subfolder
            plotDisparityResults(exp, subfolder, true);
        end
    end
end


function exp = runExp(IL_path, IR_path, lambda, kernel)
    disp('==== Experiment Start ====')
    fprintf('IL path: %s\n',IL_path);
    fprintf('IR path: %s\n',IR_path);
    fprintf('lambda: %f\n', lambda);
    fprintf('kernel:\n');
    disp(kernel);

    % Default kernel if not provided
    if nargin < 4 || isempty(kernel)
        kernel = ones(1,3);
    end

    % Read and normalize
    IL = im2double(imread(IL_path));
    IR = im2double(imread(IR_path));

    % Dimensions and grid
    [h,w,~] = size(IL);
    [X,Y]   = meshgrid(1:w, 1:h);

    % Run disparity estimation
    disp('======= GrayScale Mode ========')
    [d_gray, E_gray]   = depthMapPatch(IL, IR, lambda, 'grayscale', kernel);
    disp('======= Color Mode ========')
    [d_color, E_color] = depthMapPatch(IL, IR, lambda, 'color',     kernel);

    % Compute final warped images
    IR_warp_gray = interp2(rgb2gray(IR), X - d_gray, Y, 'linear', 0);
    IR_warp_color = zeros(h,w,3);
    for c = 1:3
        IR_warp_color(:,:,c) = interp2(IR(:,:,c), X - d_color, Y, 'linear', 0);
    end

    % Package results in struct
    exp.IL_path       = IL_path;
    exp.IR_path       = IR_path;
    exp.lambda        = lambda;
    exp.kernel        = kernel;
    exp.d_gray        = d_gray;
    exp.d_color       = d_color;
    exp.E_gray        = E_gray;
    exp.E_color       = E_color;
    exp.IR_warp_gray  = IR_warp_gray;
    exp.IR_warp_color = IR_warp_color;
    disp('==== Experiment End ====')
end

function plotDisparityResults(exp, outputDir, save)
    if nargin < 2 || isempty(outputDir)
        outputDir = pwd;
    end
    % Read left image for reference
    IL = im2double(imread(exp.IL_path));

    % Unpack results
    d_gray        = exp.d_gray;
    d_color       = exp.d_color;
    IR_warp_gray  = exp.IR_warp_gray;
    IR_warp_color = exp.IR_warp_color;
    E_gray        = exp.E_gray;
    E_color       = exp.E_color;

    %% Figure 1: Disparity & Warps
    f1 = figure('Name','Disparity & Warps','NumberTitle','off','Position',[100 100 800 600]);
    subplot(3,2,1);
    imshow(d_color, []); title('Color Disparity'); colormap jet; colorbar;
    subplot(3,2,3);
    imshow(IL, []); title('Color IL (Reference)');
    subplot(3,2,5);
    imshow(IR_warp_color, []); title('Color IR (Warped)');
    subplot(3,2,2);
    imshow(d_gray, []); title('Grayscale Disparity'); colormap jet; colorbar;
    subplot(3,2,4);
    imshow(rgb2gray(IL), []); title('Grayscale IL (Reference)');
    subplot(3,2,6);
    imshow(IR_warp_gray, []); title('Grayscale IR (Warped)');
    drawnow;

    % Save figure if requested
    if save
        try
        saveas(f1, fullfile(outputDir, sprintf('DisparityAndWarps_lambda_%.2f.png', exp.lambda)));
        catch
            warning('Could not save Disparity & Warps figure.');
        end
    end
    
    %% Figure 2: Energy Loss History
    f2 = figure('Name','Energy Loss History','NumberTitle','off','Position',[200 200 600 400]);
    plot(E_color, 'LineWidth',2); hold on;
    plot(E_gray,  '--','LineWidth',2); hold off;
    xlabel('Iteration'); ylabel('Energy');
    title('Energy Loss: Color vs. Grayscale');
    legend({"Color","Grayscale"}, 'Location','northeast');
    grid on;
    drawnow;

    % Save energy plot
    if save
        try
            saveas(f2, fullfile(outputDir, sprintf('EnergyLoss_lambda_%.2f.png', exp.lambda)));
        catch
            warning('Could not save Energy Loss figure.');
        end
    end
    
end

