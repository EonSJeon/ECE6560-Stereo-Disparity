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

pairs = [3 5; 1 5; 4 5];  

for p = 1:size(pairs,1)
    leftIdx  = pairs(p,1);
    rightIdx = pairs(p,2);
    IL_path  = I_paths{leftIdx};
    IR_path  = I_paths{rightIdx};
    
    subfolder = fullfile(baseDir, sprintf('L%02d_R%02d', leftIdx, rightIdx));
    if ~exist(subfolder, 'dir'), mkdir(subfolder); end

    for k = kernelSizes
        kernel = ones(1, k);
        for lam = lambdas
            % run experiments
            exp = runExp(IL_path, IR_path, lam, kernel);

            % save experiment struct
            fname = sprintf('exp_lambda%.2f_k%dx1.mat', exp.lambda, k);
            save(fullfile(subfolder, fname), 'exp');
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

    if nargin < 4 || isempty(kernel)
        kernel = ones(1,3);
    end

    IL = im2double(imread(IL_path));
    IR = im2double(imread(IR_path));

    [h,w,~] = size(IL);
    [X,Y]   = meshgrid(1:w, 1:h);

    % Run disparity matching
    disp('======= GrayScale Mode ========')
    [d_gray, E_gray]   = depthMapPatch(IL, IR, lambda, 'grayscale', kernel);
    disp('======= Color Mode ========')
    [d_color, E_color] = depthMapPatch(IL, IR, lambda, 'color',     kernel);

    % Final warped images
    IR_warp_gray = interp2(rgb2gray(IR), X - d_gray, Y, 'linear', 0);
    IR_warp_color = zeros(h,w,3);
    for c = 1:3
        IR_warp_color(:,:,c) = interp2(IR(:,:,c), X - d_color, Y, 'linear', 0);
    end

    % Exp struct
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