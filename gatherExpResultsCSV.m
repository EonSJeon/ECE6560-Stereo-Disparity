function gatherExpResultsCSV(rootDir)
% GATHEREXPREULTSCSV  Recursively scans 'rootDir' for files named
%   'exp_lambda<val>_k<kernel>x1.mat', loads them, and writes a CSV table
%   summarizing: which left photo, lambda, kernel, final energies & iteration
%   counts for both E_gray & E_color.

  % Find all matching .mat files under rootDir
  files = dir(fullfile(rootDir, '**', 'exp_lambda*_k*x1.mat'));

  % Open CSV for writing
  csvFile = 'results_table.csv';
  fid = fopen(csvFile, 'w');
  if fid < 0
    error('Unable to open %s for writing.', csvFile);
  end

  % Write header including new LeftPhoto column
  fprintf(fid, ...
    'FileName,RelPath,LeftPhoto,Lambda,Kernel,E_gray_end,Iter_gray,E_color_end,Iter_color\n');

  for f = 1:numel(files)
    filePath = fullfile(files(f).folder, files(f).name);

    % --- Extract LeftPhoto index from folder name ---
    [~, folderName] = fileparts(files(f).folder);
    % Expect folderName like 'L01_R05'
    tk = regexp(folderName, '^L(\d+)_R\d+$', 'tokens', 'once');
    if isempty(tk)
      leftNum = NaN;
    else
      leftNum = str2double(tk{1});
    end

    % --- Parse lambda & kernel from filename ---
    [~, baseName] = fileparts(files(f).name);
    tk2 = regexp(baseName, '^exp_lambda([\d\.]+)_k(\d+)x1$', 'tokens', 'once');
    if isempty(tk2)
      continue;
    end
    lamVal    = str2double(tk2{1});    % e.g. 0.30
    kernelVal = str2double(tk2{2});    % e.g. 1

    % --- Load the experiment struct ---
    S = load(filePath);
    expData = S.exp;

    % --- Final energies & iteration counts ---
    Eg = expData.E_gray;   Ec = expData.E_color;
    Eg_end   = Eg(end);
    Ec_end   = Ec(end);
    it_g     = numel(Eg);
    it_c     = numel(Ec);

    % --- Relative path (for reference) ---
    relPath = strrep(filePath, [rootDir filesep], '');

    % --- Write one CSV row ---
    fprintf(fid, '%s,%s,%d,%.3f,%d,%.6g,%d,%.6g,%d\n', ...
      baseName, relPath, leftNum, lamVal, kernelVal, ...
      Eg_end, it_g, Ec_end, it_c);
  end

  fclose(fid);
  fprintf('Results written to %s\n', csvFile);
end
