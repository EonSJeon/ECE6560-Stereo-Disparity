% === Load results ===
T = readtable('results_table.csv');

% Unique grid values
leftVals   = unique(T.LeftPhoto);   % e.g. [1;3;4]
kernVals   = unique(T.Kernel);      % e.g. [1;3;5]

% Preallocate matrices
nL = numel(leftVals);
nK = numel(kernVals);
IterGray  = NaN(nL,nK);
IterColor = NaN(nL,nK);
Egray_end = NaN(nL,nK);
Ecolor_end= NaN(nL,nK);

% Fill in matrices
for i = 1:nL
  for j = 1:nK
    mask = T.LeftPhoto==leftVals(i) & T.Kernel==kernVals(j);
    row  = T(mask,:);
    if isempty(row), continue; end
    IterGray(i,j)   = row.Iter_gray;
    IterColor(i,j)  = row.Iter_color;
    Egray_end(i,j)  = row.E_gray_end;
    Ecolor_end(i,j) = row.E_color_end;
  end
end

% Offsets for grouping
dx = 0.15;  

% --- 1) Iteration‐count bar3 ---
figure('Name','Iterations','NumberTitle','off')
hold on

% Grayscale bars (shift left)
hG = bar3(IterGray);
for k = 1:numel(hG)
    x = get(hG(k),'XData');
    set(hG(k),'XData', x - dx);
    set(hG(k),'FaceColor',[0.5 0.5 0.5]);  % gray
end

% Color bars (shift right)
hC = bar3(IterColor);
for k = 1:numel(hC)
    x = get(hC(k),'XData');
    set(hC(k),'XData', x + dx);
    set(hC(k),'FaceColor',[1 0 0]);        % red
end

% Axes labels & ticks
set(gca, ...
  'XTick', 1:nL, ...
  'XTickLabel', string(leftVals), ...
  'YTick', 1:nK, ...
  'YTickLabel', string(kernVals))
xlabel('LeftPhoto')
ylabel('Kernel')
zlabel('Iterations')
title('Convergence Iterations (Gray vs. Color)')
legend({'Grayscale','Color'}, 'Location','best')
view(45,30)
grid on
hold off

% --- 2) Final‐energy bar3 ---
figure('Name','Energies','NumberTitle','off')
hold on

hG = bar3(Egray_end);
for k = 1:numel(hG)
    x = get(hG(k),'XData');
    set(hG(k),'XData', x - dx);
    set(hG(k),'FaceColor',[0.5 0.5 0.5]);
end

hC = bar3(Ecolor_end);
for k = 1:numel(hC)
    x = get(hC(k),'XData');
    set(hC(k),'XData', x + dx);
    set(hC(k),'FaceColor',[1 0 0]);
end

set(gca, ...
  'XTick', 1:nL, ...
  'XTickLabel', string(leftVals), ...
  'YTick', 1:nK, ...
  'YTickLabel', string(kernVals))
xlabel('LeftPhoto')
ylabel('Kernel')
zlabel('Final Energy')
title('Final Energy (Gray vs. Color)')
legend({'Grayscale','Color'}, 'Location','best')
view(45,30)
grid on
hold off
