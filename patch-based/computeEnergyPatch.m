function E_val = computeEnergyPatch(d, IL, IR, lambda, Sel, patchRadius)
% E(d) = Σ_i,j [ Sel(i,j)* (λ/2)* Σ_{(rr,cc) in N(i,j)}(IL(rr,cc)-IR(rr,cc + d(i,j)))^2]
%        + (1-λ)/2 * Σ ||∇d||^2

    [rows, cols] = size(IL);
    dataTerm = 0;

    for i = 1:rows
        for j = 1:cols

            if Sel(i,j) == 0
                continue;
            end

            disp_ij = d(i,j);

            rowStart = max(i - patchRadius, 1);
            rowEnd   = min(i + patchRadius, rows);
            colStart = max(j - patchRadius, 1);
            colEnd   = min(j + patchRadius, cols);

            sumPatch = 0;
            for rr = rowStart:rowEnd
                for cc = colStart:colEnd
                    valL = IL(rr, cc);
                    ccWarp = cc + disp_ij;
                    if ccWarp < 1 || ccWarp > cols
                        continue;
                    end
                    valR = IR(rr, round(ccWarp));
                    diffVal = (valL - valR);
                    sumPatch = sumPatch + diffVal^2;
                end
            end

            dataTerm = dataTerm + 0.5*lambda* sumPatch;
        end
    end

    % Smoothness term
    [dx, dy] = gradient(d);
    smoothTerm = 0.5 * (1 - lambda) * sum(dx(:).^2 + dy(:).^2);

    E_val = dataTerm + smoothTerm;
end
