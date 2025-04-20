function grad_data = computeGradPatch(d, IL, IR, IxR, lambda, Sel, patchRadius)
%COMPUTEGRADPATCH 
%   For each pixel (i,j), sums over local patch N(i,j), 
%   mismatch = IL(rr, cc) - IR(rr, cc + d(i,j)).
%   Then derivative wrt d is mismatch * -IxR(rr, cc + d(i,j)).
%   Multiply by λ*Sel(i,j).

    [rows, cols] = size(IL);
    grad_data = zeros(rows, cols);

    for i = 1:rows
        for j = 1:cols

            if Sel(i,j) == 0
                continue; 
            end

            disp_ij = d(i,j);
            sumPatchGrad = 0;

            % define patch boundaries
            rowStart = max(i - patchRadius, 1);
            rowEnd   = min(i + patchRadius, rows);
            colStart = max(j - patchRadius, 1);
            colEnd   = min(j + patchRadius, cols);

            for rr = rowStart:rowEnd
                for cc = colStart:colEnd

                    valL = IL(rr, cc); 
                    ccWarp = cc + disp_ij; 
                    
                    % check bounds
                    if ccWarp < 1 || ccWarp > cols
                        continue;
                    end
                    % nearest neighbor (or you can do linear interpolation)
                    valR = IR(rr, round(ccWarp));

                    mismatch = (valL - valR);

                    % derivative of IR w.r.t. x
                    % IxR is forward diff: IxR(rr, cc+1) ~ IR(rr, cc+1)-IR(rr, cc)
                    % but we want IxR(rr, cc+disp). We'll do nearest neighbor:
                    cci = round(ccWarp);
                    if cci < 1 || cci > cols
                        continue;
                    end
                    gradI = - IxR(rr, cci); 

                    sumPatchGrad = sumPatchGrad + mismatch * gradI;
                end
            end
            % multiply by λ * Sel(i,j)
            grad_data(i,j) = lambda * sumPatchGrad;
        end
    end
end
