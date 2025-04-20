function L = laplacianNeumann(U)
%LAPLACIANNEUMANN 2D Laplacian with replicated boundaries (approx Neumann)
    [rows, cols] = size(U);

    U_up    = [U(1,:);  U(1:rows-1,:)];
    U_down  = [U(2:rows,:); U(rows,:)];
    U_left  = [U(:,1), U(:,1:cols-1)];
    U_right = [U(:,2:cols), U(:,cols)];

    L = (U_up + U_down + U_left + U_right) - 4*U;
end
