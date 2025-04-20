function Fx = forwardDiffX(I)
%FORWARDDIFFX in x-direction, forward difference
    [rows, cols] = size(I);
    Fx = zeros(rows, cols);
    Fx(:,1:end-1) = I(:,2:end) - I(:,1:end-1);
    Fx(:,end) = Fx(:,end-1);  
end

