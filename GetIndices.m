function [I J]  = GetIndices(N)
I = zeros(N*(N-1)/2,1);
J = zeros(N*(N-1)/2,1);
k = 1;
for i = 1 : N
    for j = i+1 : N
        I(k) = i;
        J(k) = j;
        k = k + 1;
    end
end

        