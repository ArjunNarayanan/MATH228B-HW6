function [V,inV,Dr] = VandermondeMatrix(P,r)
%Vandermonde Matrix, V_{ij} = phi_j(r_i)

V = zeros(length(r),P+1);
Vr = zeros(length(r),P+1);

for j=1:P+1
    V(:,j) = JacobiP(r(:), 0, 0, j-1);
end
inV = inv(V);

for i=0:P
   [Vr(:,i+1)] = GradJacobiP(r(:),0,0,i);
end

Dr = Vr/V;

return
