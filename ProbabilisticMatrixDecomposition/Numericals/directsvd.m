function [S,U,V]=directsvd(A,Q)
   B=Q'*A;
   [U,S,V] = svd(B);
   U=Q*U;
end