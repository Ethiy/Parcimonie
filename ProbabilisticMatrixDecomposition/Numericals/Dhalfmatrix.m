function D=Dhalfmatrix(W)
    D=sum(W);
    D=1./D;
    D=sqrt(D);
    D=diag(D);
end