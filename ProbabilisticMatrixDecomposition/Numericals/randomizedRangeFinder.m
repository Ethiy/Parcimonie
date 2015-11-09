function Q=randomizedRangeFinder(A,l,q)
    O=random('Normal',0,1,size(A,2),l);
    Y=A*O;
    [Q ~]=qr(Y);
    for j=1:q
        Y=A'*Q;
        [Q ~]=qr(Y);
        Y=A*Q;
        [Q ~]=qr(Y);
    end
end