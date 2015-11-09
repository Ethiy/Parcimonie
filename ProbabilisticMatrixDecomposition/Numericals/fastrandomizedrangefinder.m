function Q=fastrandomizedrangefinder(A,l)
    n=length(A);
    O=sfrt(n,l);
    Y=A*O;
    [Q ~]=qr(Y);
end