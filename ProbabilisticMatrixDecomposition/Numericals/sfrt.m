function O=sfrt(n,l)
   %% R matrix
   R=zeros(n,l);
   o=1:l;
   for c=1:l
       u=rand;
       i=ceil(size(o,2)*u);
       r=o(i);
       R(r,c)=1;
       o=[o(1:i-1), o(i+1:size(o,2))];
   end
   %% D matrix
   theta=2*pi*rand(1,n);
   D=exp(-1j*theta);
   D=diag(D);
   
   %% F matrix
   F=dftmtx(n);
   O=sqrt(n/l)*D*F*R;
end