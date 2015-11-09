function Wtild=weightmatrix(image,sigma,cropsize)
    imagesize=size(image,1);
    Wtild=zeros(imagesize^2);
    [X, I]=patchvector(image,cropsize);
    for i=1:imagesize^2
        for j=1:imagesize^2
            Wtild(i,j)=exp(-norm2(X(i,:)-X(j,:))^2/sigma^2);
        end
    end
end