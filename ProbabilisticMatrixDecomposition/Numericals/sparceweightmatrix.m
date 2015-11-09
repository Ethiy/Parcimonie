function W=sparceweightmatrix(image,sigma,cropsize,nbsparce)
    imagesize=size(image,1);
    Wtild=weightmatrix(image,sigma,cropsize);
    W=zeros(imagesize^2);
    for i=1:imagesize^2
        aux=Wtild(i,:);
        for s=1:nbsparce
            [w indice]=max(aux);
            W(i,indice)=w;
            aux=aux([1:(indice-1) (indice+1):size(aux,2)]);
        end
    end
end
