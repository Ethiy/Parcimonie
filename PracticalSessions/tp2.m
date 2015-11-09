getd = @(p)path(p,path);
getd('toolbox_signal/');
getd('toolbox_general/');

% On ne fait pas l'exo 5

%============Image blurring=========================================

n = 256;
name = 'lena';
name = 'mri';
name = 'boat';
f0 = load_image(name);
f0 = rescale(crop(f0,n));


s = 3;
x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
h = exp( (-X.^2-Y.^2)/(2*s^2) );
h = h/sum(h(:));

hF = real(fft2(h));

clf;
imageplot(fftshift(h), 'Filter', 1,2,1);
imageplot(fftshift(hF), 'Fourier transform', 1,2,2);

if using_matlab()
    Phi = @(x,h)real(ifft2(fft2(x).*fft2(h)));
end


lambda_list = linspace(.03, .2, 40);
err = [];
for i=1:length(lambda_list)
    lambda = lambda_list(i);
    fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
    err(i) = snr(f0,fSob);
end
clf;
hh = plot(lambda_list, err); axis tight;
set_label('lambda', 'SNR');
if using_matlab()
    set(hh, 'LineWidth', 2);
end
[tmp,i] = max(err);
lambda = lambda_list(i);
fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );

y0 = Phi(f0,h);

clf;
imageplot(f0, 'Image f0', 1,2,1);
imageplot(y0, 'Observation without noise', 1,2,2);


sigma = .02;
y = y0 + randn(n)*sigma;

clf;
imageplot(y0, 'Observation without noise', 1,2,1);
imageplot(clamp(y), 'Observation with noise', 1,2,2);

tau = 1.9 / ( 1 + lambda * 8 / epsilon);
fTV = y;
E = []; 
for i=1:niter
    % Compute the gradient of the smoothed TV functional.
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );
    % step
    e = Phi(fTV,h)-y;
    fTV = fTV - tau*( Phi(e,h) + lambda*G);
    % energy
    E(i) = 1/2*norm(e, 'fro')^2 + lambda*sum(d(:));
end
% display energy
clf;
plot(E); axis('tight');
set_label('Iteration 3', 'Energy');
%============Deconvolution with L2 Regularisation

yF = fft2(y);
lambda = 0.02;
fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );

clf;
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fL2), strcat(['L2 deconvolution, SNR=' num2str(snr(f0,fL2),3) 'dB']), 1,2,2);

snr1=[];
aux=0;

for lambda=0:.001:.03
    fl2p=real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
    aux=snr(f0,fl2p);
    snr1=[snr1,aux];
end

plot([0:.001:.03],snr1);
    

%=========deconvolution by Sobolev Regularisation=================

S = (X.^2 + Y.^2)*(2/n)^2;
lambda = 0.2;
fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );

imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fSob), strcat(['Sobolev deconvolution, SNR=' num2str(snr(f0,fSob),3) 'dB']), 1,2,2);

snr2=[];
aux=0;

for lambda=0:.001:.03
    fSobp = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
    aux=snr(f0,fSobp);
    snr2=[snr2,aux];
end
figure(2),
plot([0:.001:.03],snr1);


%============Deconvolution by Total Variation Regularization==============

epsilon = 0.4*1e-2;
lambda = 0.06;
tau = 1.9 / ( 1 + lambda * 8 / epsilon);
fTV = y;
niter = 600;


Gr = grad(fTV);
d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
G = -div( Gr./repmat(d, [1 1 2])  );

tv = sum(d(:));

e = Phi(fTV,h)-y;
fTV = fTV - tau*( Phi(e,h) + lambda*G);


    niter = 400;
lambda_list = linspace(1e-6,.01,20);
tau = 1.9 / ( 1 + max(lambda_list) * 8 / epsilon);
fBest = y; fTV = y;
err = [];
for it=1:length(lambda_list)
    lambda = lambda_list(it);
    for i=1:niter
        % Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV);
        d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
        G = -div( Gr./repmat(d, [1 1 2])  );
        % step
        e = Phi(fTV,h)-y;
        fTV = fTV - tau*( Phi(e,h) + lambda*G);
    end
    err(it) = snr(f0,fTV);
    if err(it)>snr(f0,fBest)
        fBest = fTV;
    end
end
clf;
plot(lambda_list,err);
axis('tight');
xlabel('\lambda'); ylabel('SNR');
fTV = fBest;

