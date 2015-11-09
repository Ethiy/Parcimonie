clear; close all;

name = 'lena';

n0 = 512;

f = rescale(load_image(name,n0));

figure(1), imageplot( f, 'Image f');
%figure(2),
%imageplot( crop(f,64), 'Zoom' );


%figure(3), imageplot(-f, '-f', 1,2,1);
%imageplot(f(n0:-1:1,:), 'Flipped', 1,2,2);


k = 9; % size of the kernel
h = ones(k,k);
h = h/sum(h(:)); % normalize

fh = perform_convolution(f,h);

%figure(4),
%imageplot(fh, 'Blurred image');

% --------------------------------FFT----------------------------------

F = fft2(f) / n0;
% Conservation of energy:
disp(strcat(['Energy of Image:   ' num2str(norm(f(:)))]));
disp(strcat(['Energy of Fourier: ' num2str(norm(F(:)))]));


L = fftshift(log( abs(F)+1e-1 ));

%figure(5);
%imageplot(L, 'Log(Fourier transform)');

%============Linear Approximation:=====================================

Fm = fftshift(F);

M = n0^2/64;
q = sqrt(M);

G= zeros(n0);
G(n0/2-q/2:n0/2+q/2,n0/2-q/2:n0/2+q/2) = Fm(n0/2-q/2:n0/2+q/2,n0/2-q/2:n0/2+q/2); 
G = fftshift(G);

fM = real(ifft2(G));
figure(6); 
imageplot(fM, 'Linear approximation');

figure(7),
subplot(2,1,1);
plot(f(:,n0/2));
axis('tight'); title('f');
subplot(2,1,2);
plot(fM(:,n0/2));
axis('tight'); title('f_M');

%===========================Thresholding================================

T = .2;
F = fft2(f) / n0;
FT = F .* (abs(F)>T);
LT = fftshift(log( abs(FT)+1e-1 ));
figure(8),
imageplot(LT, 'thresholded Log(Fourier transform)');

fT = real( ifft2(FT)*n0 );
figure(9);
imageplot(clamp(fM), ['Non-linear, Fourier, SNR=' num2str(snr(f,fM), 4) 'dB']);

m = sum(FT(:)~=0);
disp(['M/N = 1/'  num2str(round(n0^2/m)) '.']);

%======================================
[a, I] = sort(abs(F(:)),'descend');
GM = zeros(n0);
GM(I(1:M)) = F(I(1:M));
gM = real(ifft2(GM));

TM=a(I(M)+1);

figure(10); 
imageplot(gM, 'NonLinear approximation');

%==============Wavelet Transform==================================

Jmin = 0;
fw = perform_wavelet_transf(f,Jmin,+1);

figure(11),
plot_wavelet(fw);


%---------------Linear-------------------------------------

j0 = 5;
M = n0^2/(2^j0);
m = sqrt(M);
FL = zeros(n0);

FL(1:m,1:m) = fw(1:m,1:m);

fl = perform_wavelet_transf(FL,Jmin,-1);
figure(12),
imageplot(clamp(fl), ['Linear approximation for j0=',int2str(j0),', M/N=' num2str(M/n0^2,2)...
    ', SNR=' num2str(snr(f,fl),3) 'dB' ]);

%=============Nonlinear approximation by thresholding======================

Th = .2;

fwT = fw.*(abs(fw)>Th);
figure(13)
subplot(1,2,1);
plot_wavelet(fw);
title('Original coefficients');

subplot(1,2,2); title('Thresholded coeffs')
plot_wavelet(fwT);

%Reconstruction:
fTh = perform_wavelet_transf(fwT,Jmin,-1);
figure(14),
imageplot(clamp(fTh), strcat(['Nonlinear Approximation with Th=' num2str(Th) ', SNR=' num2str(snr(f,fTh),3) 'dB']));

%=============Nonlinear approximation with M preset=======================

%we use the same value for M
[a,Im] = sort(abs(fw(:)),'descend');

FNL = zeros(n0);
FNL(Im(1:M)) = fw(Im(1:M));

Th_guessed = abs(FNL(Im(M)));

fM = perform_wavelet_transf(FNL,Jmin,-1);
figure(15); 
imageplot(clamp(fM), ...
    ['NonLinear approximation, M/N=' num2str(M/n0^2,2)...
    ', SNR=',num2str(snr(f,fM),3),'dB, threshold:=',num2str(Th_guessed)]);

figure(16);
subplot(2,1,1);
plot(f(:,n0/2));
axis('tight'); title('f');
subplot(2,1,2);
plot(fM(:,n0/2));
axis('tight'); title('f_M');



