getd = @(p)path(p,path);
getd('toolbox_signal/');
getd('toolbox_general/');
n = 128;
name = 'lena';
f0 = load_image(name);
f0 = rescale(crop(f0,n));
clf;
imageplot(f0, 'Image f_0');
rho = .7;
Omega = zeros(n,n);
sel = randperm(n^2);
Omega(sel(1:round(rho*n^2))) = 1;
Phi = @(f,Omega)f.*(1-Omega);
y = Phi(f0,Omega);
clf;
imageplot(y, 'Observations y');
SoftThresh = @(x,T)x.*max( 0, 1-T./max(abs(x),1e-10) );
clf;
T = linspace(-1,1,1000);
plot( T, SoftThresh(T,.5) );
axis('equal');
Jmax = log2(n)-1;
Jmin = Jmax-3;
options.ti = 0; % use orthogonality.
Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
SoftThreshPsi = @(f,T)Psi(SoftThresh(PsiS(f),T));
clf;
imageplot( clamp(SoftThreshPsi(f0,.1)) );

lambda = .03;

ProjC = @(f,Omega)Omega.*f + (1-Omega).*y;
fSpars = y;
fSpars = ProjC(fSpars,Omega);
fSpars = SoftThreshPsi( fSpars, lambda );

fSpars = y;
energy = [];
niter = 1000;
for i=1:niter
    fSpars = SoftThreshPsi( ProjC(fSpars,Omega), lambda );
    % record the energy
    fW = PsiS(fSpars);
    energy(i) = 1/2 * norm(y-Phi(fSpars,Omega), 'fro')^2 + lambda * sum(abs(fW(:)));
end
clf;
h = plot(energy);
axis('tight');
set_label('Iteration', 'E');
if using_matlab()
    set(h, 'LineWidth', 2);
end

clf;
imageplot(clamp(fSpars));
niter = 1000;
lambda_list = linspace(.03,0,niter);
err = [];
for i=1:niter
    fSpars = SoftThreshPsi( ProjC(fSpars,Omega), lambda_list(i) );    
end
clf;
imageplot(clamp(fSpars), ['Sparsity inpainting, SNR=' num2str(snr(f0,fSpars),3) 'dB']);

J = Jmax-Jmin+1;
u = [4^(-J) 4.^(-floor(J+2/3:-1/3:1)) ];
U = repmat( reshape(u,[1 1 length(u)]), [n n 1] );
lambda = .01;
options.ti = 1; % use translation invariance
Xi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
Psi = @(a)Xi(a./U);
tau = 1.9*min(u);
a = U.*PsiS(fSpars);
fTI = Psi(a);
a = a + tau*PsiS( Phi( y-Phi(fTI,Omega),Omega ) );
a = SoftThresh( a, lambda*tau );
niter = 1000;
a = U.*PsiS(fSpars);
E = [];
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI,Omega);
    E(i) = 1/2*norm(d , 'fro')^2 + lambda * sum( abs(a(:)) );   
    % step 
    a = SoftThresh( a + tau*PsiS(Phi(d,Omega)), lambda*tau );
end
clf;
plot(E); axis('tight');
fTI = Psi(a);

clf;
imageplot(clamp(fTI));
niter = 3000;
lambda_list = linspace(.03,0,niter);
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI,Omega);
    % step 
    a = SoftThresh( a + tau * PsiS( Phi(d,Omega) ) , lambda_list(i)*tau ) ;
end
clf;
imageplot(clamp(fTI), ['Sparsity inpainting TI, SNR=' num2str(snr(f0,fTI),3) 'dB']);
HardThresh = @(x,t)x.*(abs(x)>t);
t = linspace(-1,1,1000);
plot( t, HardThresh(t,.5) );
axis('equal');
niter = 500;
lambda_list = linspace(1,0,niter);
fHard = y;
fHard = ProjC(fHard,Omega);
fHard = Xi( HardThresh( PsiS(fHard), tau*lambda_list(1) ) );
niter = 500;
lambda_list = linspace(1,0,niter);
fHard = y; 
for i=1:niter
    fHard = Xi( HardThresh(PsiS(ProjC(fHard,Omega)), lambda_list(i)) );
end
clf;
imageplot(clamp(fHard), ['Inpainting hard thresh., SNR=' num2str(snr(f0,fHard),3) 'dB']);