getd = @(p)path(p,path);
getd('toolbox_signal/');
getd('toolbox_general/');
name = 'lena';
n = 256;
f0 = load_image(name);
f0 = rescale(crop(f0,n));
figure(1),
imageplot(f0);

rho = .8;
Lambda = rand(n,n)>rho;

Phi = @(f)f.*Lambda;
y = Phi(f0);

figure(2),
imageplot(y);

K  = @(f)grad(f);
KS = @(u)-div(u);

Amplitude = @(u)sqrt(sum(u.^2,3));
F = @(u)sum(sum(Amplitude(u)));

ProxF = @(u,lambda)max(0,1-lambda./repmat(Amplitude(u), [1 1 2])).*u;


t = -linspace(-2,2, 201);
[Y,X] = meshgrid(t,t);
U = cat(3,Y,X);
V = ProxF(U,1);
figure(3),
surf(V(:,:,1));
colormap jet(256);
view(150,40);

ProxFS = @(y,sigma)y-sigma*ProxF(y/sigma,1/sigma);

V = ProxFS(U,1);
figure(4),
surf(V(:,:,1));
colormap jet(256);
view(150,40);
axis('tight');
camlight; shading interp;
axis('tight');
camlight; shading interp;

ProxG = @(f,tau)f + Phi(y - Phi(f));

L = 8;
sigma = 10;
tau = .9/(L*sigma);
theta = 1;

f = y;
g = K(y)*0;
f1 = f;

fold = f;
g = ProxFS( g+sigma*K(f1), sigma);
f = ProxG(  f-tau*KS(g), tau);
f1 = f + theta * (f-fold);

niter = 200;
E = []; C = [];
for i=1:niter    
    % update
    fold = f;
    g = ProxFS( g+sigma*K(f1), sigma);
    f = ProxG(  f-tau*KS(g), tau);
    f1 = f + theta * (f-fold);
    % monitor the decay of the energy
    E(i) = F(K(f));
    C(i) = snr(f0,f);
end
figure(5),
h = plot(E);
set(h, 'LineWidth', 2);
axis('tight');

figure(6),
imageplot(f);



n = 64;
name = 'square';
f0 = load_image(name,n);

a = 4;
Lambda = ones(n);
Lambda(end/2-a:end/2+a,:) = 0;
Phi = @(f)f.*Lambda;

figure(7),
imageplot(f0, 'Original', 1,2,1);
imageplot(Phi(f0), 'Damaged', 1,2,2);


y = Phi(f0);
ProxG = @(f,tau)f + Phi(y - Phi(f));
niter = 600;
ndisp = round(linspace(1,niter, 5)); ndisp(1) = [];
E = [];
f = y;
g = K(y)*0;
f1 = f;
q = 1;
figure(8),
for i=1:niter    
    % update
    fold = f;
    g = ProxFS( g+sigma*K(f1), sigma);
    f = ProxG(  f-tau*KS(g), tau);
    f1 = f + theta * (f-fold);
    % monitor the decay of the energy
    E(i) = F(K(f));
    if i==ndisp(q)
        subplot(2,2,q);
        imageplot(f);
        q = q+1;
    end
end