getd = @(p)path(p,path);

getd('toolbox_signal/');
getd('toolbox_general/');

%% Parameters

name='lena';
imagesize=40;
cropsize=5;
sigma=50;
nbsparce=7;

%% Loading image

image=load_image(name);
image=rescale(crop(image,imagesize));

%% Data Matrix

W=sparceweightmatrix(image,sigma,cropsize,nbsparce);
A=auxilarymatrix(W);

%% Stage A
l=100;
q=0;
tic
Q=randomizedRangeFinder(A,l,q);
[S,U,V]=directsvd(A,Q);
toc
s1=diag(S);
q=4;
tic
Q=randomizedRangeFinder(A,l,q);
[S,U,V]=directsvd(A,Q);
toc
s2=diag(S); 
tic
Q=fastrandomizedrangefinder(A,l);
[S,U,V]=directsvd(A,Q);
toc
s3=diag(S); 
plot(1:l,s1(1:l),'b');