%2D example.

clear all
close all
clc

%total number of samples, coming in as a stream
n=5000;

%number of samples that can be hold in memory.
l=20;

%dimension of problem
m=100;

%create sample: random x, linear relation with y, including noise
mu = 5*ones(1,m); 
Sigma = diag(rand(1,m));
Sigma=(Sigma+Sigma')*0.5;
A = mvnrnd(mu, Sigma, n);


B=zeros(l,m);
Am=zeros(1,m);
for k=1:n
    %new A comes availables. For now from a pre-matrix, later: a model
    %realisation
    Anew=A(k,:);
    
    %update running mean and running ensemble
    Am=(((k-1)*Am)+Anew)/k;
    [B]=updateSketch(l,Anew-Am,B);
    scatter(A(:,1),A(:,2),'.');hold on;scatter(Am(1)+B(:,1),Am(2)+B(:,2),'r.');scatter(Am(1,1),Am(1,2),'g.');hold off;
    axis([0 10 0 20]);
    drawnow();
end
