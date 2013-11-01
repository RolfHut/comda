clear all
close all
clc

N=1000;
n=3;

m1=zeros(n,1);
m2=zeros(n,n,1);

ensemble=zeros(n,N);

t=1;

for ensembleCounter=1:N
    Psi_f=mvnrnd([0 0 0],[0.99 0 0 ; 0 0.99 0 ; 0 0 0.99])';
    ensemble(:,ensembleCounter)=Psi_f;
    
                m1(:,t)=m1(:,t)+Psi_f;
            m2(:,:,t)=m2(:,:,t)+(Psi_f*Psi_f');
            
end

P_f=(1/(N-1))*(m2(:,:,t)-(m1(:,t)*m1(:,t)'/N));

P_f_hat=cov(ensemble');
