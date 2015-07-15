function state_out = lorenz4D(parameters,state_in,n_timesteps,forcing)
% X (txJ) = lorenz4D(tf,forcing,perturbation)
% ensure final time is divisible by h

%%% the Lorenz model is: (cyclical)
% dX[j]/dt=(X[j+1]-X[j-2])*X[j-1]-X[j]+F
J=parameters.J; %default 40;               %the number of variables
h=parameters.h; %default 0.05;             %the time step
F=parameters.F; %default 8



for i=1:n_timesteps %for each time
  state_out=state_in+rk4(state_in,h,F); %solved via RK4
  state_in=state_out;
end


function deltay = rk4(Xold,h,F)
% X[t+1] = rk4(X[t],step)
 k1 = f(Xold,F);
 k2 = f(Xold+1/2.*h.*k1,F);
 k3 = f(Xold+1/2.*h.*k2,F);
 k4 = f(Xold+h.*k3,F);
 deltay= 1/6*h*(k1+2*k2+2*k3+k4);

function k = f(X,F)

k=zeros(size(X));
J=length(X);
%first the 3 problematic cases: 1,2,J
k(1)=(X(2)-X(J-1))*X(J)-X(1);
k(2)=(X(3)-X(J))*X(1)-X(2);
k(J)=(X(1)-X(J-2))*X(J-1)-X(J);
%then the general case
for j=3:J-1
 k(j)=(X(j+1)-X(j-2)).*X(j-1)-X(j);
end
%add the F    
k=k+F;