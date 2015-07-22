function [ensembleMean,covarianceMatrix]=...
    comDA(model,observations,transformation,settings,n_timesteps,...
    n_modelStepsPerTimestep,N)
%runs comDA, needs documentations. the steps refer to the steps in the
%article  / tex document.


%create empty output matrices
ensembleMean=zeros(model.stateVectorSize,n_timesteps);
covarianceMatrix=zeros(model.stateVectorSize,model.stateVectorSize,n_timesteps);

%create empty placeholders. Note that saving the placeholders is not
%needed, but done here to compare comDA with EnKF.
m1=zeros(model.stateVectorSize,n_timesteps);
m2=zeros(model.stateVectorSize,model.stateVectorSize,n_timesteps);


firstTimeLoop=true;

t=1;
observationCounter=1;

%loop through the observation timesteps
while t < n_timesteps
    
    
    %determine the timesteps to evaluate till the next observation.
    if observationCounter==1;
        timesteps=1:observations.timestamp(observationCounter);
    elseif observationCounter>length(observations.timestamp); %observationCounter==1;
        timesteps=(observations.timestamp(observationCounter-1)+1):...
            n_timesteps;
    else
        timesteps=(observations.timestamp(observationCounter-1)+1):...
            observations.timestamp(observationCounter);
    end %observationCounter==1;
    
    
    
    %loop through the ensemble members
    for ensembleCounter=1:N
        
        
        %Step 1: draw an ensemble member from the distributions
        if firstTimeLoop; %draw from initial distribution
            Psi_0=mvnrnd(settings.mu_psi_0,settings.cov_psi_0)';
        else %draw from previous time step distribution
            Psi_0=mvnrnd(Psi_f_mu',P_f)';
        end %if n==1;
        
        %loop through the timesteps till the observation
        for t=timesteps
            
            %Step 2: run the model
            tSelect=(t-1)*n_modelStepsPerTimestep+(1:n_modelStepsPerTimestep);
            
            Psi_f=feval(model.model,model.parameters,Psi_0,n_modelStepsPerTimestep,...
                observations.forcingEnsemble(:,ensembleCounter,tSelect));
            
            %Step 3: update placeholders
            m1(:,t)=m1(:,t)+Psi_f;
            m2(:,:,t)=m2(:,:,t)+(Psi_f*Psi_f');
            
            %set input for next timestep.
            Psi_0=Psi_f;
            
        end %for t=timesteps
        
        
        %step 4: repeat for N_ensembles times
    end %for ensemble=1:N_ensembles
    
    %switch firstTimeLoop flag
    firstTimeLoop=false;
    
    if observationCounter <= length(observations.timestamp)
        
        %step 5: calculate Psi_f_mu en P_f
        Psi_f_mu=m1(:,t)/N;
        P_f=(1/(N-1))*(m2(:,:,t)-(m1(:,t)*m1(:,t)'/N));
        
        %step 6: calculate kalman Gain
        %the observation ensemble:
        D=observations.ensemble(:,:,observationCounter);
        eps_D=D*(eye(N)-(1/N)*ones(N));
        R=(1/(N-1))*(eps_D*eps_D');
        K=(P_f*transformation.H')/(transformation.H*P_f*transformation.H'+R);
        
        
        %step 7: calculate Psi_a_mu and P_a, here called Psi_f_mu and P_f,
        %because it saves memory for the next loop-iteration not to have to
        %store it double
        Psi_f_mu=Psi_f_mu+K * (observations.obs(:,observationCounter)-...
            transformation.H * Psi_f_mu);
        P_f=(eye(model.stateVectorSize)-K*transformation.H)*P_f*...
            (eye(model.stateVectorSize)-K*transformation.H)'+(K*R*K');
        
        %for numerical stability, make P_f symetrical
        P_f=0.5*(P_f+P_f');
        
        
        ensembleMean(:,t)=Psi_f_mu;
        covarianceMatrix(:,:,t)=P_f;
        
        
        observationCounter=observationCounter+1;
        
    end %if observationCounter < length(observations.timestamp)
    
end %while t <= n_timesteps



for t=1:n_timesteps
    if sum(observations.timestamp==t)==0
        ensembleMean(:,t)=m1(:,t)/N;
        covarianceMatrix(:,:,t)=(1/(N-1))*(m2(:,:,t)-(m1(:,t)*m1(:,t)'/N));
    end %if sum(observations.timestamp==t)==0
end %for t=1:n_timesteps
