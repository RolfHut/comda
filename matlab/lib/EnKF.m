function [ensemble,covarianceMatrix]=...
    EnKF(model,observations,transformation,initial_ensemble,n_timesteps,...
    n_modelStepsPerTimestep,N)
% run EnKF
% basiccaly implements the Evenson paper. 
% model, observations and transformations are all structures. The following
% fields are assumed:
%
% For the model: 
%
% model.stateVectorSize, an integer indicating the number of parameters in
% the state vector. Note that this is the "complete statevector" from
% system engineering, not the "to be updated vectors" from
% data-assimilation. in the documentation (and code) below, this is called
% n, following the Everson paper.
%
% model.model, a function handle to the model. The function that runs the
% model should have the following structure:
% state_out=model(parameters,state_in,n_timesteps,forcing)
%
% model.parameters, the parameters that are given to the model.model
% function (see above)
%
% for the observations:
% 
% observations.timestamp, a vector with indices that indicate at which
% timestep observations are available. Note that the maximum of this vector
% must be less than (or equal to) n_timesteps
% 
% observations.ensemble, a 3D matrix with the observations, as ensemble.
% when the number of observations (in time) is m_timesteps and the number
% of observations per timestep is m, the size should be m x N x m_timesteps

% observations.forcingEnsemble, a 3D matrix with for every timestep a
% matrix containing the forcing ensemble. if k forcings are needed at every
% timestep, than the size should be k x N x n_timesteps
%
% for the transformations:
%
% transformation.H is the m x n matrix that transforms the statevector to
% the measurement space.
%
% transformation.observedStates is a masker vector, in size equal to a
% state vector, that contains 1 when a state needs to be considered in
% updating and a zero if it doesn't. This allows for efficient updating,
% where only those states that are deemed to be effected are taken into
% account. Note that in a lot of implementations, the vector containing the
% states t be updates is called "the state vector", and the states not
% updated are saved in "the start-up files".
%
% other inputs:
% 
% initial_ensemble is a n x N matrix containing the starting ensemble
% 
% n_timesteps is the number of timesteps that the model needs te be run.
% 
% N is the number of ensemble members. This must off course be consistend
% with the observations.ensemble, observations.forcingEnsemble and
% initial_ensemble inputs.
%



%create empty output matrices
ensemble=zeros(model.stateVectorSize,N,n_timesteps);
%only calculate covariance matrices if asked as output
if nargout==2
    covarianceMatrix=zeros(model.stateVectorSize,model.stateVectorSize,n_timesteps);
end %if nargout==2

observationCounter=1;


%loop through the timesteps
for t=1:n_timesteps
    
    %loop through the ensemble members
    for ensembleCounter=1:N
        %run the model
        tSelect=(t-1)*n_modelStepsPerTimestep+(1:n_modelStepsPerTimestep);
        if t==1;
            ensemble(:,ensembleCounter,t)=feval(model.model,model.parameters...
                ,initial_ensemble(:,ensembleCounter),n_modelStepsPerTimestep,observations.forcingEnsemble(:,ensembleCounter,tSelect));
        else
            ensemble(:,ensembleCounter,t)=feval(model.model,model.parameters...
                ,ensemble(:,ensembleCounter,t-1),n_modelStepsPerTimestep,observations.forcingEnsemble(:,ensembleCounter,tSelect));
        end %if n==1;
        
    end %for ensemble=1:N_ensembles
    
    %if an observations is available: update
    
    if t==observations.timestamp(observationCounter)
        
        %select matrices from ensembles
        A=ensemble(:,:,t);
        Apr=A*(eye(N)-(1/N)*ones(N));
        D=observations.ensemble(:,:,observationCounter);
        Dpr=D-transformation.H*A;
        gamma=D*(eye(N)-(1/N)*ones(N));
        
        
        %now calculate [U,S,V]=SVD(H*Apr+gamma)
        [Utemp,Sigma] = svd((transformation.H*Apr)+gamma);
        
        U=zeros(size(transformation.H,1),N);
        U(1:size(transformation.H,1),1:size(transformation.H,1))=Utemp;
        Sigma=Sigma*transpose(Sigma);
        Sigma=diag(diag(Sigma).^-1);
        LambdaInv = zeros(N,N);
        LambdaInv(1:size(Sigma,1),1:size(Sigma,1))=Sigma;
        
        X1=LambdaInv*transpose(U);
        X2=X1*Dpr;
        X3=U*X2;
        X4=transpose(transformation.H*Apr)*X3;
        ensemble(:,:,t)=A+Apr*X4;
        
        if observationCounter<length(observations.timestamp);
            observationCounter=observationCounter+1;
        end
    end %if t==observations.timestamp(observationCounter)
    
    %only calculate covariance matrices if asked as output
    if nargout==2
        covarianceMatrix(:,:,t)=cov(ensemble(:,:,t)');
    end
    
end %for t=1:n_timesteps