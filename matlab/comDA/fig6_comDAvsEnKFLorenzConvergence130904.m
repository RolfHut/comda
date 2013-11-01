%% doc
% This script shows the similarities and differences between KF, EnKF and
% comDA.
% This is done by running two models: a linear AR1 model and a Lorenz
% model. Both are run within the KF, EnKF and comDA schemes. This is done
% for different numbers of ensembles.
% the models are functions located in the libDir. Models are functions that
% take a statevector and a time-step as inputs and returns an updated
% statevector

%% prelim
clc
clear all

%% settings
projectDir='/Users/rwhut/Documents/TU/eWaterCycle/matlab/comDA';
libDir='/Users/rwhut/Documents/TU/eWaterCycle/matlab/lib';
figdir=[projectDir filesep 'fig'];

addpath(libDir);
%% parameters



%total number of timesteps to run
n_timesteps=500;
n_modelStepsPerTimestep=[25];

%time axis (for plotting)
dt=1;
tAxis=dt*(1:n_timesteps);

%observation timestamps
observations.timestamp=50:50:n_timesteps;

%the actual model
model.model=@Lorenz;
model.stateVectorSize=3;
model.parameters.sigma = 10;
model.parameters.rho = 28;
model.parameters.beta = 8.0/3;
model.parameters.dt=1e-3;

%the structure relating model space to measurement space
% in this simple example: states 1 and 2 out of 3 are observed, so the model space
% has 3 dimensions, the observations space has 2. The upper left square of
% the H matrix is diagonal, because states 1 and 2 are directly observed.
transformation.observedStates=[1;2];%ones(model.stateVectorSize,1);
transformation.H=zeros(2,3);
transformation.H(1:2,1:2)=eye(length(transformation.observedStates));



%number of ensemble members
N=[5 10 20 50 100 200 500];

%number of runs per combinations of
runNr=50;

%which of the states to print in the figures
plotParameter=1;

%% settings/assumptions needed by the different schemes
%standard deviation (error) in observations
settings.sigma_d=[1;1]*[1];

%forcing error, standard deviation of observations of the forcings
observations.forcingError=[1;1;1];



%% derived size quantities, following Everson

%N=N
m=length(transformation.observedStates);
n=model.stateVectorSize;

%and derived by me
m_timesteps=length(observations.timestamp);


%% loop through the variables that change for the different sub-plots:
% n_modelStepsPerTimestep and settings.sigma_d

subPlotCounter=0;

%variables to store results to be plotted
EnKFImprovement=zeros(length(N),runNr);
comDAImprovement=zeros(length(N),runNr);


for N_counter=1:length(N);
    
    subPlotCounter=subPlotCounter+1;
    
    
    for runCounter=1:runNr; %repeat a few times to get more data.
        
        %% settings that change every run
        
        %the starting state vector
        psi_0=10*randn(model.stateVectorSize,1);
        
        %mean in starting state vector
        settings.mu_psi_0=psi_0;
        
        %covariance in starting state vector. The line below assumes equal
        %variance for all the observations and applies that variance to the
        %initial ensemble.
        settings.cov_psi_0=settings.sigma_d(1,1)*eye(model.stateVectorSize);
        
        %% create truth
        
        
        %assume that the model describes the true proces
        truth.model=model.model;
        truth.parameters=model.parameters;
        
        %true forcing
        truth.forcing=20*randn(n,n_timesteps);
        
        %true states, using true model and true forcing.
        truth.state=zeros(n,n_timesteps);
        t=0;
        for t_step=1:n_timesteps
            t=t+1;
            if t==1;
                truth.state(:,t)=feval(truth.model,truth.parameters,psi_0,n_modelStepsPerTimestep);
            else
                truth.state(:,t)=feval(truth.model,truth.parameters,truth.state(:,t-1),...
                    n_modelStepsPerTimestep,truth.forcing(:,t));
            end %if n==1;
        end %for t_step=1:n_timesteps
        
        %% create observations from truth
        
        %the actual observations (ie, not an ensemble based on the observations)
        observations.obs=truth.state(transformation.observedStates,observations.timestamp)+...
            (settings.sigma_d(transformation.observedStates)*ones(1,m_timesteps).*randn(m,m_timesteps));
        
        %the covariance of the measurement errors (ie. gamma matric)
        % this is either a dim2 matrix if the covariance is constant for all
        % (observation) timesteps, or is a dim3 matrix if it varies per
        % timestep.
        observations.obsErrorCov=eye(m);
        
        %observed forcing
        observations.forcing=truth.forcing;
        
        
        
        
        %% run EnKF
        
        %create initial ensemble
        initial_ensemble=zeros(n,N(N_counter));
        for ensembleCounter=1:N(N_counter)
            initial_ensemble(:,ensembleCounter)=mvnrnd(settings.mu_psi_0,settings.cov_psi_0);
        end %for ensembleCounter=1:N
        
        %create observation ensemble
        observations.ensemble=zeros(m,N(N_counter),m_timesteps);
        for t_step=1:m_timesteps;
            observations.ensemble(:,:,t_step)=observations.obs(:,t_step)*ones(1,N(N_counter))+...
                (settings.sigma_d(transformation.observedStates)*ones(1,N(N_counter))).*randn(m,N(N_counter));
        end %for t_step=1:length(observations.timestamp);
        
        %create forcing ensemble
        observations.forcingEnsemble=zeros(n,N(N_counter),n_timesteps);
        for t_step=1:n_timesteps;
            observations.forcingEnsemble(:,:,t_step)=observations.forcing(:,t_step)*ones(1,N(N_counter))+...
                (observations.forcingError*ones(1,N(N_counter))).*randn(n,N(N_counter));
        end %for t_step=1:length(observations.timestamp);
        
        %run the EnKF
        
        if (N_counter>1)
            ensemble=EnKF(model,observations,transformation,initial_ensemble,...
                n_timesteps,n_modelStepsPerTimestep,N(N_counter));
            
            %calculate statistics
            EnKFEnsembleMean=permute(mean(ensemble,2),[1 3 2]);
            EnKFEnsembleStd=permute(std(ensemble,[],2),[1 3 2]);
            
            %decrease in std after update
            %improveStatistic=EnKFEnsembleStd(plotParameter,observations.timestamp(2:end)-1)-EnKFEnsembleStd(plotParameter,observations.timestamp(2:end));
            
            %RMS of entire run
            improveStatistic=sqrt(mean((EnKFEnsembleMean(plotParameter,:)-...
                truth.state(plotParameter,:)).^2));
            
            EnKFImprovement(N_counter,runCounter)=improveStatistic;
        else
            EnKFImprovement(N_counter,runCounter)=NaN;
        end
        
        
        %% run comDA
        
        [comDAEnsembleMean,comDACovarianceMatrix]=...
            comDA(model,observations,transformation,settings,n_timesteps,...
            n_modelStepsPerTimestep,N(N_counter));
        
        comDAStd=zeros(n,n_timesteps);
        for t=1:n_timesteps
            comDAStd(:,t)=sqrt(diag(comDACovarianceMatrix(:,:,t)));
        end %for t=1:n_timesteps
        
        %decrease in std
        %improveStatistic=comDAStd(plotParameter,observations.timestamp(2:end)-1)-comDAStd(plotParameter,observations.timestamp(2:end));
        
        %RMS of entire run
        improveStatistic=sqrt(mean((comDAEnsembleMean(plotParameter,:)...
            -truth.state(plotParameter,:)).^2));
        comDAImprovement(N_counter,runCounter)=improveStatistic;
        
        
    end %for runCounter=1:runNr;
    
end %for N_counter=1:length(N);

%% add info to figure

tic;while (toc<5);beep;end

%% make figure


figure(2);
subplot(1,2,1);
boxplot(log10(EnKFImprovement'),'position',log10(N));
set(gca,'XTick',log10(N))
set(gca,'XTickLabel',N)
set(gca,'YTick',log10([0.5 1 2 5 10 20]))
set(gca,'YTickLabel',[0.5 1 2 5 10 20])
set(gca,'YLim',log10([0.2 20]))
xlabel('number of ensemble members')
ylabel('RMS in EnKF')

subplot(1,2,2);
boxplot(log10(comDAImprovement'),'position',log10(N));
set(gca,'XTick',log10(N))
set(gca,'XTickLabel',N)
set(gca,'YTick',log10([0.5 1 2 5 10 20]))
set(gca,'YTickLabel',[0.5 1 2 5 10 20])
set(gca,'YLim',log10([0.2 20]))
xlabel('number of ensemble members')
ylabel('RMS in comDA')

print(gcf,[figdir filesep 'fig6_comDAvsEnKFLorenzConvergence.eps'],'-depsc');
