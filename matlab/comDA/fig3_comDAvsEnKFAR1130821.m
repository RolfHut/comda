%% doc
% This script shows the differences (almost none) between EnKF and comDA
% when applied to a 3 dimensional AR1 model with 2 observed states

%% prelim
clc
clear all
close all

%% settings
projectDir='/Users/rwhut/Documents/TU/eWaterCycle/matlab/comDA';
libDir='/Users/rwhut/Documents/TU/eWaterCycle/matlab/lib';
figdir=[projectDir filesep 'fig'];

addpath(libDir);
%% parameters

%total number of timesteps to run
n_timesteps=500;
n_modelStepsPerTimestep=[1 4];

%time axis (for plotting)
dt=1;
tAxis=dt*(1:n_timesteps);

%observation timestamps
observations.timestamp=50:50:500;

%the actual model
model.model=@AR1Filter;
model.stateVectorSize=3;
model.parameters.AR1coefs=[0.5 0.2 0.3 ; 0.2 0.5 0.3 ; 0.3 0.3 0.4];

%the structure relating model space to measurement space
% in this simple example: diagonal matrix: all states are observed
transformation.observedStates=[1;2];%ones(model.stateVectorSize,1);
transformation.H=zeros(2,3);
transformation.H(1:2,1:2)=eye(length(transformation.observedStates));

%the starting state vector
psi_0=[0;0;0];


%number of ensemble members
N=100;

%number of runs per combinations of
runNr=20;

%which state to plot
plotParameter=1;

%% settings/assumptions needed by the different schemes
%mean in starting state vector
settings.mu_psi_0=psi_0;
%covariance in starting state vector
settings.cov_psi_0=model.parameters.AR1coefs;
%standard deviation (error) in observations
settings.sigma_d=[1;1]*[1 10];

%forcing error, standard deviation of observations of the forcings
observations.forcingError=[1;1;1];



%% derived size quantities, following Everson

%N=N
m=length(transformation.observedStates);
n=model.stateVectorSize;

%and derived by me
m_timesteps=length(observations.timestamp);


subPlotCounter=0;

for n_modelStepsPerTimestepCounter=1:length(n_modelStepsPerTimestep);
    for sigma_dCoutner=1:size(settings.sigma_d,1);
        
        EnKFImprovement=[];
        comDAImprovement=[];
        subPlotCounter=subPlotCounter+1;
        
        
        for runCounter=1:runNr; %repeat a few times to get more data.
            
            %% create truth
            
            
            %assume that the model describes the true proces
            truth.model=model.model;
            truth.parameters=model.parameters;
            
            %true forcing
            truth.forcing=randn(n,n_timesteps*n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter));
            
            %true states, using true model and true forcing.
            truth.state=zeros(n,n_timesteps);
            t=0;
            for t_step=1:n_timesteps
                t=t+1;
                if t==1;
                    truth.state(:,t)=feval(truth.model,truth.parameters,psi_0...
                        ,n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter)...
                        ,truth.forcing(:,t*n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter)+...
                        (0:n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter)-1)));
                else
                    truth.state(:,t)=feval(truth.model,truth.parameters,truth.state(:,t-1),...
                        n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter),truth.forcing(:,...
                        t*n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter)+...
                        (-(n_modelStepsPerTimestep(n_modelStepsPerTimestepCounter)-1):0)));
                end %if n==1;
            end %for t_step=1:n_timesteps
            
            %% create observations from truth
            
            %the actual observations (ie, not an ensemble based on the observations)
            observations.obs=truth.state(transformation.observedStates,observations.timestamp)+...
                (settings.sigma_d(transformation.observedStates,sigma_dCoutner)*ones(1,m_timesteps).*randn(m,m_timesteps));
            
            %the covariance of the measurement errors (ie. gamma matric)
            % this is either a dim2 matrix if the covariance is constant for all
            % (observation) timesteps, or is a dim3 matrix if it varies per
            % timestep.
            observations.obsErrorCov=eye(m);
            
            %observed forcing
            observations.forcing=truth.forcing;
            
            
            
            
            %% run EnKF
            
            %create initial ensemble
            initial_ensemble=zeros(n,N);
            for ensembleCounter=1:N
                initial_ensemble(:,ensembleCounter)=mvnrnd(settings.mu_psi_0,settings.cov_psi_0);
            end %for ensembleCounter=1:N
            
            %create observation ensemble
            observations.ensemble=zeros(m,N,m_timesteps);
            for t_step=1:m_timesteps;
                observations.ensemble(:,:,t_step)=observations.obs(:,t_step)*ones(1,N)+...
                    (settings.sigma_d(transformation.observedStates,sigma_dCoutner)*ones(1,N)).*randn(m,N);
            end %for t_step=1:length(observations.timestamp);
            
            %create forcing ensemble
            observations.forcingEnsemble=zeros(n,N,n_timesteps);
            for t_step=1:n_timesteps;
                observations.forcingEnsemble(:,:,t_step)=observations.forcing(:,t_step)*ones(1,N)+...
                    (observations.forcingError*ones(1,N)).*randn(n,N);
            end %for t_step=1:length(observations.timestamp);
            
            %run the EnKF
            
            ensemble=EnKF(model,observations,transformation,initial_ensemble,n_timesteps,...
                n_modelStepsPerTimestep,N);
            
            %calculate statistics
            EnKFEnsembleMean=permute(mean(ensemble,2),[1 3 2]);
            EnKFEnsembleStd=permute(std(ensemble,[],2),[1 3 2]);
            %RMS of entire run
            improveStatistic=sqrt(mean((EnKFEnsembleMean(plotParameter,:)-...
                truth.state(plotParameter,:)).^2));
            
            EnKFImprovement=[EnKFImprovement improveStatistic];
            
            
            %% run comDA
            
            [comDAEnsembleMean,comDACovarianceMatrix]=...
                comDA(model,observations,transformation,settings,n_timesteps,...
                n_modelStepsPerTimestep,N);
            
            comDAStd=zeros(n,n_timesteps);
            for t=1:n_timesteps
                comDAStd(:,t)=sqrt(diag(comDACovarianceMatrix(:,:,t)));
            end %for t=1:n_timesteps
            
            %RMS of entire run
            improveStatistic=sqrt(mean((comDAEnsembleMean(plotParameter,:)...
                -truth.state(plotParameter,:)).^2));
            comDAImprovement=[comDAImprovement improveStatistic];
            
            
        end %for runCounter=1:runNr;
        
        %scatter EnKF results vs comDA results, make the figure.
        figure(1);
        subplot(length(n_modelStepsPerTimestep),size(settings.sigma_d,1),subPlotCounter);
        scatter(EnKFImprovement,comDAImprovement);
        axis([0 10 0 10]);
        hold on
        plot(0:1:10,0:1:10,'r')
        drawnow
        
        
    end %for sigma_dCoutner=1:size(settings.sigma_d,1);
end %for n_modelStepsPerTimestepCounter=1:length(n_modelStepsPerTimestep);



%% add info to figure

subplot(2,2,1);
title('(a)');
ylabel('RMS in comDA')
subplot(2,2,2);
title('(b)');
subplot(2,2,3);
ylabel('RMS in comDA')
xlabel('RMS in EnKF');
title('(c)');
subplot(2,2,4);
xlabel('RMS in EnKF');
title('(d)');
print(gcf,[figdir filesep 'fig3_comDAvsEnKFAR1.eps'],'-depsc');

