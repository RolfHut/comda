%% doc
% This script shows the differences (almost none) between EnKF and comDA
% when applied to a 3 parameter Lorenz model with 2 observed states

%% prelim
clc
clear all
close all

%% settings
projectDir='/Users/rwhut/Documents/TU/eWaterCycle/github/eWaterCycle-comda/matlab/comDA';
libDir='/Users/rwhut/Documents/TU/eWaterCycle/github/eWaterCycle-comda/matlab/lib';
figdir=[projectDir filesep 'fig'];

addpath(libDir);
%% parameters

%total number of timesteps to run
n_timesteps=500;
n_modelStepsPerTimestep=100;


%observation timestamps
observations.timestamp=50:50:n_timesteps;

%the actual model
model.model=@Lorenz;
model.stateVectorSize=3;
model.parameters.sigma = 10;
model.parameters.rho = 28;
model.parameters.beta = 8.0/3;
model.parameters.dt=1e-3;

%time axis (for plotting)
dt=1;
tAxis=model.parameters.dt*n_modelStepsPerTimestep*(1:n_timesteps);


%the structure relating model space to measurement space
% in this simple example: diagonal matrix: all states are observed
transformation.observedStates=[1;2];%ones(model.stateVectorSize,1);
transformation.H=zeros(2,3);
transformation.H(1:2,1:2)=eye(length(transformation.observedStates));

%the starting state vector
psi_0=10*randn(3,1);


%number of ensemble members
N=100;

%which state to plot
plotParameter=1;


%% settings/assumptions needed by the different schemes
%mean in starting state vector
settings.mu_psi_0=psi_0;
%covariance in starting state vector
settings.cov_psi_0=eye(model.stateVectorSize);
%standard deviation (error) in observations
settings.sigma_d=[1;1];

%forcing error, standard deviation of observations of the forcings
observations.forcingError=[1;1;1];



%% derived size quantities, following Everson

%N=N
m=length(transformation.observedStates);
n=model.stateVectorSize;

%and derived by me
m_timesteps=length(observations.timestamp);

%% create truth


%assume that the model describes the true proces
truth.model=model.model;
truth.parameters=model.parameters;

%true forcing
truth.forcing=20*randn(n,n_timesteps*n_modelStepsPerTimestep);

%true states, using true model and true forcing.
truth.state=zeros(n,n_timesteps);

for t=1:n_timesteps
    tSelect=(t-1)*n_modelStepsPerTimestep+(1:n_modelStepsPerTimestep);
    if t==1;
        truth.state(:,t)=feval(truth.model,truth.parameters,psi_0,n_modelStepsPerTimestep,truth.forcing(:,tSelect));
    else
        truth.state(:,t)=feval(truth.model,truth.parameters,truth.state(:,t-1),n_modelStepsPerTimestep,truth.forcing(:,tSelect));
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
initial_ensemble=zeros(n,N);
for ensembleCounter=1:N
    initial_ensemble(:,ensembleCounter)=mvnrnd(settings.mu_psi_0,settings.cov_psi_0);
end %for ensembleCounter=1:N

%create observation ensemble
observations.ensemble=zeros(m,N,m_timesteps);
for t_step=1:m_timesteps;
    observations.ensemble(:,:,t_step)=observations.obs(:,t_step)*ones(1,N)+...
        (settings.sigma_d(transformation.observedStates)*ones(1,N)).*randn(m,N);
end %for t_step=1:length(observations.timestamp);

%create forcing ensemble
observations.forcingEnsemble=zeros(n,N,n_timesteps*n_modelStepsPerTimestep);
for t_step=1:(n_timesteps*n_modelStepsPerTimestep);
    observations.forcingEnsemble(:,:,t_step)=observations.forcing(:,t_step)*ones(1,N)+...
        (observations.forcingError*ones(1,N)).*randn(n,N);
end %for t_step=1:length(observations.timestamp);

%run the EnKF

ensemble=EnKF(model,observations,transformation,initial_ensemble,...
    n_timesteps,n_modelStepsPerTimestep,N);

%calculate statistics
EnKFEnsembleMean=permute(mean(ensemble,2),[1 3 2]);
EnKFEnsembleStd=permute(std(ensemble,[],2),[1 3 2]);

%% run comDA

[comDAEnsembleMean,comDACovarianceMatrix]=...
    comDA(model,observations,transformation,settings,n_timesteps,n_modelStepsPerTimestep,N);

comDAStd=zeros(n,n_timesteps);
for t=1:n_timesteps
    comDAStd(:,t)=sqrt(diag(comDACovarianceMatrix(:,:,t)));
end %for t=1:n_timesteps

%% plot results
close all
figure(1);
subplot(2,1,1);
%plot truth
plot(tAxis,truth.state(plotParameter,:),'.')
hold on
%plot EnKF results
plot(tAxis,EnKFEnsembleMean(plotParameter,:))
plot(tAxis,EnKFEnsembleMean(plotParameter,:)+2*EnKFEnsembleStd(plotParameter,:),'r')
%plot comDA results
plot(tAxis,comDAEnsembleMean(plotParameter,:),'--')
plot(tAxis,comDAEnsembleMean(plotParameter,:)+2*comDAStd(plotParameter,:),'--r')

plot(tAxis,EnKFEnsembleMean(plotParameter,:)-2*EnKFEnsembleStd(plotParameter,:),'r')
plot(tAxis,comDAEnsembleMean(plotParameter,:)-2*comDAStd(plotParameter,:),'--r')

legend('truth','EnKF Ensemble Mean',...
    'EnKF 95% ensemble interval','comDA Ensemble Mean','comDA 95% ensemble interval',...
    'Location','SouthWest');
subplot(2,1,2)
plot(tAxis,EnKFEnsembleStd(plotParameter,:),tAxis,comDAStd(plotParameter,:),'r');
legend('standard deviation EnKF','standard deviation comDA','Location','NorthWest');
xlabel('time');

print(gcf,[figdir filesep 'fig2_comDAvsEnKFLorenz.eps'],'-depsc');


