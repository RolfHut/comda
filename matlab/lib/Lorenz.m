function state_out=AR1Filter(parameters,state_in,n_timesteps,forcing)

state_out=state_in;

for timestep=1:n_timesteps
    state_out(1) = state_in(1)+parameters.dt*(parameters.sigma * (state_in(2) - state_in(1)));
    state_out(2) = state_in(2)+parameters.dt*(state_in(1) * (parameters.rho -state_in(3)) - state_in(2));
    state_out(3) = state_in(3)+parameters.dt*(state_in(1) * state_in(2) - parameters.beta* state_in(3));
    state_in=state_out;
end %for timestep=1:n_timesteps

