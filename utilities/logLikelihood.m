% log-likelihood function, computes probability that rate vector
% produced spikes

function LL = logLikelihood(rateVec,spikeVec,dt)

prob = ((exp(-dt*rateVec).*(dt*rateVec).^spikeVec))./factorial(spikeVec);
LL = sum(log(prob));

end