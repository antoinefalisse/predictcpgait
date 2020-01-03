function act = integrateActivationDynamics(t,EMG,tact,tdeact,b)

ppe = spline(t,EMG);

e0 = EMG(1);

AbsTol = 1e-6;
RelTol = 1e-3;
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
act0 = e0;

[~,act] = ode15s(@ActivationWintersContinuousOde,t,act0,options,ppe,tact,tdeact,b);

n = size(EMG,1);
if length(act) < n
    fprintf('Warning: Numerical integration -Activation Dynamics- did not complete successfully\n')
    act = zeros(n,1);
end

% Constrain the activation between 0 and a_max
act(act < 0) = 0;
act(act > 1) = 1;
act(isnan(act)) = 0;

return

function dadt = ActivationWintersContinuousOde(t,act,ppe,tact,tdeact,b)

    e = ppval(ppe,t);
    d1 = 1/(tact*(0.5+1.5*act));
    d2 = (0.5+1.5*act)/tdeact;
    f = 0.5*tanh(b*(e-act));
    dadt = (d1*(f+0.5) + d2*(-f+0.5))*(e-act);
    
return