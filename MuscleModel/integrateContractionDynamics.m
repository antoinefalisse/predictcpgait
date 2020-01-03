function [FT,error] = integrateContractionDynamics(t,a,lMT,vMT,Parameters,Fvparam,Faparam,Kpe,integrator) 

a(a < 0.001) = 0.001;
ppa = spline(t,a);
pplMT = spline(t,lMT);
ppvMT = spline(t,vMT);

fprintf('Solving for tendon forces using tendon force numerical integration method . . .\n')

% Initial guess
[~,FT0,~, ~, ~, ~, ~, ~] = HillEquilibirium_RigidTendon(a(1),lMT(1),vMT(1),Parameters,Fvparam,Faparam,Kpe);
        
% Set integrator absolute and relative integrator error tolerances
AbsTol = 1e-5;
RelTol = 1e-3;
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);               
% Perform numerical integration with specified integrator
switch integrator                    
    case 'explicit'
        [~,FT] = ode45(@HillEquilibrium_TendonFL,t,FT0,options,ppa,pplMT,ppvMT,Parameters,Fvparam,Faparam,Kpe);
                        
    case 'implicit'
        [~,FT] = ode15s(@HillEquilibrium_TendonFL,t,FT0,options,ppa,pplMT,ppvMT,Parameters,Fvparam,Faparam,Kpe);                       
end

% Set output FT values to zero if integration did not complete successfully
n = size(lMT,1);
error =0;
if length(FT) < n
    error = error+1;
    fprintf('Warning: Numerical integration - Contraction Dynamics - did not complete successfully\n')
    FT = zeros(n,1);
end
                
% Set tendon forces < 0 or = NaN equal to zero
FT(FT < 0) = 0;
FT(isnan(FT)) = 0;
                       
end