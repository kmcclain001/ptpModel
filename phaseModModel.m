function model = phaseModModel()

%%% equations for rate %%%

% spatial input
f = @(b,x_f) exp(b(1)+((-1/2)*((x_f-b(2))/b(3)).^2));

% phase modulation
g = @(b,phs_g,x_g) exp(b(4)*(cos(phs_g-b(5))-1));

rate = @(b,inputVars) f(b,inputVars(:,1)).*g(b,inputVars(:,2),inputVars(:,1));


%%% fitting constraints %%%
% range from which to select initial point in multi start
initBounds = [1, 4;     %Ax
             .1, 1;     %x0
             .1, 1;     %sx
             -1, 3;     %kth
             .1, 2*pi]; %bth
     
% range for fitting in fmincon
fitBounds = [0, Inf;    %Ax
             0, 1;      %x0
             0, 1;      %sx
             0, Inf;    %kth
             0 2*pi];   %bth
             
         
%%% model creation %%%
model = pfModel(rate,5,'initBounds',initBounds,'fitBounds',fitBounds);

end