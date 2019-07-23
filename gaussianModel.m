function model = gaussianModel()

%%% equations for rate %%%

% spatial input
f = @(b,x_f) exp(b(1)+((-1/2)*((x_f-b(2))/b(3)).^2));

rate = @(b,inputVars) f(b,inputVars(:,1));


%%% fitting constraints %%%
% range from which to select initial point in multi start
initBounds = [1, 4;     %Ax
             .1, 1;     %x0
             .1, 1];    %sx
     
% range for fitting in fmincon
fitBounds = [0, Inf;    %Ax
             0, 1;      %x0
             0, 1];     %sx
             
         
%%% model creation %%%
model = pfModel(rate,3,'initBounds',initBounds,'fitBounds',fitBounds);

end