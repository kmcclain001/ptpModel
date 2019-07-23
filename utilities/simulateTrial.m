% Simulates basic trial where animal runs through place field
%
% iputs:
%   speed - running speed on normalized scale where 0 is stationary and 1
%           is fast
%   thetaFreq - frequency of theta oscillation
%   initPhase - initial theta phase at entrace to palce field
%   fieldWidth - width of field in cm
%   
% outputs:
%   inputVars - matrix of input variables for model in the form [pos,
%               phase, speed] where position is normalized from 0 to 1 in
%               the field and phase is 0-2pi 
%

function inputVars = simulateTrial(speed, thetaFreq, initPhase, fieldWidth)

time = linspace(0,2,2500);

dist = (100/fieldWidth)*speed*time;
time = time(dist<=1);
x = dist(dist<=1);

phase = mod(initPhase+time*2*pi*thetaFreq,2*pi);
speed_trial = ones(1,length(x))*speed;

inputVars = [x',phase',speed_trial'];

end