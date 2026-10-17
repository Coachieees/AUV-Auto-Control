function out = run_cmuauv_case(caseId, stopTime)
%RUN_CMUAUV_CASE  Initialise, select a trajectory case and simulate.
%   out = run_cmuauv_case(1)        station keeping
%   out = run_cmuauv_case(2)        straight line
%   out = run_cmuauv_case(3, 60)    circle, 60 s
%   out = run_cmuauv_case(4)        3D infinity-shape
if nargin < 1, caseId   = 3;  end
if nargin < 2, stopTime = 60; end
evalin('base', 'init_cmu_auv');          % loads all parameters (also clears base ws)
assignin('base', 'traj_case', caseId);   % override the case chosen in the script
out = sim('cmuauv_controller', 'StopTime', num2str(stopTime));
end
