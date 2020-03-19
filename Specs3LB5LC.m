% This Function is called in the ODE function to show the changes in u, w,
% p, and theta over time.
%   Author: Benjiman Smith
%   Collaborators: E. Owen, I. Quezada
%   Date: 2/20/2020
%
function dydt = Specs2LB4NLC(t, Conditions, A) % function start
    dydt = A*Conditions; % Utilize A matrix to find diferential results
 
end % end