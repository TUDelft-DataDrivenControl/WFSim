function v = vaf(y,y_est) 
%VAF         Compute the percentage Variance Accounted For (VAF) 
%            between two signals. 
% 
% Syntax: 
%            v = vaf(y,y_estimate) 
% 
% Description: 
%            The VAF is calculated as: 
% 
%                         variance(y-y_est) 
%               v = ( 1 -  ----------------  ) * 100% 
%                           variance(y) 
% 
%            The VAF of two signals that are the same is 
%            100%. If they differ, the  VAF will be lower. 
%            When y and y_est have multiple columns, the VAF 
%            is calculated for every column in y and y_est. 
%            The VAF is often used to verify the 
%            correctness of a model, by comparing the real 
%            output with the estimated output of the model. 
% 
% Inputs: 
%  y         Signal 1, often the real output. 
%  y_est     Signal 2, often the estimated output of a model. 
% 
% Output: 
%  v         VAF, computed for the two signals 
 
% Bert Haverkamp, April 1996 
% Revised by Niek Bergboer, Februari 2002 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1996-2007, Delft Center of Systems and Control 

if nargin < 2
    error('VAF requires two input arguments.');
end
if size(y,2) > size(y,1)
    y = y';
end
if size(y_est,2) > size(y_est,1)
    y_est = y_est';
end

N = size(y,1);
if size(y_est,1) ~= N
    error('Both signals should have an equal number of samples.');
end
if size(y,2) ~= size(y_est,2)
    error('Both signals should have an equal number of components.');
end

v = max(diag(100*(eye(size(y,2))-cov(y-y_est)./cov(y))),0);
 
 
 

