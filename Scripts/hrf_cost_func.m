function J = hrf_cost_func(t0,tau1,tau2,A,B,sr,hrf_l,X_mat,y)
% hrf_cost_func: generate a measure of residual error between a signal 
%   predicted from an estimated HRF and the actual measured signal
%   written by Harrison Fisher 
%
%   Arguments:
%       HRF parameters that control the shape of the hrf:
%           t0: delay time in s
%           tau1: shape parameter for the positive alpha function
%           tau2: shape parameter for the negative alpha function
%           A:  amplitude parameter for the positive alpha function
%           B:  amplitude parameter for the negative alpha function
%           hrf_l: length of the HRF in seconds
%
%       X_mat: design matrix from signal (calcium) used to predict y (Hb)
%       y: measured signal to be deconvolved (Hb)
%       sr: sampling rate in Hz
    
    [hrf, ~, ~] = modified_alpha_hrf_2(t0,tau1,tau2,A,B,sr,hrf_l);
%     J = norm(y - X_mat*hrf)^2;
    conv_result = conv(X_mat,hrf);
    conv_result = conv_result(1:length(X_mat));
    J = norm(y - conv_result)^2;
end


