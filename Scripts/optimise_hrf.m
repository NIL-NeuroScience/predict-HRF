function [conv_result, params] = optimise_hrf(A,B,t0,tau1,tau2,sr,hrf_l,jRGECO,HbT,opt_method)

conv_result = []; 
params = [];
fun = @(params)hrf_cost_func(params(1),params(2),params(3),params(4),params(5),sr,hrf_l,jRGECO,HbT);
if strcmpi(opt_method,'fminsearch')
    options = optimset('MaxFunEvals',5000);
    params = fminsearch(fun,[t0, tau1, tau2, A, B],options);
elseif strcmpi(opt_method,'fminunc')
   options = optimoptions("fminunc",'MaxFunctionEvaluations',2000);
   [params, fval] = fminunc(fun,[t0, tau1, tau2, A, B],options);
else
    return
end
[hrf, ~, ~] = modified_alpha_hrf_2(params(1),params(2),params(3),params(4),params(5),sr,hrf_l);
% sigmoid_GRAB = (1-exp(-params(6)*GRAB))./(1+exp(-params(6)*GRAB));
% sigmoid_GRAB = 1./(1+exp(-params(6)*GRAB));
% params
% sigmoid_GRAB = 1./(1+exp(-(params(6)*GRAB+params(7))));
% sigmoid_GRAB = sigmoid_GRAB/mean(sigmoid_GRAB);
% modulated_jRGECO = sigmoid_GRAB.*jRGECO;
% J = norm(y - X_mat*hrf)^2;
conv_result = conv(jRGECO,hrf);
conv_result = conv_result(1:length(jRGECO));