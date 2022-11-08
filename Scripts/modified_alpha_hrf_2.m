function [hrf, D, C] = modified_alpha_hrf_2(t0,tau1,tau2,A,B,sr,hrf_l)
% modified_alpha_hrf: generate an HRF based on the sum of 2 modified alpha
% functions
%   written by Harrison Fisher 
%   
%   Arguments: 
%       t0: delay time in s
%       tau1: shape parameter for the positive alpha function
%       tau2: shape parameter for the negative alpha function
%       A:  amplitude parameter for the positive alpha function
%       B:  amplitude parameter for the negative alpha function
%       hrf_l: length of the HRF in seconds 

    N = floor(hrf_l*sr);
    dt=1/sr;
    t_hrf=((1:N)/sr)'; % make sure the output will be a column

%     tr=t_hrf(round(t0/dt):length(t_hrf))-(t0-dt);
% %     tr = ((0:N-1)/sr)-t0;
% 
%     D = zeros(size(t_hrf));
%     D(round(t0/dt):length(t_hrf)) = ((tr)./tau1).^3 .* exp(-(tr)./tau1);
% 
%     C = zeros(size(t_hrf));
%     C(round(t0/dt):length(t_hrf)) = ((tr)./tau2).^3 .* exp(-(tr)./tau2);

%     tr=t_hrf(round(t0/dt):length(t_hrf))-(t0-dt);
    tr = (((0:N-1)/sr)-t0)';

%     D = zeros(size(t_hrf));
    D = ((tr)./tau1).^3 .* exp(-(tr)./tau1);

%     C = zeros(size(t_hrf));
    C = ((tr)./tau2).^3 .* exp(-(tr)./tau2);

    hrf = A*D + B*C;
end
