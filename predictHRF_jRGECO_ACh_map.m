clear all
% load hemodynamics, ACh, and jRGECOa1 data
load('run0002_greenFP.mat', 'led_475_sm_norm_HD')
GRAB = led_475_sm_norm_HD;
load('run0002_hemodynamics.mat', 'HbT')
load('run0002_redFP.mat', 'led_590_sm_norm_HD')
jRGECO = led_590_sm_norm_HD;

% Save results in this folder
folder_to_save = 'C:\Users\Sreekanth\Documents\Projects\IRF\Predict HRF\Results\jRGECOonlyMaP\';
if ~isfolder(folder_to_save)
    mkdir(folder_to_save); 
end

seed_size = 0; % size of the seed
fontsize =20; % font size on figures
ratio = [2.5 1 1]; % figures aspect to save images

%  load a single image to display location of the seed
img = imread('C:\Users\Sreekanth\Documents\Projects\IRF\Data\Meso_211112_NM07\run02_X1002.tif');
figure; colormap('gray'); imagesc(img);

% Experimental parameters for modified alpha function
t0 =    0.1774;
tau1 =  0.4289;
tau2 =  0.4279;
A =     -805.5;
B =     808.3;

sr = 7; % frequency of the HbT, ACh, and jRGECO
hrf_l = 5; % Length of the IRF
t_hrf = (0:floor(hrf_l*sr)-1)/sr; % time points

% Generate experimenatl IRF
[hrf, hrf1, hrf2] = modified_alpha_hrf_2(t0,tau1,tau2,A,B,sr,hrf_l);
figure; plot(t_hrf,hrf); 

size_x = size(HbT,1);
size_y = size(HbT,2);

% store correlation map and moving correlation volumes
img1 = zeros(size_x,size_y);
img2 = zeros(size_x,size_y);
img3 = zeros(size_x,size_y);
img4 = zeros(size_x,size_y);
volume1 = zeros(size_x,size_y,15);
volume2 = zeros(size_x,size_y,10);
volume3 = zeros(size_x,size_y,5);
volume4 = zeros(size_x,size_y,15);
parameters = zeros(size_x,size_y,5);

for u = 1:size_x
    u
    for v = 1:size_y
        
        x_pt = u;
        y_pt = v;
        
        % temporal smooth
        jRGECO_smoothed = smooth(squeeze(mean(mean(jRGECO(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*5,1),2)));
        HbTavg_smoothed = smooth(squeeze(mean(mean(HbT(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*10^6,1),2)));
        GRAB_smoothed = smooth(squeeze(mean(mean(GRAB(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*5,1),2)));
        
        % outside the brain region all the values are NaN
        if ~any(isnan(HbTavg_smoothed))
            
            % predict HbT with experimental HRF
            conv_result = conv(jRGECO_smoothed,hrf);
            conv_result = conv_result(1:length(jRGECO_smoothed));
            
            % find correlation between HbT and predicted HbT
            R = corrcoef(HbTavg_smoothed,conv_result);
            R = R(2);
            img1(u,v) = R;
            
            % Find moving correlation between HbT and predcited HbT
            time_step = length(HbTavg_smoothed)/85;
            for w = 1:85
                t_start = floor(1+time_step*(w-1));
                t_end = min(floor(t_start+sr*30),length(HbTavg_smoothed));
                R = corrcoef(HbTavg_smoothed(t_start:t_end),conv_result(t_start:t_end));
                R = R(2);
                volume1(u,v,w) = R;
            end
            
            % Apply high pass filter to HbT and predcited HbT and get
            % correlation and moving correlation
            HbT_highpass = highpass(HbTavg_smoothed,0.05,7);
            conv_result_highpass = highpass(conv_result,0.05,7);
            R = corrcoef(HbT_highpass,conv_result_highpass);
            R = R(2);
            img2(u,v) = R;
            
            for w = 1:85
                t_start = floor(1+time_step*(w-1));
                t_end = min(floor(t_start+sr*30),length(HbTavg_smoothed));
                R = corrcoef(HbT_highpass(t_start:t_end),conv_result(t_start:t_end));
                R = R(2);
                volume2(u,v,w) = R;
            end
           
            % Apply non linear optimisation get updated parameters for IRF.
            % Start with experimental HRF and fit the parameters using HbT
            % and jRGECO.
            % With updated IRF get predicted HbT
            [conv_result, params] = optimise_hrf(A,B,t0,tau1,tau2,sr,hrf_l,jRGECO_smoothed,HbTavg_smoothed,'fminunc');
            parameters(u,v,:) = params;
            
            % find correlation between HbT and predicted HbT
            R = corrcoef(HbTavg_smoothed,conv_result);
            R = R(2);
            img3(u,v) = R;
            
            % Find moving correlation between HbT and predcited HbT
            time_step = length(HbTavg_smoothed)/85;
            for w = 1:85
                t_start = floor(1+time_step*(w-1));
                t_end = min(floor(t_start+sr*30),length(HbTavg_smoothed));
                R = corrcoef(HbTavg_smoothed(t_start:t_end),conv_result(t_start:t_end));
                R = R(2);
                volume3(u,v,w) = R;
            end
            
            % Apply high pass filter to HbT and predcited HbT and get
            % correlation and moving correlation
            conv_result_highpass = highpass(conv_result,0.05,7);
            R = corrcoef(HbT_highpass,conv_result_highpass);
            R = R(2);
            img4(u,v) = R;
            
            for w = 1:85
                t_start = floor(1+time_step*(w-1));
                t_end = min(floor(t_start+sr*30),length(HbTavg_smoothed));
                R = corrcoef(HbT_highpass(t_start:t_end),conv_result(t_start:t_end));
                R = R(2);
                volume4(u,v,w) = R;
            end

        end
    end
end
save([folder_to_save filesep 'corr_results.mat'],'img1','img2','img3','img4','volume1','volume2','volume3','volume4','parameters');




