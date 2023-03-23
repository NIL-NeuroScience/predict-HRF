clear all
root = '/projectnb/devorlab/nfominth/2023/mesoscope/23-01-25/Thy1_0110/processed/';
root_trigger = '/projectnb/devorlab/nfominth/2023/mesoscope/23-01-25/Thy1_0110/';
data_in = '/projectnb/devorlab/nfominth/2023/mesoscope/23-01-25/Thy1_0110/';

Run = 3;
rfp_norm=h5read(sprintf('%s%srun%04i.h5',root,filesep,Run),'/rfp/norm');
gfp_normHD=h5read(sprintf('%s%srun%04i.h5',root,filesep,Run),'/gfp/normHD');
HbO=h5read(sprintf('%s%srun%04i.h5',root,filesep,Run),'/hemodynamics/HbO');
HbR=h5read(sprintf('%s%srun%04i.h5',root,filesep,Run),'/hemodynamics/Hb');
HbT=HbR+HbO;
load([root_trigger filesep 'Triggers' filesep sprintf('Run%03i.mat',Run)]);
load([data_in 'dataIn.mat']); 
jRGECO= rfp_norm;
GRAB = gfp_normHD;

% Save results in this folder
folder_to_save = '/projectnb/devorlab/skura/HRF/Analysis/run_3_230125';
if ~isfolder(folder_to_save)
    mkdir(folder_to_save); 
end

seed_size = 0; % size of the seed
fontsize =20; % font size on figures
ratio = [2.5 1 1]; % figures aspect to save images

%  load a single image to display location of the seed
img = dataIn(3).template;
img = img(:,:,1);

% Detrend GRAB and Ach signals. They have some decay in signal, in  the
% future this will be fixed on the acquistion side
figure; colormap('gray'); imagesc(img);
for u = 1:size(GRAB,1)
    for v = 1:size(GRAB,2)
        GRAB(u,v,:) = detrend(squeeze(GRAB(u,v,:)),2);
    end
end

for u = 1:size(jRGECO,1)
    for v = 1:size(jRGECO,2)
        jRGECO(u,v,:) = detrend(squeeze(jRGECO(u,v,:)),2);
    end
end


%select brain mask. 
brain_mask = zeros(size(img));
h = msgbox(' Please select one side of the brain region');
uiwait(h)
figure; colormap('gray');
while(1)
    imagesc(img);
    [brain_mask1,Xi1,Yi1] = roipoly;
    hold on;
    plot(Xi1,Yi1,'color','k');
    hold off;
    button = questdlg('Are you statisfied with ROI');
    if strcmp(button, 'Yes')
      break;
    end
end
h = msgbox(' Please select other side of the brain region');
uiwait(h)
while(1)
    imagesc(img); hold on;
    plot(Xi1,Yi1,'color','k');
    hold off;
    [brain_mask2,Xi2,Yi2] = roipoly;
    hold on;
    plot(Xi2,Yi2,'color','k');
    hold off;
    button = questdlg('Are you statisfied with ROI');
    if strcmp(button, 'Yes')
      break;
    end
end
imagesc(img); hold on;
plot(Xi1,Yi1,'color','k');
plot(Xi2,Yi2,'color','k');
hold off;
brain_mask(brain_mask1==1) = 1;
brain_mask(brain_mask2==1) = 1;

save([folder_to_save filesep 'brain_mask.mat'],'brain_mask1','Xi1','Yi1','brain_mask2','Xi2','Yi2','brain_mask')

% Apply global mean regression. I am commenting this just in case if we
% want to use this.
% [nX,nY,nT] = size(HbT);
% brainIdx = find(brain_mask == 1);
% 
% yAll = reshape(HbT(:,:,:),[nY*nX,nT]);
% y1 = yAll(brainIdx,:);
% yMean = mean(y1,1);
% a = pinv(yMean*yMean')*yMean*y1';
% ynew = y1'-yMean'*a;
% HbT_new = zeros(nY*nX,nT);
% HbT_new(brainIdx,:) = ynew';
% HbT = reshape(HbT_new,[nX nY nT]);
% 
% yAll = reshape(jRGECO(:,:,:),[nY*nX,nT]);
% y1 = yAll(brainIdx,:);
% yMean = mean(y1,1);
% a = pinv(yMean*yMean')*yMean*y1';
% ynew = y1'-yMean'*a;
% jRGECO_new = zeros(nY*nX,nT);
% jRGECO_new(brainIdx,:) = ynew';
% jRGECO = reshape(jRGECO_new,[nX nY nT]);
for u = 1:size(HbT,3)
    HbT(:,:,u) = HbT(:,:,u).*brain_mask;
end

for u = 1:size(jRGECO,3)
    jRGECO(:,:,u) = jRGECO(:,:,u).*brain_mask;
end
% Experimental parameters for modified alpha function
t0 =    0.1774;
tau1 =  0.4289;
tau2 =  0.4279;
A =     -805.5;
B =     808.3;

sr = 10; % frequency of the HbT, ACh, and jRGECO
hrf_l = 5; % Length of the IRF
t_hrf = (0:floor(hrf_l*sr)-1)/sr; % time points

% Generate experimenatl IRF
[hrf, hrf1, hrf2] = modified_alpha_hrf_2(t0,tau1,tau2,A,B,sr,hrf_l);
figure; plot(t_hrf,hrf); 

size_x = round(size(HbT,1)/3);
size_y = round(size(HbT,2)/3);

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
        
        x_pt = 1+3*(u-1);
        y_pt = 1+3*(v-1);
        
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
            
            n_steps = 300;
            % Find moving correlation between HbT and predcited HbT
            time_step = length(HbTavg_smoothed)/n_steps;
            for w = 1:n_steps
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
            
            for w = 1:n_steps
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
            time_step = length(HbTavg_smoothed)/n_steps;
            for w = 1:n_steps
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
            
            for w = 1:n_steps
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




