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

seed_size = 5; % size of the seed
fontsize =20; % font size on figures
ratio = [2.5 1 1]; % figures aspect to save images
LineWidth = 2; % width of the line in plots

%  load a single image to display location of the seed
%  load a single image to display location of the seed
img = dataIn(3).template;
img = img(:,:,1);
figure; colormap('gray'); imagesc(img);


%select brain mask
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
[nX,nY,nT] = size(HbT);
brainIdx = find(brain_mask == 1);

yAll = reshape(HbT(:,:,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
HbT_new = zeros(nY*nX,nT);
HbT_new(brainIdx,:) = ynew';
HbT = reshape(HbT_new,[nX nY nT]);

yAll = reshape(jRGECO(:,:,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
jRGECO_new = zeros(nY*nX,nT);
jRGECO_new(brainIdx,:) = ynew';
jRGECO = reshape(jRGECO_new,[nX nY nT]);

% Experimental parameters for modified alpha function
t0 =    0.1774;
tau1 =  0.4289;
tau2 =  0.4279;
A =     -805.5;
B =     808.3;

sr = 10; % frequency of the HbT, ACh, and jRGECO
hrf_l = 5; % Length of the IRF
t_hrf = (0:floor(hrf_l*sr)-1)/sr; % time points
ts = [0:size(jRGECO,3)-1]/sr;
time_display = 100;
use_time_display = true;
if use_time_display
    [~,t_end] = min(abs(ts-time_display));
else
    t_end = length(ts);
end

[hrf, hrf1, hrf2] = modified_alpha_hrf_2(t0,tau1,tau2,A,B,sr,hrf_l);
figure; plot(t_hrf,hrf); 

pts = [150 90;
    250 90];

% pts = [120 150;
%     240 150];

for u = 1:size(pts,1)
    x_pt = pts(u,1);
    y_pt = pts(u,2);
    
    jRGECO_smoothed = smooth(squeeze(mean(mean(jRGECO(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*5,1),2)));
    HbTavg_smoothed = smooth(squeeze(mean(mean(HbT(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*10^6,1),2)));
    GRAB_smoothed = smooth(squeeze(mean(mean(GRAB(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*5,1),2)));
    
    % plot and save
    FigH = figure('Position', get(0, 'Screensize'));
    yyaxis left; plot(ts(1:t_end), HbTavg_smoothed(1:t_end), 'LineWidth',LineWidth); 
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)', 'fontsize',fontsize);
    ylim([-5 5]);
    hold on; yyaxis right; plot (ts(1:t_end), jRGECO_smoothed(1:t_end),'color','r','LineWidth',LineWidth); 
    ylabel('\Delta F/F','fontsize',fontsize);
    ylim([-0.5 0.5]);
    hold off
    legend('HbT-timeseries','jRGECO1a-timeseries','fontsize',fontsize);
    title('HbT time series vs jRGECO1a time series','Color','red','FontSize',fontsize);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    saveas(FigH, [folder_to_save filesep 'HbT vs jRGECO1a timeseries seed-' num2str(u) '.png']);
    
    % predict HbT with experimental HRF
    conv_result = conv(jRGECO_smoothed,hrf);
    conv_result = conv_result(1:length(jRGECO_smoothed));
    
    % plot and save
    FigH = figure('Position', get(0, 'Screensize'));
    plot(ts(1:t_end), HbTavg_smoothed(1:t_end), 'LineWidth',LineWidth); 
    hold on;  plot (ts(1:t_end), conv_result(1:t_end),'color','r','LineWidth',LineWidth); 
    RMSE = sqrt(mean((HbTavg_smoothed - conv_result).^2));
    R = corrcoef(HbTavg_smoothed,conv_result);
    R = R(2);
    title(['HbT vs predicted HbT - experimental IRF  ' 'RMSE=' num2str(RMSE) '    CorrCoef=' num2str(R)],'Color','red','FontSize',fontsize);
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)','fontsize',fontsize);
    legend('HbT','Predicted-HbT','fontsize',fontsize);
    % plot (GCaMPcorr_smoothed*10, 'LineWidth',2); 
    hold off;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    saveas(FigH, [folder_to_save filesep 'HbT vs predicted HbT IRF seed-' num2str(u) '.png']);
    
    [conv_result, params] = optimise_hrf(A,B,t0,tau1,tau2,sr,hrf_l,jRGECO_smoothed,HbTavg_smoothed,'fminunc');
    FigH = figure('Position', get(0, 'Screensize'));
    plot(ts(1:t_end), HbTavg_smoothed(1:t_end), 'LineWidth',LineWidth); 
    hold on;  plot (ts(1:t_end), conv_result(1:t_end),'color','r','LineWidth',LineWidth); 
    RMSE = sqrt(mean((HbTavg_smoothed - conv_result).^2));
    R = corrcoef(HbTavg_smoothed,conv_result);
    R = R(2);
    title(['HbT vs predicted HbT - experimental IRF -  nonlinear opt  jRGECO' 'RMSE=' num2str(RMSE) '    CorrCoef=' num2str(R)],'Color','red','FontSize',fontsize);
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)','fontsize',fontsize);
    legend('HbT','Predicted-HbT','fontsize',fontsize);
    % plot (GCaMPcorr_smoothed*10, 'LineWidth',2); 
    hold off;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    saveas(FigH, [folder_to_save filesep 'HbT vs predicted HbT IRF  nonlinear opt  jRGECO  seed-' num2str(u) '.png']);
    
    % plot updated IRF and save
    [hrf, ~, ~] = modified_alpha_hrf_2(params(1),params(2),params(3),params(4),params(5),sr,hrf_l);
    figure;
    plot(t_hrf,hrf);
    title('IRF nonlinear opt - jRGECO','Color','red','FontSize',fontsize);
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)','fontsize',fontsize);
%     hold off;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    save([folder_to_save filesep 'IRF nonlinear opt - jRGECO seed-' num2str(u) '.mat'],'hrf','t_hrf','params')
    
    HbT_highpass = highpass(HbTavg_smoothed,0.05,7);
    conv_result_highpass = highpass(conv_result,0.05,7);
    FigH = figure('Position', get(0, 'Screensize'));
    plot(ts(1:t_end), HbT_highpass(1:t_end), 'LineWidth',LineWidth); 
    hold on;  plot (ts(1:t_end), conv_result_highpass(1:t_end),'color','r','LineWidth',LineWidth); 
    RMSE = sqrt(mean((HbT_highpass - conv_result_highpass).^2));
    R = corrcoef(HbT_highpass,conv_result_highpass);
    R = R(2);
    title(['HbT vs predicted HbT - experimental IRF -  high pass after nonlinear opt jRGECO' 'RMSE=' num2str(RMSE) '    CorrCoef=' num2str(R)],'Color','red','FontSize',fontsize);
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)','fontsize',fontsize);
    legend('HbT','Predicted-HbT','fontsize',fontsize);
    % plot (GCaMPcorr_smoothed*10, 'LineWidth',2); 
    hold off;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    saveas(FigH, [folder_to_save filesep 'HbT vs predicted HbT IRF  high pass after  nonlinear opt jRGECO  seed-' num2str(u) '.png']);
    
end