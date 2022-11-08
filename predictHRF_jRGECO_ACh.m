% load hemodynamics, ACh, and jRGECOa1 data
load('run0002_greenFP.mat', 'led_475_sm_norm_HD')
GRAB = led_475_sm_norm_HD;
load('run0002_hemodynamics.mat', 'HbT')
load('run0002_redFP.mat', 'led_590_sm_norm_HD')
jRGECO = led_590_sm_norm_HD;

% Save results in this folder
folder_to_save = 'C:\Users\Sreekanth\Documents\Projects\IRF\Predict HRF\Results\jRGECOonly\';
if ~isfolder(folder_to_save)
    mkdir(folder_to_save); 
end

seed_size = 5; % size of the seed
fontsize =20; % font size on figures
ratio = [2.5 1 1]; % figures aspect to save images
LineWidth = 2; % width of the line in plots

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
ts = [0:size(jRGECO,3)-1]/sr;

[hrf, hrf1, hrf2] = modified_alpha_hrf_2(t0,tau1,tau2,A,B,sr,hrf_l);
figure; plot(t_hrf,hrf); 

pts = [87 100;
    30 100];

for u = 1:size(pts,1)
    x_pt = pts(u,1);
    y_pt = pts(u,2);
    
    jRGECO_smoothed = smooth(squeeze(mean(mean(jRGECO(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*5,1),2)));
    HbTavg_smoothed = smooth(squeeze(mean(mean(HbT(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*10^6,1),2)));
    GRAB_smoothed = smooth(squeeze(mean(mean(GRAB(x_pt-seed_size:x_pt+seed_size,y_pt-seed_size:y_pt+seed_size,:)*5,1),2)));
    
    % plot and save
    FigH = figure('Position', get(0, 'Screensize'));
    yyaxis left; plot(ts, HbTavg_smoothed, 'LineWidth',LineWidth); 
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)', 'fontsize',fontsize);
    ylim([-5 5]);
    hold on; yyaxis right; plot (ts, jRGECO_smoothed,'color','r','LineWidth',LineWidth); 
    ylabel('\Delta F/F','fontsize',fontsize);
    ylim([-0.5 0.5]);
    hold off
    legend('HbT-timeseries','jRGECO1a-timeseries','fontsize',fontsize);
    title('HbT time series vs jRGECO1a time series','Color','red','FontSize',fontsize);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    saveas(FigH, [folder_to_save 'HbT vs jRGECO1a timeseries seed-' num2str(u) '.png']);
    
    % predict HbT with experimental HRF
    conv_result = conv(jRGECO_smoothed,hrf);
    conv_result = conv_result(1:length(jRGECO_smoothed));
    
    % plot and save
    FigH = figure('Position', get(0, 'Screensize'));
    plot(ts, HbTavg_smoothed, 'LineWidth',LineWidth); 
    hold on;  plot (ts, conv_result,'color','r','LineWidth',LineWidth); 
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
    saveas(FigH, [folder_to_save 'HbT vs predicted HbT IRF seed-' num2str(u) '.png']);
    
    [conv_result, params] = optimise_hrf(A,B,t0,tau1,tau2,sr,hrf_l,jRGECO_smoothed,HbTavg_smoothed,'fminunc');
    FigH = figure('Position', get(0, 'Screensize'));
    plot(ts, HbTavg_smoothed, 'LineWidth',LineWidth); 
    hold on;  plot (ts, conv_result,'color','r','LineWidth',LineWidth); 
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
    saveas(FigH, [folder_to_save 'HbT vs predicted HbT IRF  nonlinear opt  jRGECO  seed-' num2str(u) '.png']);
    
    % plot updated IRF and save
    [hrf, ~, ~] = modified_alpha_hrf_2(params(1),params(2),params(3),params(4),params(5),sr,hrf_l);
    figure;
    plot(t_hrf,hrf);
    title('IRF nonlinear opt - jRGECO','Color','red','FontSize',fontsize);
    xlabel('t(s)','fontsize',fontsize)
    ylabel('\Delta C (\muM)','fontsize',fontsize);
    hold off;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize)
    set(gca,'LooseInset',get(gca,'TightInset'));
    pbaspect(ratio)
    save([folder_to_save 'IRF nonlinear opt - jRGECO seed-' num2str(u) '.mat'],'hrf','t_hrf','params')
    
    HbT_highpass = highpass(HbTavg_smoothed,0.05,7);
    conv_result_highpass = highpass(conv_result,0.05,7);
    FigH = figure('Position', get(0, 'Screensize'));
    plot(ts, HbT_highpass, 'LineWidth',LineWidth); 
    hold on;  plot (ts, conv_result_highpass,'color','r','LineWidth',LineWidth); 
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
    saveas(FigH, [folder_to_save 'HbT vs predicted HbT IRF  high pass after  nonlinear opt jRGECO  seed-' num2str(u) '.png']);
    
end