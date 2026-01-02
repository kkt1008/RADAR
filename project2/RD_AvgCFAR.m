%% Data Set 불러오기 & Radar Information
clc; clear; close all;

% Chirp Config
nrx               = 4;
startFreq         = 77;   
idleTime          = 100;  
adcStartTime      = 6;    
rampEndTime       = 60;   
digOutSampleRate  = 4270; 
freqSlopeConst    = 66.143;
numADCSamples     = 256;  

% Frame Config
numLoops          = 128;
framePeriodicity  = 40.0; 

% Radar Parameters
c0            = 299792458;
nr            = numADCSamples;
nd            = numLoops;
fs            = digOutSampleRate * 1e3;
kf            = freqSlopeConst * 1e12;
adc_duration  = nr / fs;
PRI           = (idleTime + rampEndTime) * 1e-6;
PRF           = 1/PRI;
Rmax          = fs * c0 / (2 * kf);

% Load Raw Data & reshape
adc_file_name = 'hand_range.bin';  % or 'human_azimuth.bin'
data = readDCA1000(adc_file_name);
disp('data was generated!');

data_1d   = data(2, :);
numTotal  = numel(data_1d);
numFrames = floor(numTotal / (nr * nd));
data_1d   = data_1d(1 : nr * nd * numFrames);
data_2d   = reshape(data_1d, [nr, nd * numFrames]);
data_3d   = reshape(data_2d, [nr, nd, numFrames]);
[nr, nd, nf] = size(data_3d);

%% 1D FFT Spectrogram + Average Threshold + α·σ
% compute 1D FFT magnitude per chirp
mag1D = abs( fft(data_2d, [], 1) );    % [nr × nchirps]
Time_Index_All = (0:size(mag1D,2)-1) * PRI;
Range_Index    = (0:nr-1)/nr * Rmax;

% α = 0
alpha0    = 0;
thr1D_0   = mean(mag1D,1) + alpha0 * std(mag1D,[],1);
resid1D_0 = mag1D - thr1D_0;
resid1D_0(resid1D_0<0) = 0;
Residual1D_AvgSigma0_dB = 20*log10(resid1D_0 + eps);

% α = 2
alpha2    = 2;
thr1D_2   = mean(mag1D,1) + alpha2 * std(mag1D,[],1);
resid1D_2 = mag1D - thr1D_2;
resid1D_2(resid1D_2<0) = 0;
Residual1D_AvgSigma2_dB = 20*log10(resid1D_2 + eps);

% Plot & save
figure;
imagesc(Time_Index_All, Range_Index, Residual1D_AvgSigma0_dB);
axis xy tight; colormap(parula);
xlabel('Time (s)'); ylabel('Range (m)');
title('1D FFT Residual (Avg+0·σ)');
saveas(gcf,'1D_Residual_AvgSigma0.png');

figure;
imagesc(Time_Index_All, Range_Index, Residual1D_AvgSigma2_dB);
axis xy tight; colormap(parula);
xlabel('Time (s)'); ylabel('Range (m)');
title('1D FFT Residual (Avg+2·σ)');
saveas(gcf,'1D_Residual_AvgSigma2.png');

%% 2D RD_map + Average Threshold + α·σ
% compute RD_map
fftRange   = fft(data_3d, [], 1);
fftDoppler = fft(fftRange, [], 2);
RD_map     = abs( fftshift(fftDoppler, 2) );
RD_map(:,floor(nd/2)+1,:) = 0;
RD_map(255:end,:,:) = 0;

% allocate
ResidualRD_AvgSigma0 = zeros(nr,nd,nf);
ResidualRD_AvgSigma2 = zeros(nr,nd,nf);

for f = 1:nf
    M = RD_map(:,:,f);
    mu0    = mean(M(:)); sigma0    = std(M(:));
    thrRD0 = mu0 + alpha0*sigma0;
    R0     = M - thrRD0; R0(R0<0)=0;
    ResidualRD_AvgSigma0(:,:,f) = R0;
    mu2    = mu0; sigma2    = sigma0;
    thrRD2 = mu2 + alpha2*sigma2;
    R2     = M - thrRD2; R2(R2<0)=0;
    ResidualRD_AvgSigma2(:,:,f) = R2;
end

% save variables
RD_AvgSigma0_Residual    = ResidualRD_AvgSigma0;    %#ok<NASGU>
RD_AvgSigma2_Residual    = ResidualRD_AvgSigma2;    %#ok<NASGU>

% animate & save for α=0
vRD0 = VideoWriter('2D_RD_AvgSigma0.avi','Motion JPEG AVI'); vRD0.FrameRate=10; open(vRD0);
figure; hImg2 = imagesc(1:nd,1:nr,ResidualRD_AvgSigma0(:,:,1));
axis xy tight; colormap(jet); colorbar; caxis([0 max(RD_map(:))]);
xlabel('Chirp Index'); ylabel('Range Bin');
title('RD Residual (Avg+0·σ) Frame 1'); drawnow; writeVideo(vRD0,getframe(gcf));
for k=2:nf
    set(hImg2,'CData',ResidualRD_AvgSigma0(:,:,k));
    title(sprintf('RD Residual (Avg+0·σ) Frame %d',k)); drawnow;
    writeVideo(vRD0,getframe(gcf));
end
close(vRD0);

% animate & save for α=2
vRD2 = VideoWriter('2D_RD_AvgSigma2.avi','Motion JPEG AVI'); vRD2.FrameRate=10; open(vRD2);
figure; hImg3 = imagesc(1:nd,1:nr,ResidualRD_AvgSigma2(:,:,1));
axis xy tight; colormap(jet); colorbar; caxis([0 max(RD_map(:))]);
xlabel('Chirp Index'); ylabel('Range Bin');
title('RD Residual (Avg+2·σ) Frame 1'); drawnow; writeVideo(vRD2,getframe(gcf));
for k=2:nf
    set(hImg3,'CData',ResidualRD_AvgSigma2(:,:,k));
    title(sprintf('RD Residual (Avg+2·σ) Frame %d',k)); drawnow;
    writeVideo(vRD2,getframe(gcf));
end
close(vRD2);

%% 3D 애니메이션 예시 (α=2)
[X,Y] = meshgrid(1:nd,1:nr);
v3D = VideoWriter('3D_RD_AvgSigma2.avi','Motion JPEG AVI'); v3D.FrameRate=10; open(v3D);
figure; hS = surf(X,Y,ResidualRD_AvgSigma2(:,:,1),'EdgeColor','none'); colormap(jet); colorbar;
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Residual');
title('3D RD Residual (Avg+2·σ) Frame 1'); view(45,30); drawnow; writeVideo(v3D,getframe(gcf));
for k=2:nf
    set(hS,'ZData',ResidualRD_AvgSigma2(:,:,k));
    title(sprintf('3D RD Residual (Avg+2·σ) Frame %d',k)); drawnow;
    writeVideo(v3D,getframe(gcf));
end
close(v3D);
