%% Data Set 불러오기 & Radar Information
clc; clear; close all;

% Chirp Config
nrx               = 4;    % Number of rx Antenna
startFreq         = 77;   % [GHz]
idleTime          = 100;  % [µs]
adcStartTime      = 6;    % [µs]
rampEndTime       = 60;   % [µs]
digOutSampleRate  = 4270; % [kHz]
freqSlopeConst    = 66.143; % [MHz/µs]
numADCSamples     = 256;  % Samples per chirp

% Frame Config
chirpStartIdx     = 0;
chirpEndIdx       = 1;
numLoops          = 128;  % Chirps per frame
framePeriodicity  = 40.0; % [ms]

% ADC Config
numADCBits = 16;
isReal     = 1;

% Radar Parameters
c0            = 299792458;
nr            = numADCSamples;
nd            = numLoops;
fs            = digOutSampleRate * 1e3;    % [Hz]
kf            = freqSlopeConst   * 1e12;   % [Hz/s]
adc_duration  = nr / fs;
BW            = adc_duration * kf;
PRI           = (idleTime + rampEndTime)*1e-6;
PRF           = 1/PRI;
Rmax          = fs * c0 / (2*kf);

% Load Raw Data & reshape
adc_file_name = 'hand_range.bin';  % or 'human_azimuth.bin'
data = readDCA1000(adc_file_name);
disp('data was generated!');

data_1d = data(2, :);
numTotal = numel(data_1d);
numFrames = floor(numTotal / (nr*nd));
data_1d = data_1d(1 : nr*nd*numFrames);
data_2d = reshape(data_1d, [nr, nd*numFrames]);
data_3d = reshape(data_2d, [nr, nd, numFrames]);
[nr, nd, nf] = size(data_3d);

%% 2D Static CFAR Thresholding on data_3d
blockSize = 8;
alpha     = 2;
thr2D_map = zeros(nr, nd, nf);

% compute RD_map
fftRange   = fft(data_3d, [], 1);
fftDoppler = fft(fftRange, [], 2);
RD_map     = abs( fftshift(fftDoppler, 2) );
RD_map(:, floor(nd/2)+1, :) = 0;
RD_map(255:end,:,:) = 0;

for f = 1:nf
    RD = RD_map(:,:,f);
    for br = 1:(nr/blockSize)
        for bc = 1:(nd/blockSize)
            rows = (br-1)*blockSize+1 : br*blockSize;
            cols = (bc-1)*blockSize+1 : bc*blockSize;
            blk  = RD(rows,cols);
            thr2D_map(rows,cols,f) = alpha * mean(blk(:));
        end
    end
end

residAll = RD_map - thr2D_map;
residAll(residAll < 0) = 0;

[X, Y]      = meshgrid(1:nd, 1:nr);
frameToPlot = 1;
nFrames     = nf;
maxResid    = max(residAll(:));

%% Figure 1: 3D RD Map
figure(1);
surf(X, Y, RD_map(:,:,frameToPlot), 'EdgeColor','none','FaceAlpha',0.8);
colormap(parula); colorbar; view(45,30);
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Amplitude');
title(sprintf('3D RD Map (Frame %d)', frameToPlot));

%% Figure 2: 3D Static CFAR Threshold
figure(2);
surf(X, Y, thr2D_map(:,:,frameToPlot), 'EdgeColor','none','FaceAlpha',0.9);
shading interp; colormap(parula); colorbar; view(45,30);
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Threshold');
title(sprintf('3D Static CFAR Threshold (Frame %d)', frameToPlot));

%% Figure 3: RD Map vs Static Threshold
figure(3); clf;
subplot(1,2,1);
surf(X, Y, RD_map(:,:,frameToPlot), 'EdgeColor','none','FaceAlpha',0.8);
colormap(parula); colorbar; view(45,30);
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Amplitude');
title('RD Map');

ax2 = subplot(1,2,2);
surf(X, Y, thr2D_map(:,:,frameToPlot), 'EdgeColor','none','FaceAlpha',0.9);
shading interp; colormap(parula); colorbar;
caxis([0, max(thr2D_map(:,:,frameToPlot),[],'all')]); view(45,30);
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Threshold');
title('Static CFAR Threshold');
set(ax2,'Color',[0.85 0.9 1]);

% static residual & video 저장용 변수 설정
nSave = nFrames;
RD_StaticCFAR_2D_Residual      = residAll(:,:,1:nSave);
RD_StaticCFAR_Target_Tracking  = zeros(nr, nd, nSave);

%% Figure 4: 2D Residual (Static CFAR) 애니메이션 + 저장
pauseTime4 = 0.033;
v4 = VideoWriter('2D_Residual_StaticCFAR.avi','Motion JPEG AVI'); v4.FrameRate = 10; open(v4);
figure(4); clf;
hImg = imagesc(1:nd, 1:nr, residAll(:,:,1));
axis xy tight; colormap(parula); colorbar; caxis([0, maxResid]);
xlabel('Chirp Index'); ylabel('Range Bin');
title('2D Residual (Static CFAR, Frame 1)');
drawnow; writeVideo(v4, getframe(gcf)); pause(pauseTime4);

for k = 2:nSave
    set(hImg, 'CData', residAll(:,:,k));
    title(sprintf('2D Residual (Static CFAR, Frame %d)', k));
    drawnow; writeVideo(v4, getframe(gcf)); pause(pauseTime4);
end
close(v4);

%% Figure 5: 3D RD_StaticCFAR_2D_Residual 애니메이션 + 저장
v5 = VideoWriter('3D_RD_StaticCFAR_2D_Residual.avi','Motion JPEG AVI'); v5.FrameRate = 10; open(v5);
figure(5); clf;
hSurf = surf(X, Y, residAll(:,:,1), 'EdgeColor','none','FaceAlpha',0.8);
colormap(parula); colorbar;
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Amplitude');
title('3D Residual Static CFAR (Frame 1)'); view(45,30);
drawnow; writeVideo(v5, getframe(gcf)); pause(pauseTime4);

for k = 2:nSave
    set(hSurf, 'ZData', residAll(:,:,k));
    title(sprintf('3D Residual Static CFAR (Frame %d)', k));
    drawnow; writeVideo(v5, getframe(gcf)); pause(pauseTime4);
end
close(v5);

%% Figure 6: Static Residual 기반 Target Tracking Enhancement + 저장
window_size = 15; half_w = floor(window_size/2); maxJump = window_size;
% Center-symmetry kernel
weight_kernel = zeros(window_size);
for i = 1:window_size
    for j = 1:window_size
        weight_kernel(i,j) = 1/(1 + abs(i-(half_w+1)) + abs(j-(half_w+1)));
    end
end
pauseTime6 = 0.1;
v6 = VideoWriter('StaticCFAR_Target_Tracking.avi','Motion JPEG AVI'); v6.FrameRate = 10; open(v6);
figure(6); clf;
firstFrame = residAll(:,:,1);
avg0       = conv2(firstFrame, weight_kernel, 'valid');
[~, idx0]  = max(avg0(:));
[r0, d0]   = ind2sub(size(avg0), idx0);
prev_r     = r0 + half_w;
prev_d     = d0 + half_w;
target_weight = 1;

for frameIdx = 1:nSave
    RD_Frame = residAll(:,:,frameIdx);
    avg_RD   = conv2(RD_Frame, weight_kernel, 'valid');
    [~, idx] = max(avg_RD(:));
    [mr, md]  = ind2sub(size(avg_RD), idx);
    target_r = mr + half_w;
    target_d = md + half_w;
    % noise suppression
    dr = target_r - prev_r; dd = target_d - prev_d;
    if sqrt(dr^2 + dd^2) > maxJump
        dr = sign(dr)*maxJump; dd = sign(dd)*maxJump;
        target_r = prev_r + dr; target_d = prev_d + dd;
    end
    % weight update
    if abs(target_r-prev_r)<=half_w && abs(target_d-prev_d)<=half_w
        target_weight = target_weight + 1;
    else
        target_weight = 1;
    end
    Tker = weight_kernel * target_weight;
    % apply mask & save
    r_rng = (target_r-half_w):(target_r+half_w);
    d_rng = (target_d-half_w):(target_d+half_w);
    RD_Frame(r_rng,d_rng) = RD_Frame(r_rng,d_rng) .* Tker;
    RD_StaticCFAR_Target_Tracking(:,:,frameIdx) = RD_Frame;
    % animate & record
    imagesc(1:nd,1:nr,RD_Frame);
    axis xy tight; colormap(parula); colorbar; caxis([0, maxResid]);
    xlabel('Chirp Index'); ylabel('Range Bin');
    title(sprintf('Tracking Frame %d', frameIdx));
    hold on;
    rectangle('Position',[d_rng(1),r_rng(1),window_size,window_size],...
              'EdgeColor','r','LineWidth',1.5);
    text(target_d, target_r-half_w-2, sprintf('Wt=%d',target_weight),...
         'Color','r','FontSize',7,'FontWeight','bold',...
         'HorizontalAlignment','center');
    hold off;
    drawnow; writeVideo(v6, getframe(gcf)); pause(pauseTime6);
    prev_r = target_r; prev_d = target_d;
end
close(v6);

%% Figure 7: 3D RD_StaticCFAR_Target_Tracking 애니메이션 + 저장
v7 = VideoWriter('3D_RD_StaticCFAR_Target_Tracking.avi','Motion JPEG AVI'); v7.FrameRate = 10; open(v7);
figure(7); clf;
hSurf2 = surf(X, Y, RD_StaticCFAR_Target_Tracking(:,:,1), 'EdgeColor','none','FaceAlpha',0.8);
colormap(parula); colorbar;
xlabel('Chirp Index'); ylabel('Range Bin'); zlabel('Residual');
title('3D Tracking Static CFAR (Frame 1)'); view(45,30);
drawnow; writeVideo(v7, getframe(gcf)); pause(pauseTime4);

for k = 2:nSave
    set(hSurf2, 'ZData', RD_StaticCFAR_Target_Tracking(:,:,k));
    title(sprintf('3D Tracking Static CFAR (Frame %d)', k));
    drawnow; writeVideo(v7, getframe(gcf)); pause(pauseTime4);
end
close(v7);
