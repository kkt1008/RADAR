%% Data Set 불러오기 & Radar Information
clc;
clear;
close all;

% Chirp Config
nrx = 4; % Number of rx Antenna
startFreq       =   77; % Start Frequency : 77 [GHz]
idleTime        =   100; % IDLE 시간 [us] : 1 Chirp Sweep 후 IDLE하는 시간
adcStartTime    =   6; % ADC 시작 시간 [us] : Chirp Sweep 시작 후 ADC 시작하는 시간
rampEndTime     =   60; % Chirp 전체 시간 [us]
digOutSampleRate=   4270; % Sampling Frequency(Fs) [kHz]
freqSlopeConst  =   66.143; % 주파수 Sweep 속도(S) [MHz/us]
numADCSamples   = 256; % number of ADC samples per chirp => Number of Samples per Chirp : Range FFT's N

% Frame Config
chirpStartIdx   =   0;
chirpEndIdx     =   1; 
numLoops        =   128; % Number of Chirps per Frame : Doppler FFT's N
framePeriodicity=   40.0; % Frame 주기 : 40 [ms]

% ADC Config
numADCBits      = 16; % number of ADC bits per sample
isReal = 1; % set to 1 if real only data, 0 if complex data0

% Radar Parameters
c0              =   299792458; % 빛의 속도 
nrx             =   1;
ntx             =   chirpEndIdx - chirpStartIdx + 1; % 1 - 0 + 1 = 2
nd              =   numLoops; % Number of Chirps per Frame : N_DopplerFFT = 128
nr              =   numADCSamples; % Number of Samples per Chirp : N_RangeFFT = 256
tfr             =   framePeriodicity * 1e-3; % Frame 주기[s]
fs              =   digOutSampleRate * 1e3; % Sampling Frequency[Hz]
kf              =   freqSlopeConst * 1e12; % 주파수 Sweep 속도 S [Hz/s]
adc_duration    =   nr / fs; % Chirp Duration Tc : Sweep 시간 per 1chirp(실제 Sampling 진행되는 시간 ramp보다 길 수 있음)
BW              =   adc_duration * kf; % Tc * S = BW
PRI             =   (idleTime + rampEndTime) * 1e-6; % Pulse Repetition Interval
PRF             =   1/PRI; % 1초당 chirp 쏘는 횟수 [Hz] -> Doppler Fs
fc              =   startFreq * 1e9 + kf * (adcStartTime * 1e-6 + adc_duration / 2); 
nfr             =   2;
frame_duration  =   PRI*numLoops;
frame_idletime  =   tfr - frame_duration;
Rmax            =   fs * c0/(2*kf);

% Load RawData
%adc_file_name = 'hand_range.bin';
adc_file_name = 'human_azimuth.bin';
data = readDCA1000(adc_file_name);
disp('data was generated!');

% channel 1's All Data Load : (1, 5734400) 행 vector
data_1d=data(2,1:end);  
% channel 1's Fist Frame First Chirp vector (1, 256)
data_sample = data(2,1:256); 

% t=[0:1/fs:(4.0e-5)-1/fs]; % ???

% data1을 (256, length(data1)/256)으로 나눠서 2차원 행렬로 바꿈 (256,22400)
data_2d=reshape(data_1d,[nr,length(data_1d)/nr]);

% (256, 22400) -> (256, 128, 175)로 바꿔야 함
% (256, 128, 175) = (256 Samples per 1Chirp, 128 Chirps, 175 Frames)
data_3d = reshape(data_2d, [nr, 128, 175]); 

[i, j]=size(data_1d); % [i, j] = [1, 5734400]
[y, x]=size(data_2d); % [y, x] = [256, 22400]
[nsamples, nchirps, nframe] = size(data_3d); % [256, 128, 175]
Range_Index = (0:nr-1) / nr * Rmax; % 0 ~ Rmax 까지 Sample 수만큼 설정


%% Range FFT 후 시간축으로 쌓음 (1D FFT Spectrogram Average Threshold + α*sigma 적용)

data_2d_fft    = fft(data_2d);
data_2d_abs    = abs(data_2d_fft);  

Time_Index_All = (0:x-1) * PRI;          % 각 chirp 발생 시간 [s]
Range_Index    = (0:nr-1)/nr * Rmax;     % 거리축 [m]

% assume data_2d_abs, Time_Index_All, Range_Index already computed

% α = 0  (mean-only threshold)
alpha0       = 0;
thr0         = mean(data_2d_abs,1) + alpha0*std(data_2d_abs,[],1);
data_dyn0    = data_2d_abs - thr0;
data_dyn0(data_dyn0<0) = 0;
data_dyn0_dB = 20*log10(data_dyn0 + eps);

% α = 2  (mean+2σ threshold)
alpha2       = 2;
thr2         = mean(data_2d_abs,1) + alpha2*std(data_2d_abs,[],1);
data_dyn2    = data_2d_abs - thr2;
data_dyn2(data_dyn2<0) = 0;
data_dyn2_dB = 20*log10(data_dyn2 + eps);

% plot side by side
figure;
imagesc(Time_Index_All, Range_Index, data_dyn0_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average Threshold only (α=0)');

figure;
imagesc(Time_Index_All, Range_Index, data_dyn2_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average +2\sigma Threshold (α=2)');

%% CFAR 동적 Th (한 처프내에서)
windowLen = 32;           % 윈도우 크기
alpha     = 2;          % thr = alpha * localMean

% Range FFT magnitude
magR = abs(fft(data_2d, [], 1));  % [nr × nchirps]

% CFAR 처리 (movmean 사용)
% movmean(vec,windowLen,[],1) 은 열방향으로 윈도우 평균을 계산해 줍니다.
localMean = movmean(magR, windowLen, 1);   
thr       = alpha * localMean;            

% thresholding
cfar_map  = magR .* ( magR >= thr );      Z
% dB 변환 및 시각화
cfar_dB = 20*log10(cfar_map + eps);
figure;
imagesc(Time_Index_All, Range_Index, cfar_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('1D FFT Spectrogram with Sliding‐Window CFAR');
% CFAR Mean Threshold vs. Single Chirp FFT Magnitude
% 1) Compute 1D FFT magnitude and dynamic threshold mean
magR      = abs(fft(data_2d,[],1));           % [nr × nchirps]
localMean = movmean(magR, windowLen, 1);      % sliding-window mean
thr       = alpha * localMean;                % dynamic threshold
thr_mean  = mean(thr, 2);                     % time-averaged threshold [nr × 1]
% 2) Select a single chirp column for comparison
colIdx  = round(size(magR,2)/2);              % middle chirp index
sig_col = magR(:, colIdx);                    % [nr × 1]
% 3) Plot both on one figure
figure;
plot(Range_Index, sig_col,  'b-', 'LineWidth', 1.5);
hold on;
plot(Range_Index, thr_mean, 'r--','LineWidth', 1.5);
xlabel('Range (m)');
ylabel('Amplitude');
legend(sprintf('FFT Mag @ chirp %d', colIdx), 'Time-Averaged Threshold');
title('Signal vs. CFAR Mean Threshold'); grid on;

%% CFAR‐Aided Dynamic Thresholding over Frames → Range–Time Map
% parameters
windowLen = 32;
alpha     = 2;

[nr, nchirps, nframes] = size(data_3d);
cluster_dyn = zeros(nr, nchirps, nframes);

% 1) per-frame CFAR
for f = 1:nframes
    X    = data_3d(:,:,f);
    mag  = abs( fft(X, [], 1) );              % [nr×nchirps]
    th   = alpha * movmean(mag, windowLen, 1);
    mask = mag >= th;
    cluster_dyn(:,:,f) = mag .* mask;
end

% 2) flatten to full range–time
RT = reshape(cluster_dyn, nr, []);          % [nr × (nchirps*nframes)]
t  = (0:size(RT,2)-1) * PRI;                % time axis
r  = (0:nr-1)/nr * Rmax;                    % range axis

% 3) plot
figure;
imagesc(t, r, 20*log10(RT + eps));
axis xy tight;
colormap(parula);
xlabel('Time (s)');
ylabel('Range (m)');
title('Dynamic‐CFAR Range–Time Map');


%% Max Clustering per Chirp + Plotting

% 1) Max clustering
mag_all       = abs( fft(data_2d, [], 1) );   % [nr × nChirpsTot]
[nr, nChirps] = size(mag_all);
cluster_max   = zeros(nr, nChirps);

for k = 1:nChirps
    vec = mag_all(:,k);
    [~, idx] = max(vec);             % peak index
    keep     = idx;         % ±1 neighbors
    keep     = keep(keep>=1 & keep<=nr);
    cluster_max(keep, k) = vec(keep);
end

% 2) Plot Range–Time Map of clustered peaks
time_all = (0:nChirps-1) * PRI;      % time axis [s]
figure;
imagesc(time_all, Range_Index, 20*log10(cluster_max + eps));
axis xy; axis tight;
colormap(parula);
xlabel('Time (s)');
ylabel('Range (m)');
title('Max Clustered Range–Time Map');



