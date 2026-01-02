%% Data Set 불러오기 & Radar Information
clc;
clear;
%close all;

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

data_2d_fft    = fft(data_2d);  % 각 열마다 FFT수행하여 주파수 영역으로 변환
data_2d_abs    = abs(data_2d_fft);  % 절대값을 통해 magnitude계산(복수수 신호 제거효과)

Time_Index_All = (0:x-1) * PRI;          % 각 chirp 발생 시간 [s]
Range_Index    = (0:nr-1)/nr * Rmax;     % 거리축 [m]

% assume data_2d_abs, Time_Index_All, Range_Index already computed

%% α = 0  (mean-only threshold) 평균만 사용하는 효과
alpha0       = 0;
thr0         = mean(data_2d_abs,1) + alpha0*std(data_2d_abs,[],1); %각 열(chirp)의 평균 계산
data_dyn0    = data_2d_abs - thr0; % magnituded에서 평균을 빼서 clutter제거 
data_dyn0(data_dyn0<0) = 0; % 음수일 경우 0으로 할당
data_dyn0_dB = 20*log10(data_dyn0 + eps); % dB단위로 변환

figure;
imagesc(Time_Index_All, Range_Index, data_dyn0_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average Threshold only (α=0)-평균값 기반 제거와 동일');

%% α = 1  (mean+σ threshold)
alpha1       = 1;
thr2         = mean(data_2d_abs,1) + alpha1*std(data_2d_abs,[],1);
data_dyn2    = data_2d_abs - thr2;
data_dyn2(data_dyn2<0) = 0;
data_dyn2_dB = 20*log10(data_dyn2 + eps);

figure;
imagesc(Time_Index_All, Range_Index, data_dyn2_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average +1\sigma Threshold (α=1)');

%% α = 2  (mean+2σ threshold)
alpha2       = 2;
thr2         = mean(data_2d_abs,1) + alpha2*std(data_2d_abs,[],1);
data_dyn2    = data_2d_abs - thr2;
data_dyn2(data_dyn2<0) = 0;
data_dyn2_dB = 20*log10(data_dyn2 + eps);

figure;
imagesc(Time_Index_All, Range_Index, data_dyn2_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average +2\sigma Threshold (α=2)');

%% α = 3  (mean+3σ threshold)
alpha3       = 3;
thr2         = mean(data_2d_abs,1) + alpha3*std(data_2d_abs,[],1);
data_dyn2    = data_2d_abs - thr2;
data_dyn2(data_dyn2<0) = 0;
data_dyn2_dB = 20*log10(data_dyn2 + eps);

figure;
imagesc(Time_Index_All, Range_Index, data_dyn2_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average +3\sigma Threshold (α=3)');

%% α = 4  (mean+4σ threshold)
alpha4       = 4;
thr2         = mean(data_2d_abs,1) + alpha4*std(data_2d_abs,[],1);
data_dyn2    = data_2d_abs - thr2;
data_dyn2(data_dyn2<0) = 0;
data_dyn2_dB = 20*log10(data_dyn2 + eps);

figure;
imagesc(Time_Index_All, Range_Index, data_dyn2_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average +4\sigma Threshold (α=4)');

%% α = 5  (mean+5σ threshold)
alpha5       = 5;
thr2         = mean(data_2d_abs,1) + alpha5*std(data_2d_abs,[],1);
data_dyn2    = data_2d_abs - thr2;
data_dyn2(data_dyn2<0) = 0;
data_dyn2_dB = 20*log10(data_dyn2 + eps);

figure;
imagesc(Time_Index_All, Range_Index, data_dyn2_dB);
axis xy; axis tight;
xlabel('Time (s)'); ylabel('Range (m)');
title('Average +5\sigma Threshold (α=5)');

% Range FFT 후 시간축으로 쌓음 (1D FFT Spectrogram Average Threshold 적용)
data_2d_1D_FFT = fft(data_2d, [], 1); % data_2d에 대해 열방향으로 fft진행(for Range)
data_2d_1D_FFT = data_2d_1D_FFT / nr; % FFT 정규화
data_2d_1D_FFT_abs = abs(data_2d_1D_FFT); % 절대값 취하기
data_2d_1D_FFT_abs_dB = 20*log10(data_2d_1D_FFT_abs + eps); % dB scaling
