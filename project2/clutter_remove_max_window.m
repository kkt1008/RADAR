%% Data Set 불러오기 & Radar Information
clc;
clear;
% close all;

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
Time_Index_All = (0:x-1) * PRI;     % 시간축: 각 chirp 발생 시간 [s]

%% Range FFT 후 시간축으로 쌓음 (1D FFT Spectrogram Average Threshold 적용)
data_2d_1D_FFT = fft(data_2d,[],1); % data_2d에 대해 열방향으로 fft진행(for Range)
data_2d_1D_FFT = data_2d_1D_FFT/length(data_2d_1D_FFT); % FFT Normalize
data_2d_1D_FFT_abs = abs(data_2d_1D_FFT); % 절대값 취하기
data_2d_1D_FFT_abs_dB = 20*log10(data_2d_1D_FFT_abs + eps); % dB scaling

data_2d_1D_FFT_abs_Normalized = data_2d_1D_FFT_abs - mean(data_2d_1D_FFT_abs, 1); % 열방향으로 평균 빼기
data_2d_1D_FFT_abs_Normalized_temp = data_2d_1D_FFT_abs_Normalized;
data_2d_1D_FFT_abs_Normalized_temp(data_2d_1D_FFT_abs_Normalized < 0) = 0; % 평균 이하 값 0으로 처리
data_2d_1D_FFT_abs_Normalized_Filtered = data_2d_1D_FFT_abs_Normalized_temp;
data_2d_1D_FFT_abs_Normalized_Filtered_dB = real(20*log10(data_2d_1D_FFT_abs_Normalized_Filtered + eps)); % dB scaling

% Range 별 Amplitude 그래프 

%% 윈도우 기반 MTI 및 최대값 기반 필터링
window_size = 10; % 윈도우 크기 (chirp 수)
data_mti_processed = zeros(size(data_2d_1D_FFT_abs)); % 결과 저장 (256, 22400)

for col = 1:x % 각 열(256x1)에 대해 처리
    % 1. 슬라이딩 윈도우 평균 제거 (MTI)
    current_col = data_2d_1D_FFT_abs(:, col);
    window_mean = movmean(data_2d_1D_FFT_abs(:, col), window_size, 'Endpoints', 'fill');
    mti_data = current_col - window_mean;
    
    % 2. 음수 값을 0으로 처리
    mti_data(mti_data < 0) = 0;
    
    % 3. 최대값 인덱스 찾기
    [~, max_idx] = max(mti_data);
    
    % 4. 최대값 중심으로 좌우 5개씩 (총 11개) 선택
    new_col = zeros(nr, 1); % 256x1 벡터 초기화
    start_idx = max(1, max_idx - 5); % 좌측 경계
    end_idx = min(nr, max_idx + 5); % 우측 경계
    new_col(start_idx:end_idx) = mti_data(start_idx:end_idx);
    
    % 5. 결과 저장
    data_mti_processed(:, col) = new_col;
end

% dB 스케일링
data_mti_processed_dB = 20*log10(data_mti_processed + eps);

% 3D 데이터로 재구성 (옵션)
data_3d_mti = reshape(data_mti_processed_dB, [nr, 128, 175]);

%% 결과 시각화
figure;
imagesc(Time_Index_All, Range_Index, data_2d_1D_FFT_abs_dB);
title('클러터 제거 전 (dB)');
axis xy; % Y축을 거리(상향)로 설정
axis tight; % 축 범위 조정
xlabel('Time (s)'); ylabel('Range (m)');
colorbar;

figure;
imagesc(Time_Index_All, Range_Index, data_mti_processed_dB);
title('After MTI and Max Filtering (dB)');
axis xy; % Y축을 거리(상향)로 설정
axis tight; % 축 범위 조정
xlabel('Time (s)'); ylabel('Range (m)');
colorbar;

