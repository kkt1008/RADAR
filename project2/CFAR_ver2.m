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
adc_file_name = ['hand_range.bin'];
%adc_file_name = 'human_azimuth.bin';
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

% 3) dynamic threshold 계산 (열별 평균+α·표준편차)
alpha = 2;                                               % 문턱 곱셈 계수
thr   = mean(data_2d_abs,1) + alpha*std(data_2d_abs,[],1); % 1×x

% 4) 문턱 적용 — 음수는 0으로 클리핑
data_dyn = data_2d_abs - thr;       % thresholding
data_dyn(data_dyn < 0) = 0;

% dc_bins = find(Range_Index < 0.15);  % 0.1 m 이하 필터링
% data_dyn(dc_bins, :) = 0;

% 0.15 ~ 1 m 구간 외 모든 값 제거 (원타겟 이외 클러스트 제거)
% mask_out = Range_Index < 0.15 | Range_Index > 1.0;
% data_dyn(mask_out, :) = 0;

data_dyn_dB = 20*log10(data_dyn + eps);

figure;
imagesc(Time_Index_All, Range_Index, data_dyn_dB);
axis xy;
axis tight;
xlabel('Time(s)');
ylabel('Range(m)');
title('1D FFT Spectrogram (Dynamic Threshold + DC Filter)');

%% CFAR 파라미터
windowLen = 32;           % 윈도우 크기
alpha     = 2;          % thr = alpha * localMean

% Range FFT magnitude
magR = abs(fft(data_2d, [], 1));  % [nr × nchirps]

% CFAR 처리 (movmean 사용)
% movmean(vec,windowLen,[],1) 은 열방향으로 윈도우 평균을 계산해 줍니다.
localMean = movmean(magR, windowLen, 1);   
thr       = alpha * localMean;            

% thresholding
cfar_map  = magR .* ( magR >= thr );      

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
title('Signal vs. CFAR Mean Threshold');
grid on;

%% Hough Transform 기반 궤적 클러스터링 (추가)
% - 1D FFT Spectrogram (data_dyn_dB)에서 강한 궤적을 직선으로 검출
% - Hough Transform → houghlines로 궤적 픽셀만 cluster_hough에 남김

% 1) Binary mask 생성 (예: –40 dB 이상만 유효 신호)
BW = data_dyn_dB > -40;  % [nr × nChirpsTot]

% 2) Hough Transform
[H, theta, rho] = hough(BW);

% 3) Peak 검출
P = houghpeaks(H, 5, 'Threshold', ceil(0.3*max(H(:))));

% 4) 선분 정보 추출
lines = houghlines(BW, theta, rho, P, 'FillGap', 5, 'MinLength', 20);

% 5) 궤적만 남길 빈 맵 준비
cluster_hough = zeros(size(BW));

% 6) 각 선분을 따라 픽셀들만 보존
for k = 1:length(lines)
    p1 = lines(k).point1;  % [x1 y1]
    p2 = lines(k).point2;  % [x2 y2]
    % 선분 위의 모든 좌표 계산
    numPts = max(abs(p2 - p1)) + 1;
    x = round(linspace(p1(1), p2(1), numPts));
    y = round(linspace(p1(2), p2(2), numPts));
    idx = sub2ind(size(BW), y, x);
    cluster_hough(idx) = magR_all(idx);  % 원본 magnitude
end

% 7) 결과 플롯
figure;
imagesc(time_all, range_vec, 20*log10(cluster_hough + eps));
axis xy; axis tight;
colormap(jet);
xlabel('Time (s)');
ylabel('Range (m)');
title('Range–Time Map via Hough Transform Clustering');
colorbar;
