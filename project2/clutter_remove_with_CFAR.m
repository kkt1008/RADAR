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
adc_file_name = 'human_azimuth.bin';
data = readDCA1000(adc_file_name);
disp('data was generated!');

% channel 1's All Data Load : (1, 5734400) 행 vector
data_1d = data(2,1:end);  
% channel 1's First Frame First Chirp vector (1, 256)
data_sample = data(2,1:256); 

% data_1d을 (256, length(data_1d)/256)으로 나눠서 2차원 행렬로 바꿈 (256,22400)
data_2d = reshape(data_1d, [nr, length(data_1d)/nr]);

% (256, 22400) -> (256, 128, 175)로 변환
data_3d = reshape(data_2d, [nr, 128, 175]); 

[i, j] = size(data_1d); % [1, 5734400]
[y, x] = size(data_2d); % [256, 22400]
[nsamples, nchirps, nframe] = size(data_3d); % [256, 128, 175]
Range_Index = (0:nr-1) / nr * Rmax; % 거리축: 0 ~ Rmax 까지 Sample 수만큼 설정
Time_Index_All = (0:x-1) * PRI;     % 시간축: 각 chirp 발생 시간 [s]

%% Range FFT
data_2d_1D_FFT = fft(data_2d, [], 1); % Range FFT (열방향)
data_2d_1D_FFT = data_2d_1D_FFT / length(data_2d_1D_FFT); % FFT Normalize
data_2d_1D_FFT_abs = abs(data_2d_1D_FFT); % 절대값

%% 윈도우 기반 MTI 및 최대값 기반 필터링 (최적화)
window_size = 10; % 윈도우 크기 (chirp 수)
window = hamming(window_size)'; % Hamming 윈도우
window = window / sum(window); % 정규화
data_mti_processed = zeros(size(data_2d_1D_FFT_abs)); % 결과 저장 (256, 22400)

% 1. 윈도우 기반 MTI (Hamming 윈도우)
for col = 1:x
    window_mean = conv(data_2d_1D_FFT_abs(:, col), window, 'same');
    data_mti_processed(:, col) = data_2d_1D_FFT_abs(:, col) - window_mean;
end

% 2. 음수 값을 0으로 처리
data_mti_processed(data_mti_processed < 0) = 0;
data_mti_processed_dB = 20*log10(data_mti_processed_new + eps);

figure;
imagesc(Range_Index, 1:x, data_mti_processed_dB.');
title('apply hamming window (dB)');
xlabel('Range (m)'); ylabel('Chirp Index');
colorbar;

% 3. 최대값 인덱스 찾기 및 좌우 5개 샘플 선택
[~, max_indices] = max(data_mti_processed, [], 1); % 모든 열의 최대값 인덱스
data_mti_processed_new = zeros(size(data_mti_processed));
for col = 1:x
    start_idx = max(1, max_indices(col) - 5);
    end_idx = min(nr, max_indices(col) + 5);
    data_mti_processed_new(start_idx:end_idx, col) = data_mti_processed(start_idx:end_idx, col);
end

% 4. dB 스케일링
data_mti_processed_dB = 20*log10(data_mti_processed_new + eps);

%% 수동 CA-CFAR 적용
% CFAR 파라미터 설정
num_training_cells = 20; % 훈련 셀 수
num_guard_cells = 4; % 보호 셀 수
alpha = 10; % 스케일링 팩터 (조정 가능)

% CFAR 적용 (Range 축 방향)
cfar_output = zeros(size(data_mti_processed_dB));
for col = 1:x
    cfar_output(:, col) = ca_cfar_1d(data_mti_processed_dB(:, col), num_training_cells, num_guard_cells, alpha);
end

% CFAR 결과 적용: 타겟만 남기기
data_cfar_processed_dB = data_mti_processed_dB;
data_cfar_processed_dB(~cfar_output) = -300; % 타겟이 아닌 부분은 -300 dB로 설정

% 5. 3D 데이터로 재구성 (옵션)
data_3d_mti = reshape(data_mti_processed_dB, [nr, 128, 175]);
data_3d_cfar = reshape(data_cfar_processed_dB, [nr, 128, 175]);

%% 결과 시각화
figure;
% Before MTI
subplot(1, 3, 1);
imagesc(Range_Index, 1:x, 20*log10(data_2d_1D_FFT_abs + eps).');
title('Before MTI (dB)');
xlabel('Range (m)'); ylabel('Chirp Index');
colorbar;

% After MTI and Max Filtering
subplot(1, 3, 2);
imagesc(Range_Index, 1:x, data_mti_processed_dB.');
title('After MTI and Max Filtering (dB)');
xlabel('Range (m)'); ylabel('Chirp Index');
colorbar;

% After CFAR
subplot(1, 3, 3);
imagesc(Range_Index, 1:x, data_cfar_processed_dB.');
title('After CFAR (dB)');
xlabel('Range (m)'); ylabel('Chirp Index');
colorbar;

%% CFAR 함수 직접 구현
function cfar_output = ca_cfar_1d(data, num_training_cells, num_guard_cells, alpha)
    % data: 입력 데이터 (1D 벡터, 예: 256x1)
    % num_training_cells: 훈련 셀 수 (양쪽 합)
    % num_guard_cells: 보호 셀 수 (양쪽 합)
    % alpha: 스케일링 팩터 (임계값 조정)
    % cfar_output: 이진 출력 (1: 타겟, 0: 비타겟)

    N = length(data); % 데이터 길이 (예: 256)
    half_training = floor(num_training_cells / 2); % 한쪽 훈련 셀 수
    half_guard = floor(num_guard_cells / 2); % 한쪽 보호 셀 수
    cfar_output = zeros(size(data)); % 출력 초기화

    for i = 1:N
        % 훈련 셀 범위 계산
        start_idx = max(1, i - half_training - half_guard);
        end_idx = min(N, i + half_training + half_guard);
        
        % 보호 셀 제외한 훈련 셀 추출
        training_region = data(start_idx:end_idx);
        guard_start = max(1, i - half_guard);
        guard_end = min(N, i + half_guard);
        
        % 훈련 셀에서 보호 셀 제외
        if guard_start > start_idx
            training_region(1:(guard_start - start_idx)) = [];
        end
        if guard_end < end_idx
            training_region((guard_end - start_idx + 1):end) = [];
        end
        
        % 훈련 셀 평균 계산
        if ~isempty(training_region)
            noise_level = mean(training_region);
        else
            noise_level = 0; % 경계 처리
        end
        
        % 임계값 계산
        threshold = noise_level + alpha;
        
        % 타겟 검출
        if data(i) > threshold
            cfar_output(i) = 1;
        end
    end
end