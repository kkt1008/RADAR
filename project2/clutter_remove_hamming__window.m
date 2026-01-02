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
%% sidelobe 비교 (Hamming 윈도우 전vs후)
%% ===== Range FFT: Rect(무윈도우) vs Hamming(윈도우) =====
wR = hamming(nr, 'periodic');      % nr=256
cg = sum(wR)/nr;                  % coherent gain (메인피크 크기 보정용)

% (nr x nChirps) 에 윈도우 곱하기 (implicit expansion)
data_2d_ham = data_2d .* wR;

% Range FFT
FFT_rect = fft(data_2d,     [], 1) / nr;
FFT_ham  = fft(data_2d_ham, [], 1) / (nr * cg);  % cg로 보정(선택이지만 비교에 유리)

mag_rect = abs(FFT_rect);
mag_ham  = abs(FFT_ham);

%% ===== PSL(peak sidelobe level) 계산 =====
guard = 6;  % 메인로브 제외 구간(빈). Hamming은 메인로브가 넓어서 4~8 사이로 조절 추천.

psl_rect_db = zeros(1, x);
psl_ham_db  = zeros(1, x);

for col = 1:x
    % --- Rect ---
    [main_amp, main_idx] = max(mag_rect(:, col));
    tmp = mag_rect(:, col);
    lo = max(1, main_idx - guard);
    hi = min(nr, main_idx + guard);
    tmp(lo:hi) = 0;  % 메인로브(가드 구간) 제거
    side_amp = max(tmp);
    psl_rect_db(col) = 20*log10((side_amp + eps) / (main_amp + eps));

    % --- Hamming ---
    [main_amp2, main_idx2] = max(mag_ham(:, col));
    tmp2 = mag_ham(:, col);
    lo2 = max(1, main_idx2 - guard);
    hi2 = min(nr, main_idx2 + guard);
    tmp2(lo2:hi2) = 0;
    side_amp2 = max(tmp2);
    psl_ham_db(col) = 20*log10((side_amp2 + eps) / (main_amp2 + eps));
end

% 요약 출력 (PSL은 음수 dB가 정상: 더 작은 값일수록 sidelobe가 더 억제됨)
fprintf('\n[PSL 비교] guard=%d bins\n', guard);
fprintf('Rect  PSL: median = %.2f dB, mean = %.2f dB, worst = %.2f dB\n', ...
    median(psl_rect_db), mean(psl_rect_db), max(psl_rect_db));
fprintf('Hamm  PSL: median = %.2f dB, mean = %.2f dB, worst = %.2f dB\n', ...
    median(psl_ham_db), mean(psl_ham_db), max(psl_ham_db));

impr_med = median(psl_rect_db) - median(psl_ham_db);  % (예: -13 - (-42) = +29 dB 개선)
fprintf('Median PSL improvement: %.2f dB (≈ x%.1f amplitude reduction)\n', ...
    impr_med, 10^(impr_med/20));

%% ===== 대표 chirp 하나 골라서 before/after 레인지 컷 플롯 =====
[~, col0] = max(max(mag_rect, [], 1));  % 전체 중 피크가 가장 큰 chirp 선택(대표)
prof_rect_db = 20*log10(mag_rect(:, col0) + eps);
prof_ham_db  = 20*log10(mag_ham(:,  col0) + eps);

figure;
plot(Range_Index, prof_rect_db); hold on;
plot(Range_Index, prof_ham_db);
grid on; axis tight;
xlabel('Range (m)'); ylabel('Magnitude (dB)');
title(sprintf('Range profile (chirp #%d): Rect vs Hamming', col0));
legend('Rect (no window)', 'Hamming window');


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

% 3. 최대값 인덱스 찾기 및 좌우 5개 샘플 선택
[~, max_indices] = max(data_mti_processed, [], 1); % 모든 열의 최대값 인덱스
data_mti_processed_new = zeros(size(data_mti_processed));
for col = 1:x
    start_idx = max(1, max_indices(col));
    end_idx = min(nr, max_indices(col));
    data_mti_processed_new(start_idx:end_idx, col) = data_mti_processed(start_idx:end_idx, col);
end

% 4. dB 스케일링
data_mti_processed_dB = 20*log10(data_mti_processed_new + eps);

% 5. 3D 데이터로 재구성 (옵션)
data_3d_mti = reshape(data_mti_processed_dB, [nr, 128, 175]);

%% 결과 시각화
figure;
imagesc(Time_Index_All, Range_Index, 20*log10(data_2d_1D_FFT_abs + eps));
title('Before apply window (dB)');
axis xy; % Y축을 거리(상향)로 설정
axis tight; % 축 범위 조정
xlabel('Time (s)'); ylabel('Range (m)');
colorbar;

figure;
imagesc(Time_Index_All, Range_Index, data_mti_processed_dB);
title('After apply hamming window (dB)');
axis xy; % Y축을 거리(상향)로 설정
axis tight; % 축 범위 조정
xlabel('Time (s)'); ylabel('Range (m)');
colorbar;
