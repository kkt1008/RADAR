%% ========== 초기 설정 및 환경 정리 ==========
clc; clear; close all;  % 콘솔 출력 초기화, 모든 변수 삭제, 열린 Figure 창 닫기

%% ========== 레이더 시스템 및 파형 파라미터 설정 ==========
% 수신 안테나 수
nrx = 4;

% FMCW 파라미터: 레이더가 송신하는 chirp 신호의 설정
startFreq        = 77;        % 시작 주파수 [GHz]
idleTime         = 100;       % chirp 후 IDLE 시간 [us]
adcStartTime     = 6;         % chirp 시작 후 ADC가 샘플링 시작하는 시간 [us]
rampEndTime      = 60;        % chirp 하나의 전체 시간 [us]
digOutSampleRate = 4270;      % ADC 샘플링 주파수 [kHz]
freqSlopeConst   = 66.143;    % chirp의 주파수 증가 속도 [MHz/us]
numADCSamples    = 256;       % 1 chirp 동안 수집하는 샘플 수 → Range FFT 크기

% 프레임 구성: chirp를 묶은 단위
chirpStartIdx    = 0;
chirpEndIdx      = 1;         % 이 설정은 송신 안테나 수 추정에도 사용됨
numLoops         = 128;       % 1프레임당 chirp 수 → Doppler FFT 크기
framePeriodicity = 40.0;      % 프레임 간격 시간 [ms]

% ADC 관련 파라미터
numADCBits = 16;              % ADC의 비트 해상도
isReal     = 1;               % 실수값 수집 여부 (1이면 실수값)

% === 파생 파라미터 계산 (상수와 설정 기반) ===
c0             = 299792458;                  % 빛의 속도 [m/s]
nrx            = 1;                          % 실제 처리에서는 1채널만 사용
ntx            = chirpEndIdx - chirpStartIdx + 1; % 송신 안테나 수 (여기선 2)
nd             = numLoops;                   % Doppler FFT 길이
nr             = numADCSamples;              % Range FFT 길이
tfr            = framePeriodicity * 1e-3;    % 프레임 간격 [s]
fs             = digOutSampleRate * 1e3;     % 샘플링 주파수 [Hz]
kf             = freqSlopeConst * 1e12;      % 주파수 증가율 [Hz/s]
adc_duration   = nr / fs;                    % ADC 샘플링 시간 [s]
BW             = adc_duration * kf;          % chirp 대역폭 = 샘플 시간 × 슬로프
PRI            = (idleTime + rampEndTime) * 1e-6; % Pulse Repetition Interval [s]
PRF            = 1 / PRI;                    % Pulse Repetition Frequency [Hz]
fc             = startFreq * 1e9 + kf * (adcStartTime * 1e-6 + adc_duration / 2); % 중심 주파수 추정
frame_duration = PRI * numLoops;             % 프레임 하나에 걸리는 총 시간
frame_idletime = tfr - frame_duration;       % 프레임 간 Idle 시간
Rmax           = fs * c0 / (2 * kf);         % 최대 거리 = λ/2 계산 기반

%% ========== Raw 데이터 불러오기 (DCA1000 형식) ==========
% adc_file_name = 'hand_range.bin';            % 수집된 bin 파일 이름
adc_file_name = 'human_azimuth.bin';
data = readDCA1000(adc_file_name);        % raw binary 파일에서 데이터 읽기
disp('data was generated!');

% 수신 채널 2의 전체 데이터 추출 (1×5734400 벡터 등)
data_1d = data(2, :);

% 2D 형태로 변환: (샘플 수 × chirp 수) = (256 × 전체 chirp 수)
data_2d = reshape(data_1d, [nr, length(data_1d)/nr]);
[x, y] = size(data_2d);  % x=256 (range bin 수), y=chirp 개수

% 거리 및 시간축 인덱스 생성
Range_Index     = (0:nr-1) / nr * Rmax;       % 거리축
Time_Index_All  = (0:y-1) * PRI;              % 시간축

% CFAR 설정값 정의 
numthreshold = 32;           % 전체 FFT 결과를 몇 블록으로 나눌 것인가
weight       = 4;            % CFAR threshold 가중치
hwin         = hamming(nr);  % Hamming 윈도우 생성
CFAR_spectrogram = zeros(nr, y); % 최종 결과 저장용 행렬

% 모든 chirp에 Hamming, FFT, CFAR 적용 
for k = 1:y
    % 1. 현재 chirp에 Hamming 윈도우 적용
    sample = data_2d(:,k) .* hwin;

    % 2. FFT 수행 후 절대값으로 magnitude 계산
    fft_out = abs(fft(sample));

    % 3. CFAR 처리: 전체 벡터를 블록(열) 단위로 나눔
    reshaped = reshape(fft_out, nr/numthreshold, numthreshold);

    % 4. 각 블록의 평균을 기반으로 threshold 계산 후 가중치 적용
    threshold = mean(reshaped) * weight;

    % 5. threshold와 비교하여 해당 블록의 값이 작으면 제거
    for i = 1:numthreshold
        for j = 1:nr/numthreshold
            if reshaped(j,i) <= threshold(i)
                reshaped(j,i) = 0;
            end
        end
    end
    % 6. 다시 1D로 reshape하여 최종 결과 저장
    filtered = reshape(reshaped, nr, 1);
    CFAR_spectrogram(:,k) = filtered;
end

% 마지막 chirp의 threshold를 256포인트로 확장
threshold_full = repelem(threshold, nr/numthreshold);

% FFT 결과와 threshold를 dB 스케일 없이 시각화
figure;
plot(Range_Index, fft_out, 'b'); hold on;
plot(Range_Index, threshold_full, 'r--');
xlabel('Range (m)');
ylabel('Amplitude');
legend('FFT Magnitude', 'CFAR Threshold');
title('Last Chirp FFT and CFAR Threshold');

% 전체 chirp Spectrogram 출력 (Range-Time Map) 
CFAR_spectrogram_dB = 20 * log10(CFAR_spectrogram + eps); % dB 스케일 변환
figure;
imagesc(Time_Index_All, Range_Index, CFAR_spectrogram_dB);
axis xy; axis tight;
xlabel('Time (s)');
ylabel('Range (m)');
title(sprintf('Integrated Spectrogram (CFAR + Hamming), weight = %.1f', weight));
colorbar;
