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
isReal = 1; % set to 1 if real only data, 0 if complex data

% Radar Parameters
c0              =   299792458; % 빛의 속도 
nrx             =   1;
ntx             =   chirpEndIdx - chirpStartIdx + 1; % 1 - 0 + 1 = 2
nd              =   numLoops; % Number of Chirps per Frame : N_DopplerFFT = 128
nr              =   numADCSamples; % Number of Samples per Chirp : N_RangeFFT = 256
tfr             =   framePeriodicity * 1e-3; % Frame 주기[s]
fs              =   digOutSampleRate * 1e3; % Sampling Frequency[Hz]
kf              =   freqSlopeConst * 1e12; % 주파수 Sweep 속도 S [Hz/s]
adc_duration    =   nr / fs; % Chirp Duration Tc : Sweep 시간 per 1chirp
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
data_1d=data(2,1:end);  
% channel 1's Fist Frame First Chirp vector (1, 256)
data_sample = data(2,1:256); 

% data1을 (256, length(data1)/256)으로 나눠서 2차원 행렬로 바꿈 (256,22400)
data_2d=reshape(data_1d,[nr,length(data_1d)/nr]);

% (256, 22400) -> (256, 128, 175)로 바꿔야 함
% (256, 128, 175) = (256 Samples per 1Chirp, 128 Chirps, 175 Frames)
data_3d = reshape(data_2d, [nr, 128, 175]); 

[i, j]=size(data_1d); % [i, j] = [1, 5734400]
[y, x]=size(data_2d); % [y, x] = [256, 22400]
[nsamples, nchirps, nframe] = size(data_3d); % [256, 128, 175]
Range_Index = (0:nr-1) / nr * Rmax; % 0 ~ Rmax 까지 Sample 수만큼 설정

%% 1 Chirp Samples Sample Plot (Average Threshold 적용)
first_chirp_range = fft(data_sample); % FFT
first_chirp_range_abs = abs(first_chirp_range/length(first_chirp_range)); % FFT Normalization -> 절대값
first_chirp_range_abs_dB = 20*log10(first_chirp_range_abs + eps); % dB scaling, eps to avoid log(0)

first_chirp_range_abs_Normalized = first_chirp_range_abs - mean(first_chirp_range_abs); % 정규화(평균 뺌)
first_chirp_range_abs_Normalized_temp = first_chirp_range_abs_Normalized;
first_chirp_range_abs_Normalized_temp(first_chirp_range_abs_Normalized_temp < 0) = 0; % 평균 이하 값 0으로 처리
first_chirp_range_abs_Normalized_Filtered = first_chirp_range_abs_Normalized_temp;
first_chirp_range_abs_Normalized_Filtered_dB = 20*log10(first_chirp_range_abs_Normalized_Filtered + eps); % dB scaling

figure;
subplot(4, 2, 1);
plot(Range_Index, first_chirp_range_abs);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude')
title('1Chirp Samples Range Plot (window x / dB x / Avg Filter x)')

subplot(4, 2, 3);
plot(Range_Index, first_chirp_range_abs_dB);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude (dB)')
title('1Chirp Samples Range Plot (window x / dB o / Avg Filter x)')

subplot(4, 2, 5);
plot(Range_Index, first_chirp_range_abs_Normalized_Filtered);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude')
title('1Chirp Samples Range Plot (window x / dB x / Avg Filter o)')

subplot(4, 2, 7);
plot(Range_Index, first_chirp_range_abs_Normalized_Filtered_dB);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude (dB)')
title('1Chirp Samples Range Plot (window x / dB o / Avg Filter o)')

data_sample_windowed = data_sample .* hamming(256)'; % Hamming Window 적용
first_chirp_range = fft(data_sample_windowed); % FFT
first_chirp_range_abs = abs(first_chirp_range/length(first_chirp_range)); % FFT Normalization -> 절대값
first_chirp_range_abs_dB = 20*log10(first_chirp_range_abs + eps); % dB scaling

first_chirp_range_abs_Normalized = first_chirp_range_abs - mean(first_chirp_range_abs); % 정규화(평균 뺌)
first_chirp_range_abs_Normalized_temp = first_chirp_range_abs_Normalized;
first_chirp_range_abs_Normalized_temp(first_chirp_range_abs_Normalized_temp < 0) = 0; % 평균 이하 값 0으로 처리
first_chirp_range_abs_Normalized_Filtered = first_chirp_range_abs_Normalized_temp;
first_chirp_range_abs_Normalized_Filtered_dB = 20*log10(first_chirp_range_abs_Normalized_Filtered + eps); % dB scaling

subplot(4, 2, 2);
plot(Range_Index, first_chirp_range_abs);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude')
title('1Chirp Samples Range Plot (window o / dB x / Avg Filter x)')

subplot(4, 2, 4);
plot(Range_Index, first_chirp_range_abs_dB);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude (dB)')
title('1Chirp Samples Range Plot (window o / dB o / Avg Filter x)')

subplot(4, 2, 6);
plot(Range_Index, first_chirp_range_abs_Normalized_Filtered);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude')
title('1Chirp Samples Range Plot (window o / dB x / Avg Filter o)')

subplot(4, 2, 8);
plot(Range_Index, first_chirp_range_abs_Normalized_Filtered_dB);
axis tight;
xlabel('Range(m)')
ylabel('Amplitude (dB)')
title('1Chirp Samples Range Plot (window o / dB o / Avg Filter o)')

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

figure;
subplot(4, 2, 1);
Time_Index_All_Chirps = (1:x)/x*tfr*nframe; 
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window x / dB x / Avg Filter x)')

subplot(4, 2, 3);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs_dB);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window x / dB o / Avg Filter x)')

subplot(4, 2, 5);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs_Normalized_Filtered);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window x / dB x / Avg Filter o)')

subplot(4, 2, 7);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs_Normalized_Filtered_dB);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window x / dB o / Avg Filter o)')

data_2d_windowed = data_2d .* repmat(hamming(256), 1, x); % FFT전 Hamming Window 적용
data_2d_1D_FFT=fft(data_2d_windowed,[],1); % data_2d에 대해 열방향으로 fft진행(for Range)
data_2d_1D_FFT = data_2d_1D_FFT/length(data_2d_1D_FFT); % FFT Normalize
data_2d_1D_FFT_abs = abs(data_2d_1D_FFT); % 절대값 취하기
data_2d_1D_FFT_abs_dB = 20*log10(data_2d_1D_FFT_abs + eps); % dB scaling

data_2d_1D_FFT_abs_Normalized = data_2d_1D_FFT_abs - mean(data_2d_1D_FFT_abs, 1); % 열방향으로 평균 빼기
data_2d_1D_FFT_abs_Normalized_temp = data_2d_1D_FFT_abs_Normalized;
data_2d_1D_FFT_abs_Normalized_temp(data_2d_1D_FFT_abs_Normalized < 0) = 0; % 평균 이하 값 0으로 처리
data_2d_1D_FFT_abs_Normalized_Filtered = data_2d_1D_FFT_abs_Normalized_temp;
data_2d_1D_FFT_abs_Normalized_Filtered_dB = real(20*log10(data_2d_1D_FFT_abs_Normalized_Filtered + eps)); % dB scaling

subplot(4, 2, 2);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window o / dB x / Avg Filter x)')

subplot(4, 2, 4);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs_dB);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window o / dB o / Avg Filter x)')

subplot(4, 2, 6);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs_Normalized_Filtered);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window o / dB x / Avg Filter o)')

subplot(4, 2, 8);
imagesc(Time_Index_All_Chirps, Range_Index, data_2d_1D_FFT_abs_Normalized_Filtered_dB);
axis xy; 
axis tight;
xlabel('Time(s)')
ylabel('Range(m)')
title('1D FFT Spectrogram (window o / dB o / Avg Filter o)')

%% RD-Map (클러터 제거 추가)
frame_number = 1;
rd = (0:nr-1) / nr * fs * c0/(2*kf); % Range 축
vd = (-nd/2:1:nd/2-1)/nd*PRF*c0/2/fc; % Doppler 축
RD_Map_3D = zeros(256, 128, 175); % RD-Map 저장 배열

% 클러터 제거를 위한 평균 계산 (모든 프레임에 대한 평균)
mean_frame_data = mean(data_3d, 3); % 프레임 차원(3번째 차원)에서 평균 계산

% High-Pass 필터 설계 (Doppler 차원에서 저주파 제거)
[b, a] = butter(2, 0.1, 'high'); % 2차 Butterworth High-Pass 필터 (정규화된 주파수 0.1)

figure;
while frame_number < 176
    % 1개 Frame Data(256, 128) 갖고 오기
    frame_data = data_3d(:,:,frame_number);
    
    % 평균 제거 (클러터 제거)
    frame_data_clutter_removed = frame_data - mean_frame_data;
    
    % 열방향 FFT 진행 (Range FFT)
    frame_data_1D_FFT = fft(frame_data_clutter_removed, [], 1);
    frame_data_1D_FFT_Normalized = frame_data_1D_FFT / length(frame_data_1D_FFT(:, 1)); % FFT Normalize
    
    % 행방향 FFT 진행 (Doppler FFT)
    frame_data_2D_FFT = fftshift(fft(frame_data_1D_FFT_Normalized, [], 2), 2);
    
    % MTI 필터링: Doppler 차원에서 High-Pass 필터 적용
    for r = 1:nr
        frame_data_2D_FFT(r, :) = filtfilt(b, a, frame_data_2D_FFT(r, :));
    end
    
    % RD-Map 생성
    RD_Map = abs(frame_data_2D_FFT);
    RD_Map(:, nd/2+1) = 0; % DC 성분 추가 제거
    
    % % 선택 사항: CA-CFAR 적용 (주석 처리)
    % guard_cells = 4; % 보호 셀
    % training_cells = 8; % 학습 셀
    % PFA = 1e-3; % False Alarm 확률
    % threshold = zeros(size(RD_Map));
    % for r = guard_cells+training_cells+1:nr-(guard_cells+training_cells)
    %     for d = guard_cells+training_cells+1:nd-(guard_cells+training_cells)
    %         training_region = RD_Map(r-training_cells-guard_cells:r+training_cells+guard_cells, ...
    %                                  d-training_cells-guard_cells:d+training_cells+guard_cells);
    %         threshold(r, d) = mean(training_region(:)) * (PFA^(-1/(2*training_cells)) - 1);
    %     end
    % end
    % RD_Map = RD_Map > threshold; % 임계값 이상인 신호만 남김
    
    % RD-Map 저장
    RD_Map_3D(:,:,frame_number) = RD_Map;
    
    % 시각화
    imagesc(vd, rd, RD_Map);
    axis xy;
    axis tight;
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
    title(sprintf('RD-MAP Frame %d (Clutter Removed)', frame_number));
    colorbar;
    
    drawnow; % 창 업데이트
    % pause(0.1); % 프레임 전환 시간
    
    frame_number = frame_number + 1;
end

%% Range Profiling (클러터 제거 추가)
frame_number = 1;
rd = (0:nr-1) / nr * fs * c0/(2*kf); % Range 축
vd = (-nd/2:1:nd/2-1)/nd*PRF*c0/2/fc; % Doppler 축
td = (1:nframe)*tfr; % 시간축
RD_Map_3D = zeros(256, 128, 175); % RD-Map 저장 배열

% 클러터 제거를 위한 평균 계산
mean_frame_data = mean(data_3d, 3);

% High-Pass 필터 설계
[b, a] = butter(2, 0.1, 'high');

figure;
while frame_number < 176
    % 1개 Frame Data(256, 128) 갖고 오기
    frame_data = data_3d(:,:,frame_number);
    
    % 평균 제거 (클러터 제거)
    frame_data_clutter_removed = frame_data - mean_frame_data;
    
    % 열방향 FFT 진행 (Range FFT)
    frame_data_1D_FFT = fft(frame_data_clutter_removed, [], 1);
    frame_data_1D_FFT_Normalized = frame_data_1D_FFT / length(frame_data_1D_FFT(:, 1)); % FFT Normalize
    
    % 행방향 FFT 진행 (Doppler FFT)
    frame_data_2D_FFT = fftshift(fft(frame_data_1D_FFT_Normalized, [], 2), 2);
    
    % MTI 필터링
    for r = 1:nr
        frame_data_2D_FFT(r, :) = filtfilt(b, a, frame_data_2D_FFT(r, :));
    end
    
    % RD-Map 생성
    RD_Map = abs(frame_data_2D_FFT);
    RD_Map(:, nd/2+1) = 0; % DC 성분 제거
    RD_Map_3D(:,:,frame_number) = RD_Map;
    
    % Range Profiling
    Range_Profiling_1D = sum(RD_Map, 2); % Doppler 방향으로 합산
    Doppler_Profiling_2D(:, frame_number) = Range_Profiling_1D;
    
    % 시각화
    imagesc(td, rd, Doppler_Profiling_2D);
    axis xy;
    axis tight;
    xlabel('Time(s)');
    ylabel('Range (m)');
    title(sprintf('Range Profiling Frame %d (Clutter Removed)', frame_number));
    colorbar;
    
    drawnow; % 창 업데이트
    % pause(0.1); % 프레임 전환 시간
    
    frame_number = frame_number + 1;
end

%% Doppler Profiling (클러터 제거 추가)
frame_number = 1;
rd = (0:nr-1) / nr * fs * c0/(2*kf); % Range 축
vd = (-nd/2:1:nd/2-1)/nd*PRF*c0/2/fc; % Doppler 축
td = (1:nframe)*tfr; % 시간축
RD_Map_3D = zeros(256, 128, 175); % RD-Map 저장 배열

% 클러터 제거를 위한 평균 계산
mean_frame_data = mean(data_3d, 3);

% High-Pass 필터 설계
[b, a] = butter(2, 0.1, 'high');

figure;
while frame_number < 176
    % 1개 Frame Data(256, 128) 갖고 오기
    frame_data = data_3d(:,:,frame_number);
    
    % 평균 제거 (클러터 제거)
    frame_data_clutter_removed = frame_data - mean_frame_data;
    
    % 열방향 FFT 진행 (Range FFT)
    frame_data_1D_FFT = fft(frame_data_clutter_removed, [], 1);
    frame_data_1D_FFT_Normalized = frame_data_1D_FFT / length(frame_data_1D_FFT(:, 1)); % FFT Normalize
    
    % 행방향 FFT 진행 (Doppler FFT)
    frame_data_2D_FFT = fftshift(fft(frame_data_1D_FFT_Normalized, [], 2), 2);
    
    % MTI 필터링
    for r = 1:nr
        frame_data_2D_FFT(r, :) = filtfilt(b, a, frame_data_2D_FFT(r, :));
    end
    
    % RD-Map 생성
    RD_Map = abs(frame_data_2D_FFT);
    RD_Map(:, nd/2+1) = 0; % DC 성분 제거
    RD_Map_3D(:,:,frame_number) = RD_Map;
    
    % Doppler Profiling
    Doppler_Profiling_1D = sum(RD_Map, 1); % Range 방향으로 합산
    Doppler_Profiling_2D(:, frame_number) = Doppler_Profiling_1D;
    
    % 시각화
    imagesc(td, vd, Doppler_Profiling_2D);
    axis xy;
    axis tight;
    xlabel('Time(s)');
    ylabel('Velocity (m/s)');
    title(sprintf('Doppler Profiling Frame %d (Clutter Removed)', frame_number));
    colorbar;
    
    drawnow; % 창 업데이트
    pause(0.1); % 프레임 전환 시간
    
    frame_number = frame_number + 1;
end

% Doppler Profiling 시각화 (전체 프레임)
figure;
imagesc(td, vd, Doppler_Profiling_2D);
colormap('jet'); % 컬러맵을 'jet'으로 설정하여 신호 강도 대비를 높임
colorbar;
axis xy;
axis tight;
xlabel('Time(s)');
ylabel('Velocity (m/s)');
title('Doppler Profiling (Clutter Removed)');

% 특정 프레임에서의 Doppler 변화 분석 (프레임 113)
figure;
specific_frame = 113;
plot(vd, Doppler_Profiling_2D(:, specific_frame), 'b-', 'LineWidth', 2);
grid on;
xlabel('Velocity (m/s)');
ylabel('Amplitude');
title(sprintf('Doppler Profile at Frame %d (Clutter Removed)', specific_frame));
xlim([-5 5]); % Velocity 축 범위 설정