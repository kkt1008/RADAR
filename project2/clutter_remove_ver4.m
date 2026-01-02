% Chirp Config
nrx = 4; % Number of rx Antenna
startFreq = 77; % Start Frequency : 77 [GHz]
idleTime = 100; % IDLE 시간 [us]
adcStartTime = 6; % ADC 시작 시간 [us]
rampEndTime = 60; % Chirp 전체 시간 [us]
digOutSampleRate = 4270; % Sampling Frequency(Fs) [kHz]
freqSlopeConst = 66.143; % 주파수 Sweep 속도(S) [MHz/us]
numADCSamples = 256; % Number of ADC samples per chirp

% Frame Config
chirpStartIdx = 0;
chirpEndIdx = 1;
numLoops = 128; % Number of Chirps per Frame
framePeriodicity = 40.0; % Frame 주기 [ms]

% ADC Config
numADCBits = 16; % Number of ADC bits per sample
isReal = 1; % Real only data

% Radar Parameters
c0 = 299792458; % Speed of light
nrx = 1;
ntx = chirpEndIdx - chirpStartIdx + 1; % Number of tx Antenna
nd = numLoops; % Number of Doppler FFT bins
nr = numADCSamples; % Number of Range FFT bins
tfr = framePeriodicity * 1e-3; % Frame period [s]
fs = digOutSampleRate * 1e3; % Sampling frequency [Hz]
kf = freqSlopeConst * 1e12; % Frequency slope [Hz/s]
adc_duration = nr / fs; % Chirp duration [s]
BW = adc_duration * kf; % Bandwidth [Hz]
PRI = (idleTime + rampEndTime) * 1e-6; % Pulse Repetition Interval [s]
PRF = 1 / PRI; % Pulse Repetition Frequency [Hz]
fc = startFreq * 1e9 + kf * (adcStartTime * 1e-6 + adc_duration / 2); % Center frequency [Hz]
nfr = 2;
frame_duration = PRI * numLoops;
frame_idletime = tfr - frame_duration;
Rmax = fs * c0 / (2 * kf); % Maximum range [m]

% Load RawData
try
    adc_file_name = 'human_azimuth.bin';
    data = readDCA1000(adc_file_name);
    if size(data, 2) ~= 5734400
        error('Data size mismatch!');
    end
    disp('data was generated!');
catch
    disp('Error loading data!');
    return;
end

% Data reshaping
data_1d = data(2, 1:end); % Channel 1 data
data_2d = reshape(data_1d, [nr, length(data_1d) / nr]); % 256x22400
data_3d = reshape(data_2d, [nr, 128, 175]); % 256x128x175
[y, x] = size(data_2d); % [256, 22400]
Range_Index = (0:nr-1) / nr * Rmax; % Range axis

%% Range FFT
data_2d_1D_FFT = fft(data_2d, [], 1) / length(data_2d_1D_FFT); % Range FFT
data_2d_1D_FFT_abs = abs(data_2d_1D_FFT); % Absolute value

%% MTI and Filtering (Optimized)
target_velocity = 1; % Target velocity [m/s]
window_size = round(PRF / (2 * target_velocity / (c0 / fc))); % Dynamic window size
window = hamming(max(5, min(20, window_size)))'; window = window / sum(window);
window_mean = conv2(data_2d_1D_FFT_abs, window(:), 'same');
data_mti_processed = data_2d_1D_FFT_abs - window_mean;
data_mti_processed(data_mti_processed < 0) = 0;
[~, max_indices] = max(data_mti_processed, [], 1);
data_mti_processed_new = zeros(size(data_mti_processed));
for col = 1:x
    start_idx = max(1, max_indices(col) - 5);
    end_idx = min(nr, max_indices(col) + 5);
    data_mti_processed_new(start_idx:end_idx, col) = data_mti_processed(start_idx:end_idx, col);
end
data_mti_processed_dB = 20*log10(data_mti_processed_new + eps);
data_3d_mti = reshape(data_mti_processed_dB, [nr, 128, 175]);

%% Doppler FFT
data_3d_doppler = fft(data_3d_mti, [], 2); % Doppler FFT
data_3d_doppler_abs = abs(data_3d_doppler);
data_3d_doppler_dB = 20*log10(data_3d_doppler_abs + eps);

%% Visualization
figure;
imagesc(Time_Index_All, Range_Index, 20*log10(data_2d_1D_FFT_abs + eps));
title('Before MTI (dB)');
axis xy; % Y축을 거리(상향)로 설정
axis tight; % 축 범위 조정
xlabel('Time (s)'); ylabel('Range (m)');
colorbar;

figure;
imagesc(Time_Index_All, Range_Index, data_mti_processed_dB);
title('After MTI (dB)');
axis xy; % Y축을 거리(상향)로 설정
axis tight; % 축 범위 조정
xlabel('Time (s)'); ylabel('Range (m)');
colorbar;
