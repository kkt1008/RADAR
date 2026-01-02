% ================================
% SAR SLC .tiff 파일 불러오기
% ================================
clear;
close all;
clc;

% (1) TIFF 파일 경로
tiffFile = 'tiff_data1.tiff';  

% (2) Tiff 객체로 열기
t = Tiff(tiffFile, 'r');

% (3) TIFF 데이터 읽기 (대개 실수 or 복소수 magnitude)
data_raw = double(fopen(tiffFile,'r'));  % 보통 uint16 또는 int16 → double로 변환
close(t);

% 복소수면 abs 취해서 진폭 영상 생성
amplitude = abs(data_raw);             % 진폭 영상

% ================================
% 조건부 진폭 변환
% ================================
amplitude_mod = amplitude;  % 원본 보존

% 50~100: 2배6
mask_50_100 = amplitude >= 50 & amplitude < 100;
amplitude_mod(mask_50_100) = amplitude(mask_50_100) * 50;

% 100~150: 4배
mask_100_150 = amplitude >= 100 & amplitude < 150;
amplitude_mod(mask_100_150) = amplitude(mask_100_150) * 3000;

% 150~200: 6배
mask_150_200 = amplitude >= 150 & amplitude < 200;
amplitude_mod(mask_150_200) = amplitude(mask_150_200) * 10000;

% 200 이상: 8배
mask_200up = amplitude >= 200;
amplitude_mod(mask_200up) = amplitude(mask_200up) * 20000;

% 0~50은 그대로(별도 처리 불필요)

% 3. 로그 스케일 변환 (SAR 영상은 동적 범위가 넓기 때문에 로그 씌움)
log_amp = 20 * log10(amplitude_mod + eps);  % eps: 0 방지용

%% 
% 4. 영상 시각화
imagesc(flipud(log_amp));
colormap gray;
axis image;
colorbar;
title('SAR Intensity Image (Log Scale, Conditional Amplitude)');
