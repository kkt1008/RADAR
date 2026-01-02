% ================================
% SAR SLC .tiff 파일 불러오기
% ================================
clc;
% (2) Tiff 객체로 열기
t = Tiff('tiff_data2.tiff','r');

% tiff 파일 구조 확인
info = imfinfo('tiff_data2.tiff'); 
disp(info(1))

rows = info.Height;       % 11952
cols = info.Width;        % 20771
bitdepth = info.BitDepth; % 32
datatype = 'int16';       % Sentinel-1 SLC는 대개 int16로 저장됨

fid = fopen('tiff_data2.tiff', 'r');
raw = fread(fid, [2*cols, rows], datatype);  % 2*cols, 각 픽셀마다 실수/허수 interleaved
fclose(fid);

raw = raw.'; % 행렬 transpose (rows x 2*cols)
real_part = raw(:,1:2:end);    % 홀수 인덱스: real
imag_part = raw(:,2:2:end);    % 짝수 인덱스: imag
data_cpx = double(real_part) + 1j*double(imag_part);  % 복소수 어레이 (행: row, 열: col)
%data_cpx = data_cpx(1:2048, 1:2048);
clear raw real_part imag_part

% (데이터 크기가 너무 크면, 부분 추출로 메모리 부담 완화)
% 예: data_cpx = data_cpx(1:2048, 1:2048);

%% 2. 시스템 파라미터 입력 (Sentinel-1 공식에 맞게 조정)
c = 3e8;                  % [m/s] 빛속도
Tr = 337e-6;               % [s] 펄스폭
BW = 36e6;              % [Hz] 대역폭
Fr = 2.160e6;               % [Hz] 레인지 샘플링 속도
Vr = 7.1e3;                 % [m/s] 레이다 플랫폼 속도
f_0 = 5.405e9;              % [Hz] 중심 주파수
lambda = c / f_0;         % [m] 파장
R_0 = 850e3;               % [m] slant range
Ka = 2*Vr^2/(lambda*R_0); % [Hz/sec] 방위 FM rate

[N_az, N_rg] = size(data_cpx); % Azimuth, Range 크기

Kr = BW/Tr;                % Range FM rate [Hz/sec]
N_S_rg = fix(Tr*Fr);       % 유효 거리 샘플 개수
pulse_start = fix((N_rg-N_S_rg)/2);
pulse_end = pulse_start + N_S_rg;
tau = (-Tr/2):1/Fr:(Tr/2-1/Fr); % Range time
%% 실제 slc 복소수 영상 시각화
% 진폭(amplitude)영상
figure;
imagesc(log10(abs(data_cpx) + eps));
colormap(gray); axis image; colorbar;
title('Sentinel-1 SLC Intensity (log10 scale)');
xlabel('Range (pixel)');
ylabel('Azimuth (pixel)');

% 위상 (phase)영상
figure;
imagesc(angle(data_cpx));
colormap(hsv); axis image; colorbar;
title('Sentinel-1 SLC Phase');
xlabel('Range (pixel)');
ylabel('Azimuth (pixel)');

%% 3. 거리(레인지) 압축 (Matched Filtering)
% Range 방향 matched filter impulse response
h_rg = zeros(1,N_rg);
h_rg(1, pulse_start+1:pulse_end) = exp(-1j*pi*Kr*(tau).^2);
H_RG = fft(h_rg, N_rg);

% 데이터의 각 azimuth line에 대해 range matched filtering
data_range = zeros(N_az, N_rg);
for I = 1:N_az
    temp = fft(data_cpx(I,:), N_rg);
    data_range(I,:) = ifft(temp .* H_RG, N_rg);
end

%% 4. 중심 Range Line의 방위 위상 추출 및 Linear FM 유사성 확인
% 중심 Range Line(거리 중심, col = N_rg/2)
center_rg = round(N_rg/2);
azimuth_line = data_range(:, center_rg);
az_phase = angle(azimuth_line);

figure;
plot(az_phase, '-');
xlabel('Azimuth sample');
ylabel('Phase [rad]');
title('Azimuth Direction Phase Profile (Range Compressed)');
grid on;

% (참고: 위상이 포물선/2차함수(Parabolic)이면 Linear FM 신호와 위상 유사)

%% 5. 방위(에지무스) 압축 및 점표적 영상 생성
% Azimuth matched filter impulse response
eta_start = -N_az/2 * 1/104;   % 104Hz: 예시 PRF, 실제 PRF로 대체
eta = eta_start : 1/104 : eta_start + (N_az-1)/104;
h_az = exp(1j*pi*Ka*eta.^2).';  % 방위 replica
H_AZ = fft(h_az, N_az);

S_rD = fft(data_range, [], 1);
data_az = zeros(N_az, N_rg);
for J = 1:N_rg
    data_az(:,J) = ifft(S_rD(:,J) .* H_AZ, N_az);
end

% 점표적 영상
point_image = abs(fftshift(data_az));
figure;
imagesc(point_image); colormap gray; colorbar;
xlabel('Range'); ylabel('Azimuth');
title('Compressed Point Target');

%% 6. 3dB 해상도 측정 (거리/방위 해상도)
% 최댓값 위치 찾기
[maxval, maxidx] = max(point_image(:));
[row_max, col_max] = ind2sub(size(point_image), maxidx);

% 거리 해상도
range_profile = point_image(row_max, :);
halfmax_r = max(range_profile)/sqrt(2); % -3dB
ind_r = find(range_profile >= halfmax_r);
range_res_pix = ind_r(end)-ind_r(1)+1;
range_resolution = c/(2*BW); % 이론 해상도 [m]

% 방위 해상도
az_profile = point_image(:, col_max);
halfmax_a = max(az_profile)/sqrt(2); % -3dB
ind_a = find(az_profile >= halfmax_a);
az_res_pix = ind_a(end)-ind_a(1)+1;
az_resolution = lambda/2; % 이론 해상도 [m] (stripmap 기준)

fprintf('Measured range -3dB width (pixels): %d\n', range_res_pix);
fprintf('Measured azimuth -3dB width (pixels): %d\n', az_res_pix);
fprintf('Theoretical range resolution: %.2f m\n', range_resolution);
fprintf('Theoretical azimuth resolution: %.2f m\n', az_resolution);


% % (3) TIFF 데이터 읽기 (대개 실수 or 복소수 magnitude)
% data_raw = double(read(t));  % 보통 uint16 또는 int16 → double로 변환
% close(t);
% 
% % 복소수면 abs 취해서 진폭 영상 생성
% amplitude = abs(data_raw);             % 진폭 영상
% 
% % 3. 로그 스케일 변환 (SAR 영상은 동적 범위가 넓기 때문에 로그 씌움)
% log_amp = 20 * log10(amplitude + eps);  % eps: 0 방지용
% 
% % 4. 영상 시각화
% imagesc(log_amp);
% colormap gray;
% axis image;
% colorbar;
% title('SAR Intensity Image (Log Scale)');
