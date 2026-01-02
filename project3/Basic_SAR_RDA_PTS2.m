%--------------------------------------------------------------------------
%   제목 : 초보자를 위한 Sentinel-1 SAR 원시 데이터 및 처리 기본 시뮬레이션
%--------------------------------------------------------------------------
%   %    
%               1. 원시 레이더 데이터 생성
%               2. 거리 압축
%               3. 방위 압축
%
%    Sentinel-1 IW 모드 기반 SAR 포인트 타겟 시뮬레이션
%
%--------------------------------------------------------------------------
%   소유권: SATECS, KAU
%         교육 목적으로만 사용
%--------------------------------------------------------------------------
clear all
close all
clc

%--------------------------------------------------------------------------
% 파라미터 로드 (Sentinel-1 IW 모드)
R_0         = 8.0e5;            % [m]       장면 중심의 사선 거리: R(n), 약 800 km
Vr          = 7600;             % [m/s]     효과적인 레이더 속도 (근사치)
Kr          = 1.0e10;           % [Hz/sec]  거리 FM 비율 (대역폭/펄스 지속 시간)
Tr          = 5e-6;             % [sec]     전송 펄스 지속 시간
BW          = 50e6;             % [Hz]      대역폭 (IW 모드 기준)
Fr          = 30e6;             % [Hz]      거리 샘플링 비율
Fa          = 1678;             % [Hz]      방위 샘플링 비율 또는 PRF
N_az        = 1024;             % [smp]     거리 라인 수
N_rg        = 1024;             % [smp]     각 거리 라인당 샘플 수
f_0         = 5.405e9;          % [Hz]      중앙 주파수 (C-대역)
c           = 3e8;              % [m/s]     빛의 속도
lamda       = c / f_0;          % [m]       파장
Ka          = (2 * Vr^2 * cos(0)^2) / (lamda * R_0); % [Hz/sec]  방위 FM 비율
La          = 12.0;             % [m]       안테나 길이 (근사치)
BeamWidth   = 0.886 * lamda / La; % 3-dB 빔폭 (라디안)

% 기하학 및 안테나 파라미터 로드 

eta_start = - N_az/2 * 1/Fa;         % 방위 시작 시간
eta_end = N_az/2 * 1/Fa;             % 방위 종료 시간

eta       = eta_start : 1/Fa : eta_end - 1/Fa;  % 방위 시간 변수
R_eta     = sqrt(R_0^2 + Vr^2 * eta.^2);    % 거리 방정식

% 들어오는 레이더 신호 데이터 축적 
tau_data = - Tr/2 + (2*R_0)/c;           % 거리 시간
N_S_rg = fix(Tr * Fr);

for K = 1 : N_S_rg                  % 거리 방향의 샘플 수
    s_0(:,K) = exp(-j*4*pi*f_0*(R_eta)./c) ...
            .* exp(j*pi*Kr*(tau_data - 2*R_eta/c).^2);
    tau_data = tau_data + 1/Fr;               % 1/Fr: 샘플링 간격
end

pulse_start = fix((N_rg-N_S_rg)/2);
pulse_end = pulse_start + N_S_rg;

temp = zeros(N_az,N_rg);
temp(:,pulse_start:pulse_end-1) = s_0;
s_0 = temp;

%--------------------------------------------------------------------------
% 매치드 필터 생성
tau_start = -Tr/2;                         % 펄스 시작 시간
tau_end = +Tr/2;                         % 펄스 종료 시간
tau = tau_start : 1/Fr : tau_end - 1/Fr;    % 거리 시간 변수

eta_start = - N_az/2 * 1/Fa;               % 방위 시작 시간
eta_end = + N_az/2 * 1/Fa;               % 방위 종료 시간
eta    = eta_start : 1/Fa : eta_end - 1/Fa; % 방위 시간 변수

M_filter_az = zeros(N_az,1);                   % 방위 매치드 필터
M_filter_rg = zeros(1,N_rg);                   % 거리 매치드 필터

M_filter_rg(1     , pulse_start:pulse_end-1) = exp(- j*pi*Kr*( tau ).^2 );
M_filter_az(1:N_az, 1                  ) = exp(  j*pi*Ka*( eta ).^2);

%--------------------------------------------------------------------------
% 주파수 영역에서 거리 방향 압축 과정
S_rc = zeros(N_az,N_rg);        % 매치드 필터링된 신호
H_RG = fft(M_filter_rg);        % 거리 레플리카

temp_Ra = fft(s_0, [], 2);      % 거리-주파수 (방위-시간) 영역

for I = 1:N_az
    S_rc(I,:) = temp_Ra(I,:) .* H_RG;   % 거리 매치드 필터링
end

s_rc = ifft(S_rc, [], 2);       % 2차원 시간 영역

%--------------------------------------------------------------------------
% 주파수 영역에서 방위 방향 압축 과정
S_ac     = zeros(N_az,N_rg);    % 방위 매치드 필터링된 신호
H_AZ    = fft(M_filter_az);            % 방위 레플리카

S_rD = fft(s_rc, [], 1);

for I=1:N_rg
    S_ac(:,I) = S_rD(:,I).*H_AZ;    % 방위 매치드 필터링
end

s_ac = fftshift(ifft(S_ac, [], 1)); % 방위 압축 신호

%--------------------------------------------------------------------------
% 생성된 원시 데이터 플로팅

figure(1)
colormap(gray(256))
image(abs(s_0));
xlabel('거리 [샘플]');
ylabel('방위 [샘플]'); title('원시 데이터 크기');

figure(2)
s_0_ang = angle(s_0);
plot(s_0_ang(10,:))
xlabel('거리 [샘플]');
ylabel('위상 [라디안]'); title('거리 방향 원시 데이터 위상');

figure(3)
plot(s_0_ang(:,512))
xlabel('방위 [샘플]');
ylabel('위상 [라디안]'); title('방위 방향 원시 데이터 위상');

%--------------------------------------------------------------------------
% 거리 압축 데이터
s_rc_comp = abs(fftshift(s_rc));
figure(4)
colormap(gray(256))
image(s_rc_comp);
xlabel('거리 시간 [샘플]');  ylabel('방위 시간 [샘플]');
title('거리 압축 신호 - 크기');

%--------------------------------------------------------------------------
% 압축된 원시 데이터 플로팅
s_ac_comp = abs(s_ac);
figure(5)
colormap(gray(256))
image(s_ac_comp);
xlabel('거리 시간 [샘플]');  ylabel('방위 시간 [샘플]');
title('압축된 포인트 타겟');
