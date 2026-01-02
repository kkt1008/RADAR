%--------------------------------------------------------------------------
%   TITLE : Basic Simulation of SAR raw data and processing for beginners
%--------------------------------------------------------------------------
%   %    
%               1. Generation of RAW radar data
%               2. Range compression
%               3. Azimuth compression
%
%    Basic satellite SAR point target with zero squint angle
%
%--------------------------------------------------------------------------
%   Property of SATECS,KAU
%         only for educational purpose
%--------------------------------------------------------------------------
clear all
close all
clc

%--------------------------------------------------------------------------
% Paramater Load (Satellite)
R_0         = 20e3;         % [m]       Slant range of scene center: R(n)
Vr          = 150;          % [m/s]     Effective radar velocity
Kr          = 0.25e12;      % [Hz/sec]  Range FM rate
Tr          = 25e-6;        % [sec]     Transmitted pulse duration
BW          = 6.25e6        % [MHz]     Bandwidth
Fr          = 7.5e6;        % [Hz]      Range sampling rate
Fa          = 104;          % [Hz]      Azimuth sampling rate or PRF
N_az        = 256;          % [smp]     Number of range lines
N_rg        = 256;          % [smp]     Samples per range line
f_0         = 5.3e9;        % [Hz]      Center frequency
c           = 3e8;          % [m/s]     Velocity of light
lamda       = c / f_0;      % [m]       Wavelength
Ka          = 2*Vr^2/(lamda*R_0); % [Hz/sec]  Azimuth FM rate
BeamWidth   = 0.015;        % 3-dB beamwidth (radian)

%load geometry and antenna parameters 

eta_start = - N_az/2 * 1/Fa;         % azimuth start time
eta_end = N_az/2 * 1/Fa;             % azimuth end time

eta       = eta_start : 1/Fa : eta_end - 1/Fa;  % azimuth time variable
R_eta     = sqrt(R_0^2 + Vr^2 * eta.^2);    % range equation

% Data accumuation of incoming radar signals 
tau_data = - Tr/2 + (2*R_0)/c;           % range time
N_S_rg = fix(Tr * Fr);

for K = 1 : N_S_rg                  % Number of Sample in range direction
    s_0(:,K) = exp(-j*4*pi*f_0*(R_eta)./c) ...
            .* exp(j*pi*Kr*(tau_data - 2*R_eta/c).^2);
    tau_data = tau_data + 1/Fr;               % 1/Fr : sampling interval
end

pulse_start = fix((N_rg-N_S_rg)/2);
pulse_end = pulse_start + N_S_rg;

temp = zeros(N_az,N_rg);
temp(:,pulse_start:pulse_end-1) = s_0;
s_0 = temp;

%--------------------------------------------------------------------------
% Generation of Matched filter
tau_start = -Tr/2;                         % Start time of pulse
tau_end = +Tr/2;                         % End time of pulse
tau = tau_start : 1/Fr : tau_end - 1/Fr;    % range time variable

eta_start = - N_az/2 * 1/Fa;               % start time of azimuth
eta_end = + N_az/2 * 1/Fa;               % end time of aizmuth
eta    = eta_start : 1/Fa : eta_end - 1/Fa; % azimuth time variable

M_filter_az = zeros(N_az,1);                   % matched filter - azimuth
M_filter_rg = zeros(1,N_rg);                   % matched filter - range

M_filter_rg(1     , pulse_start:pulse_end-1) = exp(- j*pi*Kr*( tau ).^2 );
M_filter_az(1:N_az, 1                  ) = exp(  j*pi*Ka*( eta ).^2);

%--------------------------------------------------------------------------
% Compression process in frequency domain - range direction
S_rc = zeros(N_az,N_rg);        % Matched filtered signal
H_RG = fft(M_filter_rg);        % Range replica

temp_Ra = fft(s_0, [], 2);      % range-frequency (azimuth-time) domain

for I = 1:N_az
    S_rc(I,:) = temp_Ra(I,:) .* H_RG;   % Range Matched Filtering
end

s_rc = ifft(S_rc, [], 2);       % 2-D time domain

%--------------------------------------------------------------------------
% Compression process in frequency domain - azimuth direction
S_ac     = zeros(N_az,N_rg);    % Matched filtered signal - Azimuth
H_AZ    = fft(M_filter_az);            % Azimuth replica

S_rD = fft(s_rc, [], 1);

for I=1:N_rg
    S_ac(:,I) = S_rD(:,I).*H_AZ;    % Azimuth matched filtering
end

s_ac = fftshift(ifft(S_ac, [], 1)); % Azimuth compressed signal

%--------------------------------------------------------------------------
% Ploting the generated raw data

figure(1)
colormap(gray(256))
image(abs(s_0));
xlabel('range [sample]');
ylabel('azimuth [sample]'); title('Raw Magnitude');

figure(2)
s_0_ang = angle(s_0);
plot(s_0_ang(10,:))
xlabel('range [sample]');
ylabel('Phase [Rad]'); title('Raw Data Phase in range');

figure(3)
plot(s_0_ang(:,128))
xlabel('Azimuth [sample]');
ylabel('Phase [Rad]'); title('Raw Data Phase in Azimuth');

%--------------------------------------------------------------------------
% Range Compressed data
s_rc_comp = abs(fftshift(s_rc));
figure(4)
colormap(gray(256))
image(s_rc_comp);
xlabel('range time [sample]');  ylabel('azimuth time [sample]');
title('Range compressed signal - magnitude');

%--------------------------------------------------------------------------
% Ploting the compressed raw data
s_ac_comp = abs(s_ac);
figure(5)
colormap(gray(256))
image(s_ac_comp);
xlabel('range time [sample]');  ylabel('azimuth time [sample]');
title('Compressed point target');


