

% Chirp Config
nrx = 4; % number of receivers
startFreq       =   77;
idleTime        =   100;
adcStartTime    =   6;
rampEndTime     =   60;
digOutSampleRate=   4270;
freqSlopeConst  =   66.143;
numADCSamples   = 256; % number of ADC samples per chirp

% Frame Config
chirpStartIdx   =   0;
chirpEndIdx     =   1;
numLoops        =   128;
framePeriodicity=   40.0;

% ADC Config
numADCBits      = 16; % number of ADC bits per sample
isReal = 1; % set to 1 if real only data, 0 if complex data0

% Radar Parameters
c0              =   299792458;
nrx             =   1;
ntx             =   chirpEndIdx - chirpStartIdx + 1;
nd              =   numLoops;
nr              =   numADCSamples;    
tfr             =   framePeriodicity * 1e-3;
fs              =   digOutSampleRate * 1e3;
kf              =   freqSlopeConst * 1e12;
adc_duration    =   nr / fs;
BW              =   adc_duration * kf;
PRI             =   (idleTime + rampEndTime) * 1e-6;
fc              =   startFreq * 1e9 + kf * (adcStartTime * 1e-6 + adc_duration / 2); 
nfr             =   2;

RMin            =   1;
RMax            =   10;
NFFT            =   2^(nextpow2(nr)+2);
NFFTA           =   2^(nextpow2(nrx*ntx)+2);
NFFTD           =   2^(nextpow2(nd)+2);

vRange          =   [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);
vAngDeg         =   linspace(-60,60,NFFTA);
vFreqVel        =   [-NFFTD./2:NFFTD./2-1].'./NFFTD.*(1/PRI);
vVel            =   vFreqVel*c0/(2.*fc);

[Val RMinIdx]   =   min(abs(vRange - RMin));
[Val RMaxIdx]   =   min(abs(vRange - RMax));
vRangeExt       =   vRange(RMinIdx:RMaxIdx);


vU                  =   sin(vAngDeg/180*pi);
[mRange , mU]       =   ndgrid(vRangeExt,vU);
mX                  =   mRange.*mU;
mY                  =   mRange.*cos(asin(mU));

adc_file_name = 'human_azimuth.bin';
data = readDCA1000(adc_file_name);
disp('data was generated!')

data1=data(2,1:end);  % channel 1

data_sample = data(2,1:256);
t=[0:1/fs:(4.0e-5)-1/fs];

data_2d=reshape(data1,[nr,length(data1)/nr]); %data1을 256, length(data1)/256으로 나눠서 2차원행렬로 바꿈(256,8192)
[i,j]=size(data1);
[y,x]=size(data_2d);

figure(1)
yy=(1:1:y)*(vRange(2)-vRange(1));
plot(yy,20*log((abs(fft(data_sample))))); %data의 2행 1번 sample하나 빼서 dB Scale로 변경 후 plot
nframe=j/(numADCSamples^2);
axis tight;
xlabel('range(m)')
ylabel('amplitude')

data_2d_abs = abs(data_2d); 
data_fft=fft(data_2d,[],1); %data_2d에 대해 열방향으로 fft진행

figure(2)
data_fft_abs = abs(data_fft);
xx=(1:1:x)*256/fs*10;
mesh(xx,yy,data_fft_abs);view(0,90);
axis tight;
xlabel('time(s)')
ylabel('range(m)')

figure(3)
data_fft_2=fftshift(fft(data_fft,[],2),2); %data_fft를 행방향으로 fft 진행하여 frequency doppler 찾기
data_fft_2_abs = abs(data_fft_2);
[m,n]=size(data_fft_2_abs);
vd=(-x/2:1:x/2-1)/x*digOutSampleRate*c0/2/fc; 
mesh(vd,yy,data_fft_2_abs);view(0,90);
axis tight;
xlabel('velocity(m/s)')
ylabel('range(m)')


            