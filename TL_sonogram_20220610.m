% Spectlogram function
[libPath,outputPath]=definePaths(1);
addpath(libPath)
%% Time frequency vectors definition
% 
% lent=2^22;
% tWind=1e-9;%second. Actual time array=+-time_window/2


lent=2^21;
tWind=100e-9;                                                              %second. Actual time array=+-time_window/2

% t=linspace(0,tWind,lent);
t=linspace(-tWind/2,tWind/2,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling signal definiton %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Time Lens Approach

% Find tl from maximum sampling rate of AWG
%  m=10;%ntl/divAmount;
% tltly=m/65e9;

% Follwing settings hold for both the TL and TAI approaches
% Design tl
 tltry=0.2e-9;%6/65e9;%e-12;
% phimax_pi_units=7;%8;%2.5;%eta*pi/4/pi
% %% find b2 from eta and tl
%  eta=phimax_pi_units*4;%10; % number of points per window
% Find tl from eta and phi2
% phi=2651e-24;%720*2.1823e-23;
% eta=16;
% tltly=sqrt(eta*phi*2*pi)

% optimize tl to be integer of dt
ntl=round(tltry/dt);
tl=ntl*dt;
nSamps=ceil(lent/ntl);
tSingleSamp=dt*(1:ntl);
lensInds=0:ntl:lent;

% Optimize eta so that eta/tl is integer of df
etaini=7;

eta=round(etaini/(tl*df))*tl*df;
nFreqMax=eta/tl/df
% nfreqMax=(eta/(ntl*dt^2
% nfreqMax=round(eta/tl/df)
% tl=eta/df/nfreqMax
fmax=nFreqMax*df;

%% Time Lens - Uncomment below for TL and comment "Talbot TAI section"
% 
%% Find eta from b2 and tl
% phi=50*2.1823e-23;
% % eta=(tl^2)/(phi*2*pi)

% C_tl=tl^-1*120e9;
C_tl=2*pi*eta/(tl^2);
singleSamp=C_tl/2*(tSingleSamp-tSingleSamp(round(ntl/2))).^2;
% sampSig=repmat(singleSamp,1,nSamps);
sampSigRaw=repmat(singleSamp,1,nSamps);
sampSigRawLent=sampSigRaw(1:lent);
sampSig=sampSigRawLent;%real(filtSG_tf(sampSigRawLent,t,f,round(800e9/df),5,0));
 sampSig=sampSig(1:lent);
sampFreq=1/tl
tl*C_tl/(2*pi)

%% Talbot TAI - Uncomment below for TAI and comment "Time Lens"
% 
% 
% m=round(eta);
% ts=tl/m
% AWG_nuq=1/ts
% p=1;
% s=-(m-1);%generateSparameter(p,m);
% GV=wrapTo2Pi(s/m*pi*((0:m-1)).^2); 
% nSampsPer_ts=ceil(ntl/m)
% GVtly=repelem(GV,1,nSampsPer_ts);
% t_samples=linspace(tSingleSamp(1),tSingleSamp(end),numel(GVtly));%(1:numel(GVtly))*dt;%/numel(GVtly)*tl;
% singleSamp=interp1(t_samples,GVtly,tSingleSamp);
% allSamps=repmat(singleSamp,1,nSamps);
% allSamps=allSamps(1:lent);
% sampSig=real(filtSG_tf(allSamps,t,f,round(200e9/df),10,1));
% 
% 

%%%%%%%%%%%%%%%%
%% Dispersion %% 
%%%%%%%%%%%%%%%%


%% Sample (obsolete)
% phi=tl^2/(2*pi);
%% Time Lens 
phi=1/C_tl % Get exact phi for simulations
% TAI
% phi=p/m*(1/2*pi)*(tl/m)^2;%p*m*(tl/m)^2/(2*pi);

 phi2perKm=   2.1823e-23;
 phi/phi2perKm

% Leave below uncommented for all cases.
  phi2perKm=   2.1823e-23; 
  dispersive_smfKm=phi/phi2perKm
phaseGVD=phi/2*(2*pi*f).^2;



%%%%%%%%%%%%%%%%%%
%% Generate SUT %%
%%%%%%%%%%%%%%%%%%
SUT=ones(1,lent);
% 
% 
% 
% % fSUT=0.25*(C_tl*tl)/(2*pi);
% fmax=eta/tl;
% fSUT=linspace(0,fmax,numel(t));
% SUT=(cos(2*pi*fSUT.*t)).*(superGauss(0,tWind/8,1,t,tWind/6)+superGauss(0,tWind/8,2,t,-tWind/6));
% SUT_f=nfft(SUT,dt);
% 
BW=Fs/20;
tExt=tl*0.3;
phi2=tExt/(2*pi*BW);
SUT_f=superGauss(0,BW,2,f,0).*exp(1j*phi2*(2*pi*f).^2/2);%+1j*phi3*(2*pi*f).^3/6);
SUT=nifft(SUT_f,Fs);

% 

tChirp=tSingleSamp;
% SUTfreqs=linspace(-2e12,2e12,numel(tChirp));
f1=0.01e12; f0=fmax*20;%0.1e12;
SUT=zeros(1,lent);
SUT(round((lent/2-ntl/2):(lent/2+ntl/2-1)))=exp(1j*(2*pi*((f1-f0)/(2*tl)*tChirp.^2+f0*tChirp)));%+ones(1,ntl);
SUT_f=nfft(SUT,dt);tl




% 
% % fSUT=0.25*(C_tl*tl)/(2*pi);
% fSUT=linspace(0,1e3*(C_tl*tl)/(2*pi)/4,numel(t));
% SUT=(cos(2*pi*fSUT.*t)).*(superGauss(0,tWind/8,1,t,tWind/6)+superGauss(0,tWind/8,2,t,-tWind/6));
% 
% fSUT=linspace(0,1e2*(C_tl*tl)/(2*pi)/4,ntl);
% SUTAction=(cos(2*pi*fSUT.*tSingleSamp));%.*(superGauss(0,tWind/8,1,t,tWind/6)+superGauss(0,tWind/8,2,t,-tWind/6));
% SUT=zeros(1,lent);
% SUT(lent/2-round(ntl/2):lent/2+round(ntl/2)-1)=SUTAction;
% SUT_f=nfft(SUT,dt);

% SUT_f=superGauss(0,3e12,10,f,0).*exp(1j*1*phi2perKm*(2*pi*f).^2);%
% SUT=nifft(SUT_f,Fs);
% 
% BW=3e12;%Fs/20;
% tExt=tl*2;
% phi2=tExt/(2*pi*BW);
% SUT_f=superGauss(0,BW,2,f,0).*exp(1j*phi2*(2*pi*f).^2/2);%+1j*phi3*(2*pi*f).^3/6);
% SUT=nifft(SUT_f,Fs);


figure;
subplot(2,1,1)
plot(f,abs(SUT_f).^2)
subplot(2,1,2)
plot(t,abs(SUT).^2)
hold on
plot(t,real(SUT))
% % 
% 
% figure;
% subplot(2,1,1)
% plot(f*1e-12,abs(SUT_f))
% hold on 
% plot(f*1e-12,real(SUT_f))
% xlabel('Frequency (THz)')
% subplot(2,1,2)
% plot(t*1e12,abs(SUT).^2)
% xlabel('Time (ps)')
% % 
%%%%%%%%%%%%%%%%
%% Processing %%
%%%%%%%%%%%%%%%%
% 
% 


 sampSUT=exp(1j*sampSig).*SUT   ;
sampSUT_f=nfft(sampSUT,dt);

sampSUTdisp_f=(sampSUT_f).*exp(1j*phaseGVD);
sampSUTdisp=nifft(sampSUTdisp_f,Fs);

%%% short inverse fourier tlansform
sampSUT_sift_disp=exp(1j*phaseGVD).*SUT_f;
sampSUT_sift_disp_t=nifft(sampSUT_sift_disp,Fs);

sampSUT_sift_out_t=sampSUT_sift_disp_t.*exp(1j*sampSig);
sampSUT_sift_out=nfft(sampSUT_sift_out_t,dt);







%% Filter Spectlogram Signal
% sampSUTdisp=(filtSG_tf(sampSUTdisp,t,f,round((10/tl)/df),5,0));
% sono=reshape(circshift(sampSUT_sift_out(1:floor(lent/round(nFreqMax))*round(nFreqMax)),round(nFreqMax/2)-5500),round(nFreqMax),floor(lent/nFreqMax));
sono=reshape(circshift(sampSUT_sift_out(1:floor(lent/round(nFreqMax))*round(nFreqMax)),200),round(nFreqMax),floor(lent/nFreqMax));

%%% Reshape output signal to get 2D representation
% sono=reshape(circshift(sampSUT_sift_out(1:floor(lent/round(nFreqMax))*round(nFreqMax)),0),round(nFreqMax),floor(lent/nFreqMax));
figure;imagesc(abs(sono).^2)
figure;plot(f,abs(sampSUT_sift_out).^2)
f_spec=(tSingleSamp-max(tSingleSamp)/2)/phi/(2*pi);
t_spec=linspace(-tWind/2,tWind/2,numel(lensInds)-1);

%%% Reshape output signal to get 2D representation
% sptgm=reshape(circshift(sampSUTdisp(1:lensInds(end)),0),ntl,(numel(lensInds)-1));
% f_spec=(tSingleSamp-max(tSingleSamp)/2)/phi/(2*pi);
% t_spec=linspace(-tWind/2,tWind/2,numel(lensInds)-1);


t_spec=linspace(-fmax/2,fmax/2,nFreqMax)*2*pi/C_tl;
f_spec=linspace(f(1),f(floor(lent/nFreqMax)*nFreqMax),floor(lent/nFreqMax));%:(floor(lent/nFreqMax))*
% f_spec=(tSingleSamp-max(tSingleSamp)/2)/phi/(2*pi);
% t_spec=linspace(-tWind/2,tWind/2,numel(lensInds)-1);
figure;imagesc(f_spec*1e-12,t_spec*1e9,abs(sono).^2)
set(gca,'YDir','normal')
xlabel('Frequency (THz)')
ylabel('Time (ns)')
figure;plot(f,abs(sampSUT_sift_out).^2)
figure;imagesc(abs(sono).^2)


figure;
subplot(3,1,1)
imagesc(f_spec*1e-12,t_spec*1e12,abs(sono).^2)
set(gca,'YDir','normal')
% legend('SSIFT')
xlabel('Frequency (THz)')
ylabel('Time (ps)')
subplot(3,1,2)
plot(f*1e-12,abs(sampSUT_sift_out))
hold on
plot(f*1e-12,real(sampSUT_sift_out))
legend('SUT','Spectlal tlace SSIFT')
xlabel('Frequency (THz)')
subplot(3,1,3)
plot(t*1e12,abs(SUT).^2)
yyaxis right
plot(t*1e12,abs(sampSUT_sift_out).^2)
legend('SUT','Temporal tlace SSIFT')
xlabel('Time (ps)')
% figure;imagesc(abs(sono).^2)


% 
% %% Reshape SUT for FFT
% fSingleSamp=linspace(-Fs/2,Fs/2,numel(tSingleSamp));
% padLen=2^17;
%  SUT2D=reshape(circshift(SUT(1:lensInds(end)),0),ntl,(numel(lensInds)-1));
% SUT2DPad=[zeros(padLen,numel(SUT2D(1,:))); SUT2D; zeros(padLen,numel(SUT2D(1,:)))];
% tPad=((1:numel(SUT2DPad(:,1)))-numel(SUT2DPad))*dt;
% 
% onesPad=ones(padLen,numel(SUT2D(1,:)));
% SUTiniBits=onesPad+SUT2D(1,:);
% SUTfinBits=onesPad+SUT2D(end,:);
% SUT2DPad=[SUTiniBits; SUT2D; SUTfinBits];
% figure; subplot(2,1,1)
% plot(tPad,real(SUT2DPad)); hold on; plot(tPad,imag(SUT2DPad))
% subplot(2,1,2)
% fSingleSamp=linspace(-Fs/2,Fs/2,numel(SUT2DPad(:,1)));
% plot(fSingleSamp,abs(fftshift(fft(ifftshift(SUT2DPad)))))

%% Plotting section

tlims=[-200 200]*tl+7.511e-8;
ylims=([-0.1 1.1]);
ylimsNegPos=([-1.1 1.1]);
numtlWindLims=20;
tlims='auto';%(t(end)/2+[-numtlWindLims*tl numtlWindLims*tl])*1e9;
% tlims=[t(1) t(end)];

tps=1e12; tns=1e9; tus=1e6; 
lbps='Time (ps)'; lbns='Time (ns)'; lbus='Time (us)';
fG=1e-9; lbG='Frequency (GHz)';
figure('Renderer', 'painters', 'Position', [50 50 1200 700])
% sampling signal in the time domain





subplot(4,4,1:4)
plot(t*tns,abs(SUT).^2/(max(abs(SUT).^2)))
hold on
plot(t*tns,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
ylabel('Intensity')
yyaxis right
plot(t*tns,sampSig)
%  xlim([t(1) t(end)]*tns)
 xlim(tlims);%ylim(ylims);
xlabel(lbns)
legend('sampled SUT','Dispersed wvf','Time Lens')
ylabel('Phase (rad)')



subplot(4,4,5:8)
plot(t*tns,abs(SUT));%/(max(abs(SUT))))
hold on
plot(t*tns,angle(SUT));%/(max(angle(SUT))))

plot(t*tns,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
yyaxis right
plot(t*tns,sampSig)

ylabel('Intensity')
% yyaxis right
% plot(t*tns,sampSig)
 xlim(tlims);%ylim(ylimsNegPos);
xlabel(lbns)
legend('abs','angle','spectlogram')
ylabel('Phase (rad)')

subplot(4,4,9:12)
% plot(t*tns,real(SUT)/(max(real(SUT))))
% hold on
% plot(t*tns,imag(SUT)/(max(imag(SUT))))

plot(t*tns,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
ylabel('Intensity')
% yyaxis right
% plot(t*tns,sampSig)
 xlim(tlims);%ylim(ylimsNegPos);
xlabel(lbns)
legend('real','imag','spectlogram')
ylabel('Phase (rad)')

subplot(4,4,13:16)

    imagesc(t_spec*tns,f_spec*fG,(abs(sono).^2));
    hcb=colorbar();
     xlim(tlims)

    ylabel(hcb,'Power (dB)')
    xlabel('Time (ns)')
    
%% fold over spectlogram for "eye diagram"



% subplot(4,4,5:6)
% plot(t*tns,abs(sampSUT).^2/(max(abs(sampSUT).^2)))
% hold on
% plot(t*tns,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
% ylabel('Intensity')
% yyaxis right
% plot(t*tns,sampSig)
% xlimProp=1/8;
%  xlim([xlimProp*t(end)*tns xlimProp*t(end)*tns+1.2*tl*tns])
% xlabel(lbns)
% 
% % legend('sampled SUT','Dispersed wvf','Time Lens')
% 
% subplot(4,4,7:8)
% plot(t*tns,abs(sampSUT).^2/(max(abs(sampSUT).^2)))
% hold on
% plot(t*tns,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
% yyaxis right
% plot(t*tns,sampSig)
% xlimProp=5/6;
%  xlim([xlimProp*t(end)*tns xlimProp*t(end)*tns+1.2*tl*tns])
% xlabel(lbns)
% ylabel('Phase (rad)')
% 
% 
% subplot(4,4,9:12)
% % if logPlot
%     imagesc(t_spec*tus,f_spec*fG,10*log10(abs(sptgm).^2));
%     caxis([-35 -0.5])
%     hcb=colorbar();
%     ylabel(hcb,'Power (dB)')
% % else
% %     imagesc(t_spec*tus,f_spec*fG,(abs(sptgm).^2));
% %     hcb=colorbar()
% %     ylabel(hcb,'Power (mW)')
% % end
% % xlim(tlims*tus)
% xlabel(lbus); ylabel(lbG);
% set(findall(gcf,'-property','FontSize'),'FontSize',16)
% 
% 
% subplot(4,4,13:16)
% % if logPlot
% %     imagesc(t_spec*tus,f_spec*fG,10*log10(abs(sptgm).^2));
% %     caxis([-35 -0.5])
% %     hcb=colorbar()
% %     ylabel(hcb,'Power (dB)')
% % else
%     imagesc(t_spec*tus,f_spec*fG,10*log10(abs(sptgm).^2));
%     hcb=colorbar();
%     ylabel(hcb,'Power (dB)')
% % %     imagesc(t_spec*tus,f_spec*fG,(abs(sptgm).^2));
% %     hcb=colorbar();
% %     ylabel(hcb,'Power (mW)')
% % end
% % xlim(tlims*tus)
% xlabel(lbus); ylabel(lbG);
% set(findall(gcf,'-property','FontSize'),'FontSize',16)
% 


function waveform=superGauss(C,t0,m,xs,center)
waveform=exp(-(1+1j*C)/2*((xs-center)/t0).^(2*m));
end

function fftout=nfft(sig,varargin)
%fft gives the swaped spectlum.
unnormdfft=fftshift(fft(ifftshift(sig)));

scale=nargin-1;

if scale==1
    dt=varargin{1};
fftout=dt*unnormdfft;
%     else if scale==2
%     fftout=1/sqrt(length(sig))*unnormdfft;
        else
        fftout=1/(max(abs(unnormdfft)))*unnormdfft;
%         end
end
%fftout=dt*unnormdfft;
end

function ifftout=nifft(sig,varargin)
%%%Need to first ifftshift to swap for ifft (Matlab assumes swaped fft)
%normalize if varargin is given

unnormdifft=fftshift(ifft(ifftshift(sig)));

scale=nargin-1;

if scale==1
    Fs=varargin{1};
ifftout=Fs*unnormdifft;
% else if scale==2
%     ifftout=1/sqrt(length(sig))*unnormdifft;
else
ifftout=1/(max(abs(unnormdifft)))*unnormdifft;

end
%ifftout=df*unnormdifft;
end

function s=generateSparameter(p,q)

if mod(q,2)==1
    s=2*mulinv(2*p,q);
else
    s= mulinv(p,2*q) ;
end

end


function y=mulinv(x,p)
% 
% if ~isprime(p)
%     disp('The field order is not a prime number');
%     return
% elseif sum(x>=p)
%     disp('All or some of the numbers do not belong to the field');
%     return
% elseif sum(~x)
%     disp('0 does not have a multiplicative inverse');
%     return
% end

k=zeros(size(x));   %set the counter array
m=mod(k*p+1,x);     %find reminders
while sum(m)        %are all reminders zero?
    k=k+sign(m);    %update the counter array
    m=mod(k*p+1,x); %caculate new reminders 
end
y=(k*p+1)./x;       %find the multiplicative inverses of X
end

function S_OS=digiMod()
%% Digital modulation parameters

P = 15;                                                  % PRBS order
L = 2^P-1;                                              % PRBS length

digimod_format = 'SP-PAM';                                 % Digital modulation format
M = 4;                                                 % Number of symbols (modulation order)

K = log2(M);                                            % Number of bits per symbol
N = L*K;                                                % Number of bits
B = idinput([N,2],'prbs',[1,1],[0,1]).'; B = B(1,:);    % Bit sequence
D = bi2de(reshape(B,L,K)).';                            % Symbol indices


%% Symbol sequence

switch digimod_format
    case 'SP-PAM'           % Signle-polarity PAM
        S = D;
    case 'DP-PAM'           % Dual-polarity PAM
        S = pammod(D,M);
    case 'PSK'              % PSK
        S = pskmod(D,M);
    case 'DPSK'             % Differential PSK
        S = dpskmod(D,M);
    case 'QAM'              % QAM
        S = qammod(D,M);
end

S_I = real(S);              % In-phase channel
S_Q = imag(S);              % Quadrature channel

S_M = abs(S);               % Magnitude
S_P = angle(S);             % Phase


%% Oversampling

OSR = 1;                                   % Oversampling rate

S_OS = S(ones(1,OSR),:); S_OS = S_OS(:).';  % Oversampled symbol sequence

S_OS_I = real(S_OS);                        % Oversampled In-phase channel
S_OS_Q = imag(S_OS);                        % Oversampled Quadrature channel

S_OS_M = abs(S_OS);                         % Oversampled Magnitude
S_OS_P = angle(S_OS);                       % Oversampled Phase


%% Results

figure
subplot(2,1,1)
bar(S_OS_I)
subplot(2,1,2)
bar(S_OS_Q)

figure
subplot(2,1,1)
bar(S_OS_M)
subplot(2,1,2)
bar(S_OS_P/pi)

figure
plot(S_I,S_Q,'o')
axis square


end
