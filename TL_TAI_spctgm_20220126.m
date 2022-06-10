% Spectrogram function
addpath( '/Users/ben/Documents/MATLAB/library' )
%% Time frequency vectors definition

lent=2^18;                      % Signal length
tWind=100e-9;                   % Time window span

% t=linspace(0,tWind,lent);
t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;

spectrogramTypes={'T-TAI','Time-lens'};
processUsingType=1; % Make =1 for T-TAI, =2 for Time lens
spectrgramNow=spectrogramTypes{processUsingType};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temporal phase signal definiton %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These settings are used for both the T-TAI and the time-lens spectrogram
trTry=70/65e9;        % Time-lens apperture
eta=70;                % number of analysis points
fmax=eta/trTry;         % full bandwidth that can be processed with the spectrogram

ntr=round(trTry/dt);    % tr is adjusted to be an integer amount of points
tr=ntr*dt;
nSamps=ceil(lent/ntr);
tSingleSamp=dt*(1:ntr);
lensInds=0:ntr:lent;


switch spectrgramNow
    
    case 'Time-lens'
        %% Time Lens Approach
        
        C_tl=2*pi*eta/(tr^2);
        singleSamp=C_tl/2*(tSingleSamp-tSingleSamp(round(ntr/2))).^2;
        sampSigRaw=repmat(singleSamp,1,nSamps);
        sampSigRawLent=sampSigRaw(1:lent);
        sampSig=real(filtSG_tf(sampSigRawLent,t,f,round(800e9/df),5,0));
        sampSig=sampSig(1:lent);
        sampFreq=1/tr
        tr*C_tl/(2*pi)
        
        %%%%%%%%%%%%%%%%
        %% Dispersion %%
        %%%%%%%%%%%%%%%%
        phi=1/C_tl;
        
    case 'T-TAI'
        
        %% Talbot TAI - Uncomment below for TAI and comment "Time Lens"
        
        m=eta; % number of analysis points is equivalent to "amplification factor"
        ts=tr/m;
        AWG_nuq=1/ts;
        p=1;
        s=generateSparameter(p,m);
        GV=wrapTo2Pi(s/m*pi*((0:m-1)).^2);
        nSampsPer_ts=ceil(ntr/m);
        GVtry=repelem(GV,1,nSampsPer_ts);
        t_samples=(1:numel(GVtry))/numel(GVtry)*tr;
        singleSamp=interp1(t_samples,GVtry,tSingleSamp);
        allSamps=repmat(singleSamp,1,nSamps);
        allSamps=allSamps(1:lent);
        sampSig=real(filtSG_tf(allSamps,t,f,round(200e9/df),10,1));
        
        %%%%%%%%%%%%%%%%
        %% Dispersion %%
        %%%%%%%%%%%%%%%%
        
        phi=p*m*(tr/m)^2/(2*pi);% this is equal to p*m*ts^2/(2*pi)
        
end


% Leave below uncommented for all cases.

phi2perKm=   2.1823e-23;
totalDispersion_equivalent_SMF_km=phi/phi2perKm
phaseGVD=phi/2*(2*pi*f).^2;



%%%%%%%%%%%%%%%%%%
%% Generate SUT %%
%%%%%%%%%%%%%%%%%%
SUT=ones(1,lent);



% % % % % % Linearly chirped signal
fini=0; ffin=fmax/2;
c=(ffin-fini)/tWind;
SUT=sin(2*pi*(c/2*t.^2+fini*t));


fSUT1=1/(ts*4);

% fSUT2=1/(ts*2);

% fSUT3=1/(ts);

% SUT=SUT+sin(2*pi*(fSUT1)*t);%+sin(2*pi*(fSUT2)*t)+sin(2*pi*(fSUT3)*t);

%%%%%%%%%%%%%%%%
%% Processing %%
%%%%%%%%%%%%%%%%


sampSUT=exp(1j*sampSig).*SUT   ;
sampSUT_f=nfft(sampSUT,dt,scale);

sampSUTdisp_f=(sampSUT_f).*exp(1j*phaseGVD);
sampSUTdisp=nifft(sampSUTdisp_f,Fs,scale);


%% Filter Spectrogram Signal
sampSUTdisp=(filtSG_tf(sampSUTdisp,t,f,round((10/tr)/df),5,0));

%%% Reshape output signal to get 2D representation
sptgm=reshape(circshift(sampSUTdisp(1:lensInds(end)),0),ntr,(numel(lensInds)-1));
f_spec=(tSingleSamp-max(tSingleSamp)/2)/phi/(2*pi);
t_spec=linspace(-tWind/2,tWind/2,numel(lensInds)-1);



%% Plotting section

ylims=([-0.1 1.1]);
ylimsNegPos=([-1.1 1.1]);
numtrWindLims=20;
tlims='auto';%(t(end)/2+[-numtrWindLims*tr numtrWindLims*tr])*1e9;
% tlims=[t(1) t(end)];

tps=1e12; tns=1e9; tus=1e6;
lbps='Time (ps)'; lbns='Time (ns)'; lbus='Time (us)';
fG=1e-9; lbG='Frequency (GHz)';
figure('Renderer', 'painters', 'Position', [50 50 1200 700])




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
legend('abs','angle','spectrogram')
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
legend('real','imag','spectrogram')
ylabel('Phase (rad)')

subplot(4,4,13:16)

imagesc(t_spec*tns,f_spec*fG,(abs(sptgm).^2));
hcb=colorbar();
xlim(tlims)

ylabel(hcb,'Power (dB)')
xlabel('Time (ns)')





function waveform=superGauss(C,t0,m,xs,center)
waveform=exp(-(1+1j*C)/2*((xs-center)/t0).^(2*m));
end




function fftout=nfft(sig,varargin)
%fft gives the swaped spectrum.
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
