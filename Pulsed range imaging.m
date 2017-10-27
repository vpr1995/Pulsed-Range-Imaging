%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PULSED RANGE IMAGING        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
cj=sqrt(-1);
pi2=2*pi;
c=3e8;                        % Propagation speed
%B0=100e6;                % Baseband bandwidth is plus/minus B0
B0 = 80e6;
w0=pi2*B0;
fc=2.4e9;                       % Carrier frequency
wc=pi2*fc;
Xc=2;                     % Range distance to center of target area
Tp = 1e-6;                  % Range swath echo time period % Y: Xc-X0<=x<=Xc+X0
X0=1                    
Tx=X0*4/c;
alpha=w0/Tp;                 % Chirp rate % Y
beta=wc-alpha*Tp;         % % Y: start point of frequency beta = wcm = wc -w0.
                                       % carrier freq. of the chirp is the
                                       % mid-freq. wc = beta + alpha*Tp,
                                       % and the spectral support band is
                                       % [beta, beta+2alpha*Tp]
                                       % p(t) = a(t)*exp(j*beta*t+j*alpha*t^2).
                                       % wi(t) = beta*t+2*alpha*t, t>=0,
                                       % alpha>0.

dx=c/(4*B0);                 % Range resolution % Y
dt=pi/(2*alpha*Tp);      % Time domain sampling (guard band plus minus
                                     % 50 per % Y: ws>=4w0) or use dt=1/(2*B0) for a general
                                     % radar signal % dt = 1/(2*B0)=pi/w0. dt means
                                     % Delta_t
dtc=pi/(2*alpha*Tx);    % Time domain sampling for compressed signal
                                     % (guard band plus minus 50 per) % Y: dtc means
                                     % Delta_{tc}
Ts=(2*(Xc-X0))/c;          % Start time of sampling % Y
Tf=(2*(Xc+X0))/c+Tp;   % End time of sampling % Y

% If Tx < Tp, choose compressed signal parameters for measurement
flag=0;                  % flag=0 indicates that Tx > Tp
if Tx < Tp,
   flag=1;                 % flag=1 indicates that Tx < TP
   dt_temp=dt;             % Store dt % Y
   dt=dtc;                 % Choose dtc (dtc > dt) for data acquisition
end;

% Measurement parameters
n=2*ceil((.5*(Tf-Ts))/dt);        % Number of time samples % Y: n is even
t=Ts+(0:n-1)*dt;                    % Time array for data acquisition
dw=pi2/(n*dt);                       % Frequency domain sampling % Y: fs = 1/dt is the bandwidth
w=wc+dw*(-n/2:n/2-1);        % Frequency array (centered at carrier) 
x=Xc+.5*c*dt*(-n/2:n/2-1);   % range bins (array); reference signal is
                                               % for target at x=Xc. % Y: t=2x/c, so x =
                                               % tc/2.
kx=(2*w)/c;                            % Spatial (range) frequency array % Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntarget=2;                        % number of targets
%%%%%%%%%%%%% Targets' parameters  %%%%%%%%%%%%%%%%%%

% xn: range;              fn: reflectivity

xn(1)=-0.9*X0;                   fn(1)=1;
% xn(2)=0.7*X0;               fn(2)=0.8;
% xn(3)=0.8*X0;              fn(3)=1;
% xn(4)=-0.5*X0;              fn(4)=0.8;
xn(2)=1*X0;               fn(2)=0.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=zeros(1,n);              % Initialize echoed signal array

for i=1:ntarget;
     td=t-2*(Xc+xn(i))/c; % Y: td are the actual delayed time points, t includes the delay of the closest reflector
     pha=beta*td+alpha*(td.^2);         
     % Chirp (LFM) phase 
     % Y: beta = wc-w0 is the modified chirp carrier freq.
     s=s+fn(i)*exp(cj*pha).*(td >= 0 & td <= Tp); 
     % Y: td is the new time points, only when 0<=td<=Tp, p(td)~=0.
end

% If flag=1, i.e., Tx < Tp, perform upsampling
if     flag == 1
       td0=t-2*(Xc+0)/c;
       pha0=beta*td0+alpha*(td0.^2); % Reference chirp phase 
       % Baseband compressed signal                               
       scb = conj(s).*exp(cj*pha0).*exp(cj*2*beta*Xc/c-cj*4*alpha*Xc^2/(c^2));
       
       scb=[scb,scb(n:-1:1)];  % Append mirror image in time to reduce wrap 
                             % around errors in interpolation (upsampling)
       fscb=fty(scb);          % F.T. of compressed signal w.r.t. time
       dt=dt_temp;                     % Time sampling for echoed signal
       n_up=2*ceil((.5*(Tf-Ts))/dt);   % Number of time samples for upsampling
       nz=n_up-n;                      % number of zeros for upsmapling is 2*nz % Y
       fscb=(n_up/n)*[zeros(1,nz),fscb,zeros(1,nz)]; 

       scb=ifty(fscb);
       
        scb=scb(1:n_up);            % Remove mirror image in time % Y: now scb is the upsampled signal.
        scb1=scb;
       % Upsampled parameters
       n=n_up;
       t=Ts+(0:n-1)*dt;             % Time array for data acquisition
       dw=pi2/(n*dt);                % Frequency domain sampling
       w=wc+dw*(-n/2:n/2-1);        % Frequency array (centered at carrier)
       x=Xc+.5*c*dt*(-n/2:n/2-1);   % range bins (array); reference signal is
                                              % for target at x=Xc.
       kx=(2*w)/c;                    % Spatial (range) frequency array
       s=conj(scb).*exp(cj*beta*t+cj*alpha*t.^2-cj*4*alpha*Xc*t/c);
end

td0=t-2*(Xc+0)/c;
pha0=beta*td0+alpha*(td0.^2);         % Chirp (LFM) phase
s0=exp(cj*pha0).*(td0 >= 0 & td0 <= Tp);
% *******************************************************************************************************************************************
% Baseband conversion
sb=s.*exp(-cj*wc*t); 
figure(1)
plot(t,real(sb),'b'); grid on;
xlabel('t(s)');
ylabel('Real Part');
title('Baseband Echoed Signal s_b(t)');   % 1st question
% ********************************************************************************************************************************************
sb0=s0.*exp(-cj*wc*t);
figure(2)
plot(t,real(sb0),'b'); grid on;
xlabel('t(s)');
ylabel('Real Part');
title('Baseband reference Echoed Signal s_0b(t)'); % 2nd question
% *********************************************************************************************************************************************
% Fourier transform
Fsb=fty(sb);
figure(3)
plot((w-wc)/pi2, abs(Fsb), 'b'); grid on;
xlabel('Frequency');
ylabel('Magnitude');
title('Baseband Echoed Signal Spectrum'); % 3rd question

% ***********************************************************************************************************************************************
fsb0=fty(sb0);
figure(4)
plot((w-wc)/pi2, abs(fsb0), 'b'); grid on;
xlabel('Frequency');
ylabel('Magnitude');
title('Baseband Reference Echoed Signal Spectrum'); % 4th question
% ************************************************************************************************************************************************
Fsmb=Fsb.*conj(fsb0);
figure(5)
subplot(1,2,1);
plot((w-wc)/pi2,abs(Fsmb)); grid on;
xlabel('Frequency');
ylabel('Magnitude');
title('Baseband Matched filtered Signal Spectrum');
subplot(1,2,2);
plot(kx,abs(Fsmb),'b');
xlabel('Frequency');
ylabel('Magnitude');
title('Fourier transform of target function');% 5th question
% *************************************************************************************************************************************************
smb=ifty(Fsmb);
%  E= sum(abs(smb).*abs(smb));
figure(6)
subplot(1,2,1);
plot(x,abs(smb)./max(abs(smb)),'b'); grid on;
xlabel('Range');
ylabel('magnitude');
title('reconstruction via matched filter w.r.t Range');
axis([Xc-X0 Xc+X0 0 1.1]);
subplot(1,2,2);
plot(t,real(smb),'b');
xlabel('time');
ylabel('real part');
title('reconstruction via matched filter w.r.t Time');%6th question 
%***************************************************************************************************************************************************
figure(7)
plot(t,real(scb1),'b'); grid on;
xlabel('t(s)');
ylabel('Real Part');
title('time domain baseband compressed signal scb(t)');   % 7th question

% *****************************************************************************************************************************************************
fscb1=fty(scb1);
figure(8)
X=Xc+(dx*((w-wc)/pi2)*Tp);
plot(X,abs(fscb1)./max(abs(fscb1)),'b');
axis([Xc-X0 Xc+X0 0 1.1]);
grid on;
xlabel('Range');
ylabel('|scb(w)|');
title('Range reconstruction via time domain compression w.r.t range'); % 8th question
% ****************************************************************************************************************************************************
