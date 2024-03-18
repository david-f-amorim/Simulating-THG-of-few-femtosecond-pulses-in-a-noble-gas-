function [nu,four,amplitude,phase]=DFT(Nf,time,y_t,numin,numax)
% DFT  Calculates the discrete fourier transform of a dataset.
%
% [nu,four,amplitude,phase]=DFT(Nf,time,y_t,numin,numax)
% [nu,four,amplitude,phase]=DFT(0,time,y_t,0,0);
% [nu,four,amplitude,phase]=DFT(0,time,y_t,numin,numax);
% [nu,four,~,~]=DFT(0,time,y_t,numin,numax);
% [nu,~,amplitude,phase]=DFT(0,time,y_t,numin,numax);
%
% INPUTS
% - Nf        --> # of output frequency points. If set to 0 it is
%                 calculated from the time axis.
% - time      --> Time axis
% - y_t       --> Signal in the time domain
% - numin     --> Starting frequency
% - numax     --> End frequency. If set to 0 it is 
%                 calculated from the time axis
% 
% OUTPUTS
% - nu        --> Frequency axis
% - four      --> Fourier transform in complex a+ib notation
% - amplitude --> spectral amplitude
% - phase     --> spectral phase


N       = length(time);             % Number of frequency points;
dt      = abs(time(2)-time(1));     % Time resolution;
nu      = linspace(0,(1/dt)/2,N);   % Creates full frequency vector;
dnu     = abs(nu(2)-nu(1));         % Frequency resolution;

if numax == 0
    numax=1/dt/2;
end

if Nf == 0
Nf      = round((numax-numin)/dnu); % Frequency points in the region of...
end                                 % interest

clear nu dt dnu;

Dt      = diff(time);               % Sampling step in the time domain
Dt(N)   = Dt(N-1);

nu=linspace(numin,numax,Nf);        % Define the spectral window

four=zeros(1,Nf);

for ii=1:Nf
    
    four(ii)=sum(Dt.*y_t.*exp(-1i*2*pi*time*nu(ii)));
    
end; 

[phase,amplitude] = cart2pol(real(four),imag(four));
phase=unwrap(phase);


