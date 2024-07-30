
% -----------------------------------------------------------------
%  signal_psd_1sided.m
%
%  This function computes the onesided power spectral density (PSD)
%  of a given signal x(t) in time domain.
%
%  input:
%  x    - (Ndt x 1) signal in time domain
%  fs   - sampling rate (Hz)
%  Nfft - FFT number of points (optional)
%  win  - (Ndt x 1) window (optional)
%
%  output:
%  PSD_X - one-sided power spectral density of x
%  freq  - one-sided frequency vector (Hz)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: May 31, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [PSD_X,freq] = signal_psd_1sided(x,fs,Nfft,win)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    elseif nargin == 2
    	Nfft = length(x);
        win  = ones(Nfft,1);
    elseif nargin == 3
        win = ones(length(x),1);
    end
    
    % check if fs is positive
    if fs <= 0.0
        error('fs must positive')
    end
    
    % check if Nfft is an integer greater than 1
    if Nfft <= 1
        error('Nfft must be an integer greater than 1')
    end
    
    % convert x to a column vector (if necessary)
    if find( size(x) == max(size(x)) ) > 1
        x = x';
    end
    
    
    % time step (s)
    dt = 1/fs;
    
    % compute FFT and frequency vector
    [Xfft,freq] = signal_fft_1sided(x,fs,Nfft,win);
	
    % compute PSD
    PSD_X = Xfft.*conj(Xfft);
    
    % normalize PSD
    PSD_X(2:end-1) = 2*dt*PSD_X(2:end-1);
	
return
% -----------------------------------------------------------------
