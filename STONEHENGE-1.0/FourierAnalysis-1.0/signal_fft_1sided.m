
% -----------------------------------------------------------------
%  signal_fft_1sided.m
%
%  This function computes the onesided FFT of a signal x(t) in
%  time domain and returns the corresponding transformned signal
%  X(freq) in frequency domain.
%
%  input:
%  x    - (Ndt x 1) signal in time domain
%  fs   - sampling rate (Hz)
%  Nfft - FFT number of points (optional)
%  win  - (Ndt x 1) window (optional)
%
%  output:
%  X    - (Nfft/2 x 1) transformed signal in frequency domain
%  freq - (Nfft/2 x 1) one-sided frequency vector (Hz)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: May 31, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [X,freq] = signal_fft_1sided(x,fs,Nfft,win)

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
    
    % check arguments
    if fs <= 0.0
        error('fs must positive')
    end
    
    if Nfft <= 1
        error('Nfft must be an integer greater than 1')
    end
    
    % convert x to a column vector (if necessary)
    if find( size(x) == max(size(x)) ) > 1
        x = x';
    end
    
    
    % sampling points (bins) column vector
    bins = linspace(0,Nfft-1,Nfft)';
    
    % frequency resolution (Hz)
    df = fs/Nfft;
    
	% frequency column vector (Hz)
    freq = df*bins;
    
	% compute windowed FFT
	X = fft(x.*win,Nfft);
    
	% normalize windowed FFT (to ensure Parseval identity)
	X = X/sqrt(Nfft);
    
    % one-sided frequency vector
    freq = freq(1:ceil(Nfft/2));
    
    % one-sided FFT of X
    X = X(1:ceil(Nfft/2));

return
% -----------------------------------------------------------------
