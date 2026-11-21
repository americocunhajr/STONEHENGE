
% -----------------------------------------------------------------
%  signal_psd.m
%
%  This function computes the power spectral density (PSD)
%  of a given temporal signal X(t). An algorithm based on the
%  FFT with a (natural) rectangular window is used.
%  
%  Reference:
%  C. Soize
%  Fundamentals of Random Signal Analysis
%  Application to Modal Identification in Structural Dynamics
%  Lectures Notes, 1997, pages 121--122
%
%  input:
%  X         - (Ndt x 1) temporal signal
%  freq_max  - maximun frequency in the band (Hz)
%  freq_samp - sampling rate (Hz)
%  Nfft      - number of points for FFT
%  Nsamp     - number of time steps in a signal copy
%  Ncopies   - number of independent copies of the signal
%
%  output:
%  psd_X - (Nfft x 1) power spectral density of X
%  freq  - (Nfft x 1) frequency vector (Hz)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: May 31, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [psd_X,freq] = signal_psd(X,freq_max,freq_samp,...
                                              Nfft,Nsamp,Ncopies)
    
    % check number of arguments
    if nargin < 6
        error('Too few inputs.')
    elseif nargin > 6
        error('Too many inputs.')
    end
    
    % check arguments
    if freq_max <= 0.0
        error('freq_max must positive')
    end
    
    if freq_samp <= 0.0
        error('freq_samp must positive')
    end
    
    if Nfft <= 0
        error('Nfft must positive integer')
    end
    
    if Nsamp <= 0
        error('Nsamp must positive integer')
    end
    
    if Ncopies <= 0
        error('Ncopies must positive integer')
    end
    
    % convert X to a column vector (if necessary)
    if find( size(X) == max(size(X)) ) > 1
        X = X';
    end
    
    % compute the size of X
    Ndt = length(X);
    
    if Ndt ~= Nsamp*Ncopies
        error('Ndt must be equal to Nsamp*Ncopies')
    end
    
    
    % frequency resolution or sampling frequency step (Hz)
    dfreq = freq_samp/Nfft;
    
    % time step (s)
    dt = 1/freq_samp;
    
    % period of acquisition for a signal copy (s)
    Tsamp = Nsamp*dt;
    
    % sampling points (bins) column vector
    bins_time = (0:1:(Nsamp-1))';
    bins_freq = (0:1:(Nfft -1))';
    
    % frequency vector (Hz)
    freq = -freq_max + dfreq*bins_freq;
    
    % convert the signal to a (Nsamples x Ncopies) table
    X_tab = reshape(X,Nsamp,Ncopies);

    % make the mean of X equal to zero
    X_zero_mean = X_tab - repmat(mean(X_tab,2),1,Ncopies);
    
    % define imaginary unit
    i = sqrt(-1);
    
    % complex exponentials vector
    exp_i_pi_m = exp(i*pi*bins_time);
    
    % rectangular windows coeficient
	rect_win_coef = 1/sqrt(Tsamp);
    
    % temporal signal windowed by the (natural) rectangular window
    X_win = rect_win_coef*dt*X_zero_mean.*repmat(exp_i_pi_m,1,Ncopies);

    % compute FFT of X_tab
    X_fft = fft(X_win,Nfft);

    % estimate the PSD of X
    psd_X = mean(abs(X_fft).^2,2)/(2*pi);
    
return
% -----------------------------------------------------------------
