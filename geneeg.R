# ERP template used in EJN 2024 simulations
# true onset = 160 ms, F=81, max at F=126
# true_onset <- 160
# true_offset <- 340
# Xf <- seq(0, 500, 2)
# Nf <- length(Xf)
# temp1 <- vector(mode = "numeric", length = Nf)
# erp <- dnorm(seq(-1.5,1.5,length.out=93), 0, 1)
# erp <- erp - min(erp)
# erp <- erp / max(erp)
# temp2 <- c(rep(0, 79), erp, rep(0, 79))

# Adapted from Matlab code:
# https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator

# function signal = noise (frames, epochs, srate)
# 
# % function signal = noise (frames, epochs, srate)
# %
# % Function generates noise with the power spectrum of human EEG
# % Inputs:
#   %  frames - number of signal frames per each trial
# %  epochs - number of simulated trials
# %  srate - sampling rate of simulated signal
# % Output:
#   %  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
# % Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# 
# load meanpower  
# sumsig = 50;	%number of sinusoids from which each simulated signal is composed of
# 
# signal = zeros (1, epochs * frames);
# for trial = 1:epochs
# freq=0;
# range = [(trial-1)*frames+1:trial*frames];
# for i = 1:sumsig
# freq = freq + (4*rand(1));
# freqamp = meanpower(min (ceil(freq), 125)) / meanpower(1);
# phase = rand(1)*2*pi;
# signal (range) = signal (range) + sin ([1:frames]/srate*2*pi*freq + phase) * freqamp;
# end
# end

# Generate noise with the power spectrum of human EEG
# Inputs:
#   npoints - number of time points in each trial
#   srate - sampling rate of simulated signal
#   outvar - variance of the noise
#   nsig - number of sinusoids to sum
#   ntrials - number of trials
#   meanpower - EEG power spectrum estimated by Yeung et al.
# Output:
#   matrix of simulated noise: time points x trials
# Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# R version + variance control: GAR, University of Glasgow, Oct 2023
# Running the Matlab code many times, the median variance was about 0.4, so setting that value as the default.
eeg_noise <- function(npts = 251, srate = 500, outvar = 0.4, nsig = 50, ntrials = 50, meanpower){
  out <- matrix(0, nrow = npts, ncol = ntrials)
  for(T in 1:ntrials){
    freq <- 0
    for(i in 1:nsig){
      freq <- freq + (4*runif(1))
      freqamp <- meanpower[min(ceiling(freq), 125)] / meanpower[1]
      phase <- runif(1)*2*pi
      out[,T] <- out[,T] + sin((1:npts)/srate*2*pi*freq + phase)*freqamp
    }
    # control variance
    out[,T] <- out[,T] - median(out[,T])
    out[,T] <- out[,T] / sd(out[,T])
    out[,T] <- out[,T] * sqrt(outvar)
  } # trials
  out
}

# function signal = peak (frames, epochs, srate, peakfr, position, tjitter, wave)
# 
# % function signal = peak (frames, epochs, srate, peakfr {, position, tjitter, wave})
# %
# % Function generates signal composed of a single peak in each epoch
# % Required inputs:
#   %  frames - number of signal frames per each trial
# %  epochs - number of simulated trials
# %  srate - sampling rate of simulated signal
# %  peakfr - frequency of sinusoid whos half of the cycle is taken to form the peak
# % Optional inputs:
#   %  position - position of the peak [in frames]; default: frames/2 => in the middle
# %  tjitter - stdev of time jitter of the peak; default: 0 => no jitter
# %  wave - if defined the signal is composed not from a peak, but complete sinusoid. 
# % Output:
#   %  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
# % Implemnted by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# 
# if nargin < 5
# position = frames / 2;
# end
# if nargin < 6
# tjitter = 0;
# end
# if nargin < 7
# wave = 0;
# end
# 
# %generating peak
# signal = zeros (1, epochs * frames);
# for trial = 1:epochs
# pos = position + round(randn(1)*tjitter);
# for i=1:frames
# phase = (i-pos)/srate*2*pi*peakfr;
# if wave | (phase < pi/2 & phase > -pi/2)
# signal((trial-1)*frames+i) = cos(phase);
# end
# end
# end

# R version: GAR, University of Glasgow, Nov 2024
# INPUT:
# seql = sequence length in ms (default 501 ms)
# ntrials = number of trials (default 50)
# srate = sampling rate (default 500 Hz)
# peakfr = frequency of sinusoid (default 10 Hz)
# onset = start of the signal in ms (default 150 ms)
# onset_jitter = jitter of onset in SD (normal distribution, default 0)
# sinhalf = signal is half a sinusoid cycle, otherwise one (default TRUE)
# amp = peak amplitude
# OUTPUT:
# simulated EEG signal =  matrix npts*ntrials 
eeg_signal <- function(seql = 501, ntrials = 50, srate = 500, peakfr = 10, 
                       onset = 150, onset_jitter = 0, sinhalf = TRUE, amp = 1){
  Xf <- seq(from = 0, to = seql, by = 1000*1/srate)
  Nf <- length(Xf)
  # res <- matrix(0, nrow = ntrials, ncol = length(Xf))
  res.vec <- vector(mode = "numeric", length = Nf)
  trial_onset <- onset + round(rnorm(n=1, mean=0, sd=onset_jitter))
  sps <- 1/srate # seconds per sample 
  # full sinusoid:
  # end.epoch <- seql/1000 # seconds 
  # t.seq <- seq(from = 0, to = end.epoch, by = sps) # seconds 
  # res <- sin(2*pi*peakfr*t.seq)
  # plot(t.seq, res, type="l")
  # plot(Xf, res, type="l")
  if(sinhalf){ # get time period for half a cycle  
    time.period <- 1/peakfr/2
  } else { # get time period for one cycle
    time.period <- 1/peakfr
  }
  # time points
  t.seq <- seq(from = 0, to = time.period, by = sps)
  tmp <- sin(2 * pi * peakfr * t.seq) * amp
  # plot(t.seq, tmp, type="l")
  # onset location
  F <- which(Xf == onset)
  # shift signal by one time point because first value is zero
  # res[T,F:(F+length(tmp)-1)] <- tmp
  # res[T,(F-1):(F+length(tmp)-2)] <- tmp
  # shift signal by one time point because first value is zero
  res.vec[(F-1):(F+length(tmp)-2)] <- tmp
  # plot(Xf, res.vec, type="l")
  if(ntrials>1){ # if multiple trials
    res <- matrix(rep(res.vec, ntrials), nrow = Nf, ncol = ntrials)
  } else {
    res <- res.vec
  }
  res
}

# wrapper function that calls eeg_signal and eeg_noise to generate trials
eeg_sim <- function(seql = 501, ntrials = 50, srate = 500, peakfr = 10, 
                    onset = 150, onset_jitter = 0, sinhalf = TRUE, amp = 1,
                    outvar = 0.4, nsig = 50, ...){
  Xf <- seq(from = 0, to = seql, by = 1000*1/srate)
  Nf <- length(Xf)
  res.signal <- eeg_signal(seql = seql, ntrials = ntrials, srate = srate, peakfr = peakfr, 
                           onset = onset, onset_jitter = onset_jitter, sinhalf = sinhalf, amp = amp)
  res.noise <- eeg_noise(npts = Nf, srate = srate, outvar = outvar, 
                         nsig = nsig, ntrials = ntrials, meanpower)
  res <- res.signal + res.noise
  res
}

# Wrapper function that calls eeg_signal and eeg_noise to generate trials at multiple electrodes.
# Variable amp.w is a vector of amplitude weights for the signal -- each value corresponds to an electrode.
# Variable noise.cst determines if the same noise is used at all electrodes (TRUE) or different noise (FALSE).
# Returns a 3-D array electrodes x time points x trials.
# Individual trials are low-pass filtered at 30 Hz using 4-th order Butterworth filter.
eeg_sim_elec <- function(seql = 501, ntrials = 50, srate = 500, peakfr = 10, 
                         onset = 150, onset_jitter = 0, sinhalf = TRUE, amp.w = c(1,0.75,0.3,0),
                         outvar = 0.4, nsig = 50, noise.cst=FALSE, ...){

  # Define filter
  fs <- srate # sampling rate
  cutoff <- 30 # low-pass cut-off
  bf <- butter(n=4, w=cutoff / (fs / 2), type = "low") 

  Xf <- seq(from = 0, to = seql, by = 1000*1/srate)
  Nf <- length(Xf)
  Ne <- length(amp.w)
  res <- array(0, dim = c(Ne,Nf,ntrials))
  res.noise <- eeg_noise(npts = Nf, srate = srate, outvar = outvar, 
                         nsig = nsig, ntrials = ntrials, meanpower)
  for(E in 1:Ne){
    if(amp.w[E]==0){
      res.signal <- matrix(0, nrow = Nf, ncol = ntrials)    
    } else {
      res.signal <- eeg_signal(seql = seql, ntrials = ntrials, srate = srate, peakfr = peakfr, 
                               onset = onset, onset_jitter = onset_jitter, sinhalf = sinhalf, amp = amp.w[E])
    }
    if(noise.cst==FALSE && E>1){
      res.noise <- eeg_noise(npts = Nf, srate = srate, outvar = outvar, 
                             nsig = nsig, ntrials = ntrials, meanpower)
    }
    res[E,,] <- res.signal + res.noise
    for(T in 1:Nt){ # low-pass filter
      res[E,,T] <- filtfilt(bf, res[E,,T])
    }
  }
  res
}

