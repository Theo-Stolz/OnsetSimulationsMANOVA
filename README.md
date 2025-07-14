# OnsetSimulationsMANOVA
Code and data for simulation Masters Thesis project on MEEG onset estimation. This project is in single-electrode and single-participant simulated data.

## R code
| Notebook | Content | Figures |
| ----- | ----- | ----- |
|`amplitude.Rmd`| Simulations of two signals that differ only in amplitude. | TBD |
|`latency.Rmd`| Simulations of two signals that differ only in latency. | TBD |
|`amplitudeandlatency.Rmd`| Simulations of two signals that differ in amplitude and latency. | TBD |
|`samesignal.Rmd`| Simulations of two signals do not differ. | TBD |
|`plots.Rmd`| Plots generated for report from simulations. | TBD |

### Measures in Simulations
| Measure | Calculation method | Object name |
| ----- | ----- | ----- |
| Change Point Onset | Onsets as determined by the change point detection algorithm from the @changepoint package. | cponsetT / cponsetF |
| Onsets of MAX thresholded statistics | Onsets as determined through MAX thresholded statistics. | onsetT / onsetF |
| Signal-to-Noise Ratio | Dividing median the T^2 or F statistic from baseline until end of signal by the mean absolute deviation along the baseline (from 0 to the onset).| signoiseT / signoiseF |
| Proportion of positive tests over time | Retrieving the number of tests that exceed the MAX threshold across simulations for each time point. | sigMaxT / sigMaxF |
| Proportion of positive tests over time | Retrieving the number of tests that exceed the 95th quantile threshold across simulations for each time point. | sigUniT / sigUniF |
| Number of significant time points retrieved | Sum of positive tests in a time-series by trial number. | TBD |

## Dependencies
Extra code dependencies are available in the code folder.file require its use.

## Simulation results
All simulation results are in the data folder, such that figures can be run without running the simulations.
