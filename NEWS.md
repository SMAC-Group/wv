# wv 0.1.1

## Features

This R package provides a series of tools to compute and plot quantities related to wavelet transformation and wavelet variance (WV). More specifically, it provides the following features: 

- Calculation and plotting of coefficients for the Discrete Wavelet Transform (DWT) and Maximum Overlap Discrete Wavelet Transform (MODWT).
- Computation of WV and robust WV as well as of their respective confidence intervals.
- Visualization of the estimated WV with different options comparing classic and robust estimates as well as different times series.

## Updates 

Compared to the 0.1.0 version, we add the application of wavelet variance on IMU data. Specifically, we add the following new features: 

- New function wvar.imu() which allows the computation of wavelet variance based on IMU data;
- New function plot.imu_wvar() which allows to plot the wavelet variance computed based on IMU data;
- A few datasets of wavelet variance based on real IMU data as examples to test the new functions.



