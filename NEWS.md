# wv 0.1.0

## Features

This R package provides a series of tools to compute and plot quantities related to wavelet transformation and wavelet variance (WV). More specifically, it provides the following features: 

- Calculation and plotting of coefficients for the Discrete Wavelet Transform (DWT) and Maximum Overlap Discrete Wavelet Transform (MODWT).
- Computation of WV and robust WV as well as of their respective confidence intervals.
- Visualization of the estimated WV with different options comparing classic and robust estimates as well as different times series.

## Updates

In this submission version, we made the following changes:

- We add more details in the description field of what the package does. We also add more reference that describes the methods in our pacakge. 
- For the function compare_wvar_split(), we add the use of on.exit() to ensure the settings are reset after we change the user's options.

