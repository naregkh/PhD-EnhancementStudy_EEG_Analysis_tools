# PhD-EnhancementStudy_EEG_Analysis_tools
Here are the main tools that I've used for analyisng the EEG data. 

## Time-Frequency analysis
A lot of my EEG processing pipeline for time-frequency analysis is inspired by [Mike X Cohen](https://www.mikexcohen.com/)'s lecturelets and book ([Analyzing Neural Time Series Data: Theory and Practice](https://mitpress.mit.edu/books/analyzing-neural-time-series-data)).
Thank you Mike.

### TF_decomposition.m
Performs time-frequency decomposition of the signal using a complex Morlet Wavelet convolution.

### TF_decomposition_wrapper.m
Performs TF_decomposition.m for the different conditions. 

### Extract_power.m
Extracts power within the specified time-frequency window. 

### Extract_power_wrapper.m
uses extract_power.m, extracts power within the specified time-frequency window and plots the time-frequncy map, draws a windows around the time-window of the interest, and shows the power of specified frequency for each time point in a line plot. I'm also not sure why I called this "wrapper" it should clearly be called something else. 

## Event-related potential analysis

### ERP_stats_v3.m
Plots two ERP components with their conf intervals.
Uses cluster-size permutation test correcting for multiple comparison test.
The function takes two within participants ERP components and compares them visually also plotting their 95% confidence interval, plots the difference wave, performs cluster-size permutation test to show clusters (or time points which are significantly different between the two signals) which are significant.

