# 2PAnalyses
Code for analysing 2 photon imaging data
# Processing recordings in the Cylindrical Augmented Reality System (SCAR)
#### By: Jork de Kok 
#### Last update: 10/04/2025

## 0. Extract async pulse, run compression and spike sorting with kilosort  
For all these different steps use old script from Masa/Diao/Edd. 
These scripts are accessible here: (link to github repo)
<br> <br>

#### 0.1. main_LFP_preprocessing_and_saving.m (VR_NPX_analysis)
Run the first section (LFP -> Catgt) this extracts the async pulse and performs preprocessing 
steps on the LFP signal. 
<br> <br>

#### 0.2. Run compression and spike sorting on a Linux computer 
See https://github.com/SaleemLab/si_lab for documentation
<br> <br>

## 1. importing packages and loading in data 
Fill the mouse, session and so on and get all the quality metrics, ephys, bonsai, SLEAP and metadata.
<br> <br>
<br> <br>
## 2. extract and process async pulse 

#### 2.1. rename_cols_SLEAP (processing_async)

_Usage: rename_cols_SLEAP(SLEAP_df)_

_Output: none_
<br> <br>

#### 2.2. rename_cols_bonsai (processing_async)

_Usage: rename\_cols\_bonsai(bonsai\_df)_

_Output: none_
<br> <br>

#### 2.3) process_async (processing_async)

Converts timestamps to seconds, binairize the async signal and determines timepoints when async went on
This function will also plot the signal before binarization and after. The input must be a dataframe with 
a column called 'Timestamp'. The output of the function is a vector containing when the sync pulse went on 
in seconds

_Usage: x\_async\_on = process\_async(x\_df)_

Example: SLEAP\_async\_on = process\_async(SLEAP\_df) 

_Output: SLEAP\_async\_on_
<br> <br>

#### 2.4) restructure_ephys_async (processing_async)

Converts previously extracted async pulse to a dataframe 

_Usage: ephys\_async\_on = restructure\_ephys\_async(ephys\_async\_on)_

_output: ephys\_async\_on_
<br> <br>
<br> <br>

##  Section 3: preprocess SLEAP data 

#### 3.1) alignment_angle (circular_processing)

determines angle of rotation that is necessary to align with virtual environment 
based on if orthographic or cubemap view is used

_Usage: alingment\_angle = find\_aligment\_angle(session\_meta\_data)_

_Output: alingment\_angle_

<br> <br>

#### 3.2) process_trajectory_and_HD (circular_processing)

input are sleap data, the data containing location, alignment angle and a radius (r_lim). 
This function first plots the trajectory of the head of the animal. then applies radial 
restriction and plots the same trajectory after that. After this the trajectory is rotated
to the align with the virtual environment and the trajectory is plotted again. Lastly, 
the head direction is calculated and plotted. (contains subfunctions in the same module)

_Usage: process\_trajectory\_and\_HD(SLEAP\_df, columns\_containing\_trajectory\_data, alignment\_angle, r\_lim)_

_Output: none_
<br> <br>

##  Section 4: dataframes on ephys time 

#### 4.1) align_and_interpolate (align_and_interpolation_functions)

performs crosscorrelation on 2 time vectors (ephys + bonsoi or SLEAP) that contain the time the async pulse went on
then it determines the best lag and if one of the timevectors contains more time points, it will be 
cut to be the same length as the other timevector. After that linear interpolation is applied to the input 
dataframe based on the ephys time. 

_Usage: align_and_interpolate(x_async_on, ephys_async_on, x_df)_

Example: align_and_interpolate(SLEAP_async_on, ephys_async_on, SLEAP_df) 

_Output: none_
<br> <br>

#### 4.2) apply_interpolation_based_on_new_tvec (align_and_interpolation_functions)

creates a new time vector based on the start and end time of the SLEAP_df and bonsai_df and the sampling rate. 
After that columns in SLEAP_df containing trajectory data (contain X or Z) are linearly interpolated and other 
columns in SLEAP_df are interpolated based on the nearest value. Then for the bonsai_df the status_of_experiment
column is converted from categorical value to numerical followed by interpolation based on the nearest value. After 
that the status_of_experiment column is converted back to a categorical value. Lastly, spikes are binned based on 
the new tvec. (contains subfuctions in the same module) 
<br>

_Usage: SLEAP_df, bonsai_df, spikes_df = apply_interpolation_based_on_new_tvec(SLEAP_df,bonsai_df,spikes_raw, sampling_rate)_

_output: SLEAP_df, bonsai_df, spikes_df_
<br> <br>

##  Section 5: processing trajectory and HD 

#### 5.1) interpolate_and_smooth (align_and_interpolation_functions)

first plots trajectory of the implantX and Z then it perform interpolation over the columns containing trajectory data.
then it plots the trajectory again and then smoothing is applied followed by plotting the trajectory again. Lastly, 
the new head direction is calculated based on the interpolated and smoothed trajectory. By default the window size for 
the smoothing is set to 3. 

_Usage: interpolate_and_smooth(SLEAP_df, columns_containing_trajectory_data,window_size = 3)_

Example: SLEAP_df = interpolate_and_smooth(SLEAP_df, columns_containing_trajectory_data)

_output: SLEAP_df_
<br> <br>

#### 5.2) calc_distance_and_speed (binning_spikes)

calculated the distance and speed of trajectory and also for the head direction. 

_Usage: calc_distance_and_speed(SLEAP_df, sampling_rate)_

_output: none_
<br> <br>

#### 5.3) landmark_rotation (circular_processing) 

Calculated the global direction by equating the first baseline to 0 degrees. Then plots the trajectory of the implantX and Z
of SLEAP_df after which the trajectory is rotated relative to the global direction and then the trajectory is plotted again. 
After that a histogram of the HD is plotted followed by calculating the rotated head direction based on the landmark rotation.
Lastly, the rotated HD is plotted. Bins can adjust the amount of bins in the histogram (contains subfunction in the same module)

_Usage: landmark_rotation(SLEAP_df,global_direction, columns_containing_trajectory_data, bins =20)_

Example: SLEAP_df, global_direction  landmark_rotation(SLEAP_df,global_direction, columns_containing_trajectory_data)

_output: SLEAP_df, global_direction_
<br> <br>

##  Section 6: binning spikes based on HD 

#### 6.1) create_HD_bins (binning_spikes)

creates a new dataframe based which can be based on specific time range (list with start and end time). The new dataframe contains
columns that contian the witdh of the HD bins, the midle of bins in rad, middle of the bins in degrees, the occupancy of HD and the 
rotated (based on landmark) HD. Lastly, it plots the HD occupancy and the rotated HD occupancy. 

_Usage: create_HD_bins(SLEAP_df, sampling_rate, num_bins=20, time_range=False)_

Example: HD_bins = create_HD_bins(SLEAP_df, sampling_rate=1/30, num_bins=20, time_range=[0,900])

_output: HD_bins_
<br> <br>

#### 6.2) apply_spikes_per_HD_bin (binning_spikes)

calculculates spikes per bin which were defined HD_bins dataframe. The binning can be based on a specfic time range. (contains subfuctions)

_Usage: apply_spikes_per_HD_bin(HD_bins, SLEAP_df, spikes_df, sampling_rate, num_bins, time_range=False)_

Example: binned_spikes_HD = apply_spikes_per_HD_bin(HD_bins, SLEAP_df, spikes_df, sampling_rate=1/30,num_bins=20, time_range=[0,900])

_output: binned_spikes_HD_
<br> <br>

#### 6.3) calc_neuron_characteristics (circular_processing)

calculates max firing, absolute vector lenght, rayleigh vector lenght and preferred firing direction of all neurons. 

_Usage: calc_neuron_characteristics(binned_spikes_HD,neurons)_

Example: neuron_characteristics = calc_neuron_characteristics(binned_spikes_HD, neurons)

_output: neuron_characteristics_
<br> <br>

##  Section 7: quality metrics 

#### 7.1) quality_metrics_thresholds (metrics) 

applied threshold on quality metrics and outputs all neurons that adhere to the thresholds. 

_Usage: quality_metrics_thresholds(quality_metrics)_

Example: valid_neurons_ids = metrics.quality_metrics_thresholds(quality_metrics=quality_metrics)

_output: valid_neurons_ids_
<br> <br>

##  Section 8: plot binned HD spikes and heatmap 
