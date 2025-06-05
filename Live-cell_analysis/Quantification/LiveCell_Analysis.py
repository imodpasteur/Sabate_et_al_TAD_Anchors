# -*- coding: utf-8 -*-
"""
@author: Thomas Sabaté
For usage of the code, please cite: Sabaté et al, Uniform dynamics of cohesin-mediated loop extrusion, Nature Genetics, 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import itertools
from matplotlib.patches import Rectangle
from matplotlib import ticker
from random import sample
from Utils_Live import *


def get_closing_rates_based_on_criterion(crit, threshold_t_extr):
    """
    Parameters
    ----------
    crit : dict
        Contains BIC values, parameters and closing rates estimated from the
        1- or 3-parameter models.

    Returns
    -------
    closing_rate : float
        Closing rate of the model associated with the lowest BIC value.
    fit_used : int
        number of parameter of the model associated with the lowest BIC value.
    """
    if hasattr(crit, '__len__'):     # check if crit is not nan
        t_extr = crit['3_parameter'][1][0]
        if t_extr >= threshold_t_extr:     # Constant distance
            closing_rate = np.nan
        if t_extr < threshold_t_extr:     # Linear decrease
            closing_rate = crit['3_parameter'][2]
    else:               # Not sufficient data to do the fitting
        closing_rate = np.nan
    return closing_rate


def get_squared_distances_minus_localization_errors_from_distances_and_scores(distances, scores, mean_Std_prec):
    """
    Parameters
    ----------
    distances : Array
        Array of distances aligned on proximal states. Last row corresponds to the timepoints closest to the proximal state.
        Columns correspond to different proximal states.
    scores : Array
        Scores (based on localization precision) corresponding to distances in 'distances'.
    mean_Std_prec : array of float
        Mean and standard deviation of the Gaussian used to convert scores to localization precision.

    Returns
    -------
    norm_mean_weighted_dist : 1D array
        Mean squared distances, weighted by localization precision.
    std_weighted : 1D array
        Mean standard error of the mean, weighted by localization precision.
    """
    sc_p = Score2prec(mean_Std_prec[0], mean_Std_prec[1])
    loc_prec = sc_p.GaussianRepartionFunctionShifted2Gaussian(scores)
    squared_loc_prec = np.power(loc_prec, 2)

    # Remove localization precision contribution
    closed_states_squared_raw = np.power(distances, 2)
    closed_states_squared = np.subtract(closed_states_squared_raw, squared_loc_prec)
    # replace infinite values (score of 0)
    closed_states_squared[np.isinf(closed_states_squared)] = closed_states_squared_raw[np.isinf(closed_states_squared)]
    
    # Compute mean weighted squared distance
    weighted_dist = np.multiply(closed_states_squared, scores)
    norm_mean_weighted_dist = np.divide(np.sum(weighted_dist, axis=1), np.nansum(scores, axis=1))
    
    # Compute weighted standard error of the mean
    sh = np.shape(closed_states_squared)
    nbtrack = sh[1]
    sem_weighted = np.divide(np.sqrt(np.divide(np.sum(np.multiply(np.power(np.transpose(np.subtract(np.transpose(distances), norm_mean_weighted_dist)), 2), scores),axis=1), np.sum(scores,axis=1))), np.sqrt(nbtrack))
    return norm_mean_weighted_dist, sem_weighted


# %% Load QC data
# Here, indicate the path to the .csv filtered distance time series files, which can be downloaded on Zenodo.

Cell_line = 'L1'

if Cell_line == 'L1':
    tracks_WT = pd.read_csv('/Path_To_Filtered_Dataset/L1_NoAuxin_QC.csv', sep=';')
    tracks_deplete = pd.read_csv('/Path_To_Filtered_Dataset/L1_Auxin2h_QC.csv', sep=';')
    size_loop = 345  # in kb
    
if Cell_line == 'L2':
    tracks_WT = pd.read_csv('/Path_To_Filtered_Dataset/L2_NoAuxin_QC.csv', sep=';')
    tracks_deplete = pd.read_csv('/Path_To_Filtered_Dataset/L2_Auxin2h_QC.csv', sep=';')
    size_loop = 566  # in kb

if Cell_line == 'T1':
    tracks_WT = pd.read_csv('/Path_To_Filtered_Dataset/T1_NoAuxin_QC.csv', sep=';')
    tracks_deplete = pd.read_csv('/Path_To_Filtered_Dataset/T1_Auxin2h_QC.csv', sep=';')
    size_loop = 918  # in kb

if Cell_line == 'HalfTAD':
    tracks_WT = pd.read_csv('/Path_To_Filtered_Dataset/NoTAD_NoAuxin_QC.csv', sep=';')
    tracks_deplete = pd.read_csv('/Path_To_Filtered_Dataset/NoTAD_Auxin2h_QC.csv', sep=';')
    size_loop = 576  # in kb
    
tracks_adjacent = pd.read_csv('/Path_To_Filtered_Dataset/Adjacent_QC.csv', sep=';')


# %% Segmentation of proximal states

# Add NAN in missing frames and interpolate distances
tracks_na = Add_NA_in_missing_frames_And_interpolate(tracks_WT.copy(), True)
tracks_deplete_na = Add_NA_in_missing_frames_And_interpolate(tracks_deplete.copy(), True)
tracks_adjacent_na = Add_NA_in_missing_frames_And_interpolate(tracks_adjacent.copy(), True)

if Cell_line == 'L1':
    thresh_space = 0.19945400124787457
if Cell_line == 'L2':
    thresh_space = 0.218277245599599
if Cell_line == 'T1':
    thresh_space = 0.2363151398746736
if Cell_line == 'HalfTAD':
    thresh_space = 0.2204404763636802
thresh_space_adjacent = 0.24938339452179015   # In um
time_threshold = 6     # In frames


# Segment proximal states
data = tracks_na.copy()
masked_WT_tracks, start_closed_state_notfiltered, end_closed_state_notfiltered = get_closed_state_in_tracks(data.copy(), thresh_space, time_threshold)
filtered_WT, start_closed_state, end_closed_state, censored = mean_filter_masked_tracks(masked_WT_tracks.copy(), time_threshold)

data = tracks_deplete_na.copy()
masked_deplete_tracks, start_closed_state_deplete_notfiltered, end_closed_state_deplete_notfiltered = get_closed_state_in_tracks(data.copy(), thresh_space, time_threshold)
filtered_deplete, start_closed_state_deplete, end_closed_state_deplete, censored_deplete = mean_filter_masked_tracks(masked_deplete_tracks.copy(), time_threshold)

data = tracks_adjacent_na.copy()
masked_adjacent_tracks, start_closed_state_adjacent_notfiltered, end_closed_state_adjacent_notfiltered = get_closed_state_in_tracks(data.copy(), thresh_space_adjacent, time_threshold)
filtered_adjacent, start_closed_state_adjacent, end_closed_state_adjacent, censored_adjacent = mean_filter_masked_tracks(masked_adjacent_tracks.copy(), time_threshold)

# %% Fraction

print('Proximal fraction No Auxin = ', list(filtered_WT['Proximal_filt']).count(1) / len(filtered_WT))
print('Proximal fraction + Auxin = ', list(filtered_deplete['Proximal_filt']).count(1) / len(filtered_deplete))
print('Proximal fraction adjacent = ', list(filtered_adjacent['Proximal_filt']).count(1) / len(filtered_adjacent))

# %% Frequency

s_frame = 30    # 30s/frame
mean_frequency = get_closed_state_frequency(filtered_WT.copy(), s_frame)
mean_frequency_deplete = get_closed_state_frequency(filtered_deplete.copy(), s_frame)
mean_frequency_adjacent = get_closed_state_frequency(filtered_adjacent.copy(), s_frame)

print('Proximal fraction No Auxin = %.2f h⁻¹'% (mean_frequency * 3600))
print('Proximal fraction + Auxin = %.2f h⁻¹'% (mean_frequency_deplete * 3600))
print('Proximal fraction adjacent = %.2f h⁻¹'% (mean_frequency_adjacent * 3600))

# %% Lifetimes

closed_lifetime, closed_lifetime_list, Censored, Censored_list = get_lifetime_closed_states(censored, start_closed_state, end_closed_state)
closed_lifetime_deplete, closed_lifetime_list_deplete, Censored_deplete, Censored_list_deplete = get_lifetime_closed_states(censored_deplete, start_closed_state_deplete, end_closed_state_deplete)
closed_lifetime_adjacent, closed_lifetime_list_adjacent, Censored_adjacent, Censored_list_adjacent = get_lifetime_closed_states(censored_adjacent, start_closed_state_adjacent, end_closed_state_adjacent)

conversion_to_min = 60 / s_frame
closed_lifetime_list_in_min = [i / conversion_to_min for i in closed_lifetime_list]
closed_lifetime_list_in_min_deplete = [i / conversion_to_min for i in closed_lifetime_list_deplete]
closed_lifetime_list_in_min_adjacent = [i / conversion_to_min for i in closed_lifetime_list_adjacent]

# %%% MLE fitting of lifetimes (as in Gabriele et al, Science, 2022)

mean, low, high = MLE_censored_exponential(closed_lifetime_list_in_min, Censored_list, 0.95)
print('No Auxin mean lifetime = ' + str(round(mean, 2)) + ' min')

mean_deplete, low_deplete, high_deplete = MLE_censored_exponential(closed_lifetime_list_in_min_deplete, Censored_list_deplete, 0.95)
print('+ Auxin mean lifetime = ' + str(round(mean_deplete, 2)) + ' min')

mean_adjacent, low_adjacent, high_adjacent = MLE_censored_exponential(closed_lifetime_list_in_min_adjacent, Censored_list_adjacent, 0.95)
print('Adjacent mean lifetime = ' + str(round(mean_adjacent, 2)) + ' min')


plt.figure()
n, bins, patches = plt.hist(closed_lifetime_list_in_min_adjacent, color='grey', edgecolor='k', bins='auto',
                            density=True, label='Adjacent', alpha=0.4)
y_fitted_adjacent_not_normalized = [(1/mean_adjacent) * np.exp(-i/mean_adjacent) for i in bins]
plt.semilogy(bins, [i * np.sum(np.diff(bins) * n) for i in y_fitted_adjacent_not_normalized], color='k', lw=5)

n, bins, patches = plt.hist(closed_lifetime_list_in_min, color='red', edgecolor='k', bins='auto',
                            density=True, label='- Auxin', alpha=0.4)
y_fitted_not_normalized = [(1/mean) * np.exp(-i/mean) for i in bins]
plt.semilogy(bins, [i * np.sum(np.diff(bins) * n) for i in y_fitted_not_normalized], color='darkred', lw=5)

n, bins, patches = plt.hist(closed_lifetime_list_in_min_deplete, color='b', edgecolor='k', bins='auto',
                            density=True, label='+ Auxin', alpha=0.4)
y_fitted_deplete_not_normalized = [(1/mean_deplete) * np.exp(-i/mean_deplete) for i in bins]
plt.semilogy(bins, [i * np.sum(np.diff(bins) * n) for i in y_fitted_deplete_not_normalized], color='darkblue', lw=5)

plt.xlabel('Proximal state lifetime\n(in min)', fontsize=21)
plt.ylabel('Counts', fontsize=22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.locator_params(axis='x', nbins=6)
plt.ylim((0.0007, 0.6))
plt.xlim((0, 80))
plt.legend(fontsize=22, frameon=False)

# %% Fit closing rate
# Use a more conservative segmentation of proximal states for closing rate fitting
if Cell_line == 'L1':
    thresh_space_CR = 0.13631379864169063
if Cell_line == 'L2':
    thresh_space_CR = 0.14625971061135684
if Cell_line == 'T1':
    thresh_space_CR = 0.1589465752706419
if Cell_line == 'HalfTAD':
    thresh_space_CR = 0.14335860065048628

# Conservative segmentation of proximal states
data = tracks_na.copy()
masked_WT_tracks, start_closed_state_notfiltered, end_closed_state_notfiltered = get_closed_state_in_tracks(data.copy(), thresh_space_CR, time_threshold)
filtered_WT_closingrate, start_closed_state_closingrate, end_closed_state_closingrate, censored_closingrate = mean_filter_masked_tracks(masked_WT_tracks.copy(), time_threshold)

# Retrieve proximal states in untreated time series and align distances on proximal states
window = 50  # size of fitting window (in frames)
closed_tracks, Scores = retrieve_closed_states_without_closed_states_before_wdw_timepoints_fulltracksOnly_WithScores(filtered_WT_closingrate, start_closed_state_closingrate, window)
Closed_arr, Score_arr = make_array_of_distances_scores_from_dict(closed_tracks, Scores, window)

# Randomly shuffle data
tracks_rdm = randomize_each_timepoint_in_df(filtered_WT_closingrate.copy())
# Detect proximal states in randomly shuffled data
masked_rdm_tracks, start_closed_state_notfiltered_rdm, end_closed_state_notfiltered_rdm = get_closed_state_in_tracks(tracks_rdm.copy(), thresh_space_CR, time_threshold)
filtered_rdm, start_closed_state_rdm, end_closed_state_rdm, censored_rdm = mean_filter_masked_tracks(masked_rdm_tracks.copy(), time_threshold)
# Retrieve proximal states in randomly shuffled data and align distances on proximal states
closed_tracks_rdm, Scores_rdm = retrieve_closed_states_without_closed_states_before_wdw_timepoints_fulltracksOnly_WithScores(filtered_rdm, start_closed_state_rdm, window)
Closed_arr_rdm, Score_arr_rdm = make_array_of_distances_scores_from_dict(closed_tracks_rdm, Scores_rdm, window)

# Fit closing rate
DEB = 0
window_fit = window
Threshold_t_extr = -121    # is seconds, threshold on t_extr
frame_s = 1 / s_frame    # frames / second
list_dist = list(filtered_deplete['Distance'])
list_precision = list(filtered_deplete['precision_Distance'])
R02 = np.nanmean([list_dist[i]**2 - list_precision[i]**2 for i in range(len(list_dist))])
Mean_std_prec = np.loadtxt('/Path_To_MeanAndStd_All_Cell_Lines/240626_Mean_Std_prec_All_Cell_Lines.txt')

closing_rate_3param, R2_estimate, slope_avg, x, t_extr, R2int_fit, criterion = fit_ClosingRate_TakingIntoAccount_LocalizationPrecision(Closed_arr, Score_arr, frame_s, size_loop, R02, Mean_std_prec)
closing_rate = get_closing_rates_based_on_criterion(criterion, Threshold_t_extr)
if math.isnan(closing_rate) is False:
    print('No auxin Closing rate = %.3f' % (closing_rate) + ' kb/s')
if math.isnan(closing_rate) is True:
    print('No auxin time series exhibit constant distances')

closing_rate_rdm_3param, R2_estimate_rdm, slope_avg_rdm, x_rdm, t_extr_rdm, R2int_fit_rdm, criterion_rdm = fit_ClosingRate_TakingIntoAccount_LocalizationPrecision(Closed_arr_rdm, Score_arr_rdm, frame_s, size_loop, R02, Mean_std_prec)
closing_rate_rdm = get_closing_rates_based_on_criterion(criterion_rdm, Threshold_t_extr)
if math.isnan(closing_rate_rdm) is False:
    print('Randomly shuffled closing rate = %.3f' % (closing_rate_rdm) + ' kb/s')
if math.isnan(closing_rate_rdm) is True:
    print('Randomly shuffled time series exhibit constant distances')

# %%% Plot closing rate fit in untreated and randomly shuffled time series

s_f = SpeedFitting(t_extr, R2int_fit, R2_estimate)
fitted_extr = s_f.func2slopes(x, [t_extr, R2int_fit, R2_estimate])

s_f_rdm = SpeedFitting(t_extr_rdm, R2int_fit_rdm, R2_estimate_rdm)
fitted_extr_rdm = s_f_rdm.func2slopes(x_rdm, [t_extr_rdm, R2int_fit_rdm, R2_estimate_rdm])


# Plots
close_square = Closed_arr ** 2
close_square_rdm = Closed_arr_rdm ** 2
normalized_close_square_mean, sem = get_squared_distances_minus_localization_errors_from_distances_and_scores(Closed_arr, Score_arr, Mean_std_prec)
normalized_close_square_mean_rdm, sem_rdm = get_squared_distances_minus_localization_errors_from_distances_and_scores(Closed_arr_rdm, Score_arr_rdm, Mean_std_prec)

x_plot = np.arange(- (window_fit + 1) * (1 / frame_s), 0 + 1/frame_s, 1/frame_s) / 60


plt.figure()
plt.plot(x_plot, normalized_close_square_mean_rdm, color='olive', lw=4)
plt.fill_between(x_plot, normalized_close_square_mean_rdm - sem_rdm, normalized_close_square_mean_rdm + sem_rdm, facecolor='olive', alpha=0.4, lw=0)
if math.isnan(closing_rate_rdm) is False:
    plt.plot(x_plot[DEB:window_fit], fitted_extr_rdm[DEB:window_fit], color='olive', ls='--', label='Shuffled: %.2f' % (closing_rate_rdm) + ' kb/s', lw=4)
if math.isnan(closing_rate_rdm) is True:
    plt.plot(x_plot[DEB:window_fit], fitted_extr_rdm[DEB:window_fit], color='olive', ls='--', label='Shuffled: Constant', lw=4)

plt.plot(x_plot, normalized_close_square_mean, color='r', lw=4)
plt.fill_between(x_plot, normalized_close_square_mean - sem, normalized_close_square_mean + sem, facecolor='r', alpha=0.4, lw=0)
plt.plot(x_plot[DEB:window_fit], fitted_extr[DEB:window_fit], color='darkred', ls='--', label='No auxin: %.2f' % (closing_rate) + ' kb/s', lw=4)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r't - t$_{\rm start}$ (in min)', fontsize=24)
plt.ylabel('Weighted squared\n3D distance (in $\mu m^2$)', fontsize=23)
plt.legend(fontsize=17, frameon=True)
plt.locator_params(axis='y', nbins=4)
plt.locator_params(axis='x', nbins=6)
plt.xlim(right=0)
plt.xlim((-window/2 - 0.5, 0))
plt.ylim((0, 0.2))


# %% Plot time series
# %%% Plot individual segmented tracks

data = filtered_WT.copy()
List_track_pair = list(data['Track_pair'].unique())
fps = 1/30   # frame / s
nbr_tracks_to_plot = 5
track_idx_to_plot = random.sample(List_track_pair, k=nbr_tracks_to_plot)
for pair in track_idx_to_plot:
    plt.figure(figsize=(8, 2))
    
    plt.plot(data['Frame'][data['Track_pair']==pair] / fps / 60, data['Distance'][data['Track_pair']==pair], color='r', lw=3)

    # Find positions of start and end of proximal state
    s = list(data['Proximal_filt'][data['Track_pair'] == pair])
    positions = []
    current_position = 0
    for value, group in itertools.groupby(s):
        group_length = len(list(group))
        if value == 1:
            positions.extend([current_position, current_position + group_length - 1])
        current_position += group_length
    deb_close_idx = positions[0::2]
    end_close_idx = positions[1::2]
    sub_df = data[data['Track_pair'] == pair]
    sub_df.reset_index(drop=True, inplace=True)
    deb_close_idx = [i + sub_df.loc[0, 'Frame'] for i in deb_close_idx]
    end_close_idx = [i + sub_df.loc[0, 'Frame']  for i in end_close_idx]
    
    # Plot segmented proximal states
    c = 0
    for i in deb_close_idx:
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((deb_close_idx[c] / fps / 60, -0.15), end_close_idx[c] / fps / 60 - deb_close_idx[c] / fps / 60, 0.15, facecolor="indigo"))
        c += 1

    plt.ylim((-0.1, 1.3))
    plt.xlabel('Time (in min)', fontsize=30)
    plt.ylabel('3D distance\n(in $\mu m$)', fontsize=25)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.xlim((0, 120))
    plt.locator_params(axis='y', nbins=4)
    plt.locator_params(axis='x', nbins=7)
    plt.gca().tick_params( length=6, width=1.5)
    plt.show()


# %%% Plot randomly chosen tracks together

to_plot = filtered_WT
fps = 1/30
to_plot['Frame_in_s'] = (to_plot['Frame'] * 0.5)
all_pairs = list(to_plot['Track_pair'].unique())
nbr_sample = 5
random_choice = sample(all_pairs, nbr_sample)

fig, axs = plt.subplots(nbr_sample, 1, figsize=(35, 35), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.15, wspace=0.2)
i = 0
for pair in random_choice:
    sub_df = to_plot[to_plot['Track_pair'] == pair]
    axs[i].plot(sub_df['Frame_in_s'], sub_df['Distance'], lw=8, color='r')
    
    # Find positions of start and end of proximal states
    s = list(to_plot['Proximal_filt'][to_plot['Track_pair'] == pair])
    positions = []
    current_position = 0
    for value, group in itertools.groupby(s):
        group_length = len(list(group))
        if value == 1:
            positions.extend([current_position, current_position + group_length - 1])
        current_position += group_length
    deb_close_idx = positions[0::2]
    end_close_idx = positions[1::2]

    sub_df.reset_index(drop=True, inplace=True)
    deb_close_idx = [i + sub_df.loc[0, 'Frame'] for i in deb_close_idx]
    end_close_idx = [i + sub_df.loc[0, 'Frame']  for i in end_close_idx]
    
    # Plot proximal states
    c = 0
    for p in deb_close_idx:
        axs[i].add_patch(Rectangle((deb_close_idx[c] / fps / 60, -0.15), end_close_idx[c] / fps / 60 - deb_close_idx[c] / fps / 60, 0.15, facecolor="indigo"))
        c += 1

    axs[i].set_xlim([0, 120])
    axs[i].set_ylim([-0.15, 1.3])
    if i != nbr_sample - 1:
        plt.setp(axs[i].get_xticklabels(), visible=False)
    
    axs[i].tick_params(labelsize=80, length=10, width=2.5)
    yticks = ticker.MaxNLocator(3)
    axs[i].yaxis.set_major_locator(yticks)
    i += 1

plt.xlabel('Time (in min)', fontsize=100)
fig.text(0.03, 0.5, '3D distance (in $\mu m$)', va='center', rotation='vertical', fontsize=100)
plt.show()

# %% Estimate fractions of each loop state

from utils_coordinate_difference import *
from models_coordinate_difference import *


def compute_var(sigma_x_green, sigma_x_red, sigma_y_green, sigma_y_red, sigma_z_green, sigma_z_red):
    """
    Parameters
    ----------
    sigma_x_green : float
        Mean localization precision in X in the green channel.
    sigma_x_red : float
        Mean localization precision in X in the red channel.
    sigma_y_green : float
        Mean localization precision in Y in the green channel.
    sigma_y_red : float
        Mean localization precision in Y in the red channel.
    sigma_z_green : float
        Mean localization precision in Z in the green channel.
    sigma_z_red : float
        Mean localization precision in Z in the red channel.

    Returns
    -------
    var_x : float
        Variance of localization precision in X for both channels.
    var_y : float
        Variance of localization precision in Y for both channels.
    var_z : float
        Variance of localization precision in Z for both channels.

    """
    var_x = sigma_x_green ** 2 + sigma_x_red ** 2
    var_y = sigma_y_green ** 2 + sigma_y_red ** 2
    var_z = sigma_z_green ** 2 + sigma_z_red ** 2
    return var_x, var_y, var_z


def split_dataset_per_track_WholeData(data):
    """
    Parameters
    ----------
    data : Pandas df
        Contains coordinates of anchors.

    Returns
    -------
    df1 : Pandas df
        Bootstrapped data, same size than 'data', with replacement
    df2 : Pandas df
        Bootstrapped data, same size than 'data', with replacement
        
    Bootstrapping is done on tracks, not timepoints
    """
    data_to_sample = data.copy()
    data_to_sample.reset_index(drop=True, inplace=True)
    list_pair = list(data['Track_pair'].unique())
    idx_set1 = random.choices(list_pair, k=len(list_pair))
    idx_set2 = random.choices(list_pair, k=len(list_pair))

    list_seen1 = []
    c = 0
    bootstrapped_df1 = []
    for pair in idx_set1:
        sub_df1 = data_to_sample[data_to_sample['Track_pair'] == pair].copy()
        if pair in list_seen1:
            sub_df1.loc[:, 'Track_pair'] = pair + '_' + str(c)    # Change the name of the track if it is a duplicated track to have only unique track names.
            c += 1
        else:
            list_seen1.append(pair)
        bootstrapped_df1.append(sub_df1)
    # Rebuild the whole dataframe
    df1 = pd.concat(bootstrapped_df1, ignore_index=True).reset_index(drop=True)

    list_seen2 = []
    c = 0
    bootstrapped_df2 = []
    for pair in idx_set2:
        sub_df2 = data_to_sample[data_to_sample['Track_pair'] == pair].copy()
        if pair in list_seen2:
            sub_df2.loc[:, 'Track_pair'] = pair + '_' + str(c)    # Change the name of the track if it is a duplicated track to have only unique track names.
            c += 1
        else:
            list_seen2.append(pair)
        bootstrapped_df2.append(sub_df2)
    # Rebuild the whole dataframe
    df2 = pd.concat(bootstrapped_df2, ignore_index=True).reset_index(drop=True)
    return df1, df2


def estimate_variance_distribution(data, model):
    """
    Parameters
    ----------
    data : Pandas df
        Bootstrapped set of time series containing anchor coordinates in XYZ.
    model : Models
        Model of coordinate difference.

    Returns
    -------
    var_loop : float
        Variance of distribution of coordinate difference.

    """
    data_SetSigma = read_and_transform_data(data)
    data_hist = Utils_coordinate_difference.makeHistogram(data_SetSigma, bins)
    data_hist = np.divide(data_hist, np.sum(data_hist))

    # estimate variance
    init_vals = [np.sum(data_hist) / np.sum(model.gaussian(xyz, 1, 0.0001)), 0.1]
    best_vals, covar = curve_fit(model.gaussian, xyz, data_hist, p0=init_vals, maxfev=10000, bounds=(0, np.inf))
    var_loop = max(best_vals[1], 0)
    return var_loop


def read_and_transform_data(data_to_read):
    
    Frame_full_,X0_full_,Y0_full_,Z0_full_,X1_full_,Y1_full_,Z1_full_ = Utils_coordinate_difference.readLoc_withoutImport(data_to_read)
    XYZ0_full_ = Utils_coordinate_difference.myArrayFlatten(X0_full_,Y0_full_,Z0_full_)
    XYZ1_full_ = Utils_coordinate_difference.myArrayFlatten(X1_full_, Y1_full_, Z1_full_)
    transformed_data = Utils_coordinate_difference.transform_data(XYZ0_full_, XYZ1_full_)
    return transformed_data


def precompute_and_infer_fractions(data, model_closed, model_noAux, model_auxin, var_closed, var_auxin):
    """
    Parameters
    ----------
    data : Pandas DataFrame
        Bootstrapped set of time series containing anchor coordinates in XYZ.
    model_closed : Models
        Model of proximal coordinate differences.
    model_noAux : Models
        Model of - Auxin coordinate differences.
    model_auxin : Models
        Model of + Auxin coordinate differences.
    var_closed : float
        Variance of proximal coordinate difference distribution.
    var_auxin : float
        Variance of + Auxin coordinate difference distribution.

    Returns
    -------
    tuple
        Fraction of (proximal, extruding, open) states.
    """
    ###############################
    # Precompute models
    ###############################
    precomp_open_data = model_auxin.gaussian(xyz, 1, var_auxin)
    precomp_open_data /= np.sum(precomp_open_data)
    def precomp_free(x, amplitude):
        return np.multiply(precomp_open_data, amplitude)

    precomp_closed_data = model_closed.gaussian(xyz, 1, var_closed)
    precomp_closed_data /= np.sum(precomp_closed_data)
    def precomp_loop(x, amplitude):
        return np.multiply(precomp_closed_data, amplitude)

    precomp_extr_data = model_noAux.gaussian_integrated_fast(np.array(xyz), 1, var_closed, var_auxin, False)  # Put false if does not want to take into account the two different extrusion speeds
    precomp_extr_data /= np.sum(precomp_extr_data)
    def precomp_extr(x, amplitude):
        return np.multiply(precomp_extr_data, amplitude)

    def precomp_full(x, ampl_loop, ampl_free, ampl_extr):  # gaussian model 3d
        return (np.multiply(precomp_closed_data, ampl_loop) + np.multiply(precomp_extr_data, ampl_extr) + np.multiply(
            precomp_open_data, ampl_free))

    ################################
    # Fit 3-state model
    ###############################
    data_Estimate = read_and_transform_data(data)
    data_hist = Utils_coordinate_difference.makeHistogram(data_Estimate, bins)
    data_hist = np.divide(data_hist, np.sum(data_hist))
    init_vals = [.33, .33, .33]
    best_vals, covar = curve_fit(precomp_full, xyz, data_hist, p0=init_vals, bounds=(0, 1))

    ampl_closed_fitted = best_vals[0] / np.sum(best_vals[:])
    ampl_extr_fitted = best_vals[2] / np.sum(best_vals[:])
    ampl_open_fitted = best_vals[1] / np.sum(best_vals[:])
    return (ampl_closed_fitted, ampl_extr_fitted, ampl_open_fitted)


def estimate_fraction_loop_states(file_closed, file_open, file_wt, model_Closed, model_NoAux, model_Auxin, Columns_Keep):
    """
    Parameters
    ----------
    file_closed : Pandas df
        Segmented proximal states. Contains anchor coordinates.
    file_open : Pandas df
        Auxin-treated cells. Contains anchor coordinates.
    file_wt : Pandas df
        Untreated cells. Contains anchor coordinates.
    model_Closed : Model of proximal states
    model_NoAux :  Model
    model_Auxin : Model
    Columns_Keep : list of str
        Column names that need to be kept for the pipeline.

    Returns
    -------
    fitted_closed : tuple
        Contains the estimated fraction of (proximal, extruding, open) states in the segmented proximal states.
    fitted_WT : tuple
        Contains the estimated fraction of (proximal, extruding, open) states in the untreated cells.
    fitted_open : tuple
        Contains the estimated fraction of (proximal, extruding, open) states in the auxin-treated cells.
    """
    # Deletes unecessary columns
    file_WT = Keep_Only_Desired_Columns(file_wt, Columns_Keep)
    file_Open = Keep_Only_Desired_Columns(file_open, Columns_Keep)
    file_Closed = Keep_Only_Desired_Columns(file_closed, Columns_Keep)

    ###############################
    # Create two bootstrapped sets and estimate variance on one bootstrapped set
    ###############################
    Set_Unused_WT, Set_Estimate_WT = split_dataset_per_track_WholeData(file_WT)
    Set_Sigma_Closed, Set_Estimate_Closed = split_dataset_per_track_WholeData(file_Closed)
    Var_closed = estimate_variance_distribution(Set_Sigma_Closed, model_Closed)
    Set_Sigma_Open, Set_Estimate_Open = split_dataset_per_track_WholeData(file_Open)
    Var_auxin = estimate_variance_distribution(Set_Sigma_Open, model_Auxin)

    ###############################
    # Fit 3-state model on the second bootstrapped set
    ###############################
    fitted_WT = precompute_and_infer_fractions(Set_Estimate_WT, model_Closed, model_NoAux, model_Auxin, Var_closed, Var_auxin)
    fitted_open = precompute_and_infer_fractions(Set_Estimate_Open, model_Closed, model_NoAux, model_Auxin, Var_closed, Var_auxin)
    fitted_closed = precompute_and_infer_fractions(Set_Estimate_Closed, model_Closed, model_NoAux, model_Auxin, Var_closed, Var_auxin)
    return fitted_closed, fitted_WT, fitted_open


File_WT = tracks_WT.copy()
File_Open = tracks_deplete.copy()
File_Closed = filtered_WT[filtered_WT['Proximal_filt'] == 1]

######################################################################
# Localization precision
######################################################################
sigma_x_green_NoAuxin = np.nanmean(File_WT['Precision_2_X'])
sigma_y_green_NoAuxin = np.nanmean(File_WT['Precision_2_Y'])
sigma_z_green_NoAuxin = np.nanmean(File_WT['Precision_2_Z'])
sigma_x_red_NoAuxin = np.nanmean(File_WT['Precision_1_X'])
sigma_y_red_NoAuxin = np.nanmean(File_WT['Precision_1_Y'])
sigma_z_red_NoAuxin = np.nanmean(File_WT['Precision_1_Z'])

sigma_x_green_Auxin = np.nanmean(File_Open['Precision_2_X'])
sigma_y_green_Auxin = np.nanmean(File_Open['Precision_2_Y'])
sigma_z_green_Auxin = np.nanmean(File_Open['Precision_2_Z'])
sigma_x_red_Auxin = np.nanmean(File_Open['Precision_1_X'])
sigma_y_red_Auxin = np.nanmean(File_Open['Precision_1_Y'])
sigma_z_red_Auxin = np.nanmean(File_Open['Precision_1_Z'])

sigma_x_green_Closed = np.nanmean(File_Closed['Precision_2_X'])
sigma_y_green_Closed = np.nanmean(File_Closed['Precision_2_Y'])
sigma_z_green_Closed = np.nanmean(File_Closed['Precision_2_Z'])
sigma_x_red_Closed = np.nanmean(File_Closed['Precision_1_X'])
sigma_y_red_Closed = np.nanmean(File_Closed['Precision_1_Y'])
sigma_z_red_Closed = np.nanmean(File_Closed['Precision_1_Z'])

# Define each model, taking into account localization precisions
Var_x_NoAuxin, Var_y_NoAuxin, Var_z_NoAuxin = compute_var(sigma_x_green_NoAuxin, sigma_x_red_NoAuxin, sigma_y_green_NoAuxin, sigma_y_red_NoAuxin, sigma_z_green_NoAuxin, sigma_z_red_NoAuxin)
Var_x_Auxin, Var_y_Auxin, Var_z_Auxin = compute_var(sigma_x_green_Auxin, sigma_x_red_Auxin, sigma_y_green_Auxin, sigma_y_red_Auxin, sigma_z_green_Auxin, sigma_z_red_Auxin)
Var_x_Closed, Var_y_Closed, Var_z_Closed = compute_var(sigma_x_green_Closed, sigma_x_red_Closed, sigma_y_green_Closed, sigma_y_red_Closed, sigma_z_green_Closed, sigma_z_red_Closed)
Model_NoAuxin = Models(Var_x_NoAuxin, Var_y_NoAuxin, Var_z_NoAuxin)
Model_Auxin = Models(Var_x_Auxin, Var_y_Auxin, Var_z_Auxin)
Model_Closed = Models(Var_x_Closed, Var_y_Closed, Var_z_Closed)

# Parameters of the analytical model
Columns_To_Keep = ['Track_pair', 'Track_1_id', 'Track_1_id.1', 'Frame', 'Spot_1_X', 'Spot_1_Y', 'Spot_1_Z', 'Spot_2_X', 'Spot_2_Y', 'Spot_2_Z']
maxibin = 2.5        # total width of distributions
binstep = 0.01       # Binning size
bins = np.arange(-maxibin, maxibin + binstep, binstep)
binsub = bins[0:-1] + binstep / 2
xyz = np.array([binsub] * 3).ravel()
lenvect = len(xyz)//3

Fitted_Closed, Fitted_WT, Fitted_Open = estimate_fraction_loop_states(File_Closed, File_Open, File_WT, Model_Closed, Model_NoAuxin, Model_Auxin, Columns_To_Keep)

# Print results
print('SEGMENTED PROXIMAL STATES: Fraction Proximal=%.3f, Fraction extruding=%.2f, Fraction open=%.2f' % (Fitted_Closed[0], Fitted_Closed[1], Fitted_Closed[2]))
print('UNTREATED CELLS: Fraction Proximal=%.3f, Fraction extruding=%.2f, Fraction open=%.2f' % (Fitted_WT[0], Fitted_WT[1], Fitted_WT[2]))
print('AUXIN-TREATED CELLS: Fraction Proximal=%.3f, Fraction extruding=%.2f, Fraction open=%.2f' % (Fitted_Open[0], Fitted_Open[1], Fitted_Open[2]))


end = time.time()
print(end - start)