# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 20:58:08 2023

@author: thoma
"""

import pandas as pd
import numpy as np
from itertools import groupby
from scipy.optimize import least_squares
from random import choices
import random
import scipy
from scipy import stats, optimize


def sub_add_nan_rows_in_missing_frames_and_interpolate(x, Columns_list, interpolate):
    """
    Parameters
    ----------
    x : Pandas df
        From pd.groupby(). DataFrame of a single track pair.
    Columns_list : list of str
        List of column names.
    interpolate : boolean
        If True, interpolate distances of missing frames.

    Returns
    -------
    x : Pandas df
        Contains nan in missing frames.
    """
    x.reset_index(drop=True, inplace=True)
    idx = x[x['Diff_Frame'] != 1].index.tolist()   # index of frame after a missing frame
    if idx[0] == 0:    
        idx.remove(0)          # first timepoint is nan
    c = 0
    for t in idx:
        t = int(t + c)          # Update the index after having added rows.
        nbr_miss_frame = x['Diff_Frame'][t] - 1
        miss_frame_df = pd.DataFrame(np.nan, index=list(np.arange(0, nbr_miss_frame, 1)), columns=Columns_list)
        miss_frame_df['Track_pair'] = x['Track_pair'].unique()[0]
        miss_frame_df['Track_1_id'] = x['Track_1_id'].unique()[0]
        miss_frame_df['Track_1_id.1'] = x['Track_1_id.1'].unique()[0]
        miss_frame_df['Score'] = 0                          # Set score to 0 for interpolated distances
        miss_frame_df['Frame'] = list(np.arange(x['Frame'].loc[t-1], x['Frame'].loc[t], 1))[1:]

        if interpolate is True:    # interpolate distance values between last and next known timepoints
            miss_frame_df['Distance'] = list(np.linspace(x['Distance'].loc[t-1], x['Distance'].loc[t], len(miss_frame_df)+2, endpoint=True))[1:-1]

        x = pd.concat([x.iloc[:t], miss_frame_df, x.iloc[t:]], ignore_index=True)   # Add df of nan at missing frame location
        c += nbr_miss_frame
    return x


def Add_NA_in_missing_frames_And_interpolate(df, Interpolate):
    """
    Parameters
    ----------
    df : Pandas Dataframe
        Contains the position of spots, Frames and distances.
    Interpolate : boolean
        If True, interpolates the distance between missing frames

    Returns
    -------
    df_na : Pandas Dataframe
        Same as input df but with NAN added in missing frames.
    """
    df.reset_index(drop=True, inplace=True)
    list_pair = list(df['Track_pair'].unique())
    df['Diff_Frame'] = df['Frame'].groupby(df['Track_pair']).diff()
    df['Diff_Dist'] = df['Distance'].groupby(df['Track_pair']).diff()
    df.loc[df['Diff_Frame'] < 0, 'Diff_Frame'] = np.nan
    columns_list = list(df)

    df_na = df.groupby('Track_pair').apply(sub_add_nan_rows_in_missing_frames_and_interpolate, Columns_list=columns_list, interpolate=Interpolate).reset_index(drop=True)
    df_na.drop(columns=['Diff_Dist', 'Diff_Frame'], inplace=True)
    
    if Interpolate is True:
        print(str(1 - len(df)/len(df_na)) + ' fraction of timepoints interpolated')
    return df_na


def get_higher_number_of_consecutive_frames_under_threshold_dist(data, threshold_spatial, Quantile_Temporal):
    
    """
    Parameters
    ----------
    data : Pandas df
        Contains the position of spots, Frames and distances.
    threshold_spatial : float
        Spatial threshold (in um). Used to compute the number of consecutive frames under this threshold.
    Quantile_Temporal : float
        Quantile of the distribution of consecutive frames under spatial threshold to be used to compute the temporal threshold in frames.

    Returns
    -------
    list
        list of durations of consecutive frames under spatial threshold.
    int
        Temporal threshold in frames.
    """
    data['BelowThresh'] = data['Distance'] < threshold_spatial
    data['streak_id'] = (data['BelowThresh'] != data['BelowThresh'].shift(1)).cumsum()
    streak_groups = data.groupby(['Track_pair','BelowThresh','streak_id'])["Frame"].agg(['min','max']).reset_index()                # df containing length of intervals below or above the spatial threshold (True or False)
    positive_streaks = streak_groups[streak_groups['BelowThresh'] & (streak_groups['min'] != streak_groups['max'])].copy()          # Keeps only values below spatial threshold and where the distance is below the spatial threshold for more than one frame.
    positive_streaks["duration"] = positive_streaks["max"] - positive_streaks["min"] + 1                                            # Duration of each interval below threshold
    return list(positive_streaks["duration"]), round(np.quantile(positive_streaks['duration'], Quantile_Temporal))


def add_proximal_information(x, Positive_groups, Temporal_Threshold):
    """
    Parameters
    ----------
    x : Pandas df
        From pd.groupby(). DataFrame of a single track pair.
    Positive_groups : Pandas df
        Duration, start, end of proximal states.
    Temporal_Threshold : int
        Temporal threshold (in frames).

    Returns
    -------
    x : Pandas df
        Raw proximal state segmentation is indicated in the column 'Proximal'.
    """
    pair = x['Track_pair'].unique()[0]
    Positive_pair = Positive_groups[Positive_groups['Track_pair'] == pair]
    Positive_pair = Positive_pair[Positive_pair['duration'] > Temporal_Threshold]    # Keep only intervals where distances are below the spatial threshold for longer than the temporal threshold
    x.loc[:, 'Proximal'] = [0] * len(x)
    List_Frames_Positive = []
    for idx, row in Positive_pair.iterrows():
        List_Frames_Positive.extend(np.arange(row['min'], row['max'] + 1, 1))
    Mask_Frame = x['Frame'].isin(List_Frames_Positive)
    x.loc[:, 'Proximal'] = np.where(Mask_Frame, 1, 0)
    x.reset_index(drop=True, inplace=True)
    return x


def get_closed_state_in_tracks(data, dist_thresh, time_thresh):
    """
    Parameters
    ----------
    data : Pandas df
        Contains the position of spots, Frames and distances.
    dist_thresh : float
        Spatial threshold (in um).
    time_thresh : int
        Temporal threshold (in frames).

    Returns
    -------
    df_masked : Pandas df
        Contains raw proximal state segmentation.
    closed_state_begin : dict
        Contains the frame number of the start (first timepoint) of proximal state.
    closed_state_start : dict
        Contains the frame number of the last timepoint of proximal state.
    """
    data.reset_index(drop=True, inplace=True)
    data['BelowThresh'] = data['Distance'] < dist_thresh
    data['streak_id'] = (data['BelowThresh'] != data['BelowThresh'].shift(1)).cumsum()
    streak_groups = data.groupby(['Track_pair','BelowThresh','streak_id'])["Frame"].agg(['min','max']).reset_index()            # df containing length of intervals below or above the spatial threshold (True or False)
    positive_streaks = streak_groups[streak_groups['BelowThresh'] & (streak_groups['min'] != streak_groups['max'])].copy()      # Keeps only values below spatial threshold and where the distance is below the spatial threshold for more than one frame.
    positive_streaks["duration"] = positive_streaks["max"] - positive_streaks["min"] + 1                                        # Duration of each interval below threshold
    
    # Save start and end of closed states
    closed_state_begin = {}
    closed_state_end = {}
    Long_positive = positive_streaks[positive_streaks['duration'] > time_thresh] 
    for pair in list(data['Track_pair'].unique()):
        closed_state_begin[pair] = {}
        closed_state_end[pair] = {}
        if pair in list(Long_positive['Track_pair']):
            Long_positive_pair = Long_positive[Long_positive['Track_pair'] == pair]
            for idx, row in Long_positive_pair.iterrows():
                closed_state_begin[pair][row['min']] = [row['min']]
                closed_state_end[pair][row['max']] = [row['max']]

    # Add 0 if not proximal, 1 if proximal state
    df_masked = data.groupby('Track_pair').apply(add_proximal_information, Positive_groups=positive_streaks, Temporal_Threshold=time_thresh).reset_index(drop=True)
    df_masked.drop(columns=['streak_id', 'BelowThresh'], inplace=True)
    df_masked.reset_index(drop=True, inplace=True)
    return df_masked, closed_state_begin, closed_state_end



def start_stop(a, trigger_val):
    """
    Parameters
    ----------
    a : list
        Contains proximal state segmentation.
    trigger_val : boolean
        Value to be found.

    Returns
    -------
    list
        Returns index of start and end of consecutives sequences of 'trigger_val'
        First list contains starting index of closed states.
        Second list contains end index of closed states
    """
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False, np.equal(a, trigger_val), False]
    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])
    # Get the start and end indices with slicing along the shifting ones
    return [idx[::2], idx[1::2]-1]


def get_indices_of_closed_state_starts_and_ends_adapted_for_1_timepoint_closed_state(x):
    """
    Parameters
    ----------
    x : Pandas df
        Contains the position of spots, Frames and distances for a single track pair.

    Raises
    ------
    ValueError
        There was not the same number of proximal state start and ends.

    Returns
    -------
    output_in_frame_start : list
        Indexes of proximal state start.
    output_in_frame_end : list
        Indexes of proximal state end.
    """
    isProximal = x['Proximal_filt'].eq(1)
    output = start_stop(isProximal, True)
    indexes_start = x.index[output[0]]
    indexes_end = x.index[output[1]]
    output_in_frame_start = list(x.loc[indexes_start, 'Frame'])  # Frame number corresponding to the indexes of proximal state start.
    output_in_frame_end = list(x.loc[indexes_end, 'Frame'])     # Frame number corresponding to the indexes of proximal state end.

    # If no proximal state of 1 timepoint
    if len(output_in_frame_start) == len(output_in_frame_end):
        tmp = 0
    # If there is a proximal state of 1 timepoint
    else:
        raise ValueError('Not same number of start and end of proximal states')
    return output_in_frame_start, output_in_frame_end


def convolve_proximal(x, temporal_threshold):
    """
    Parameters
    ----------
    x : pandas df from groupby
        Unique track pair.
    temporal_threshold : int
        Temporal threshold.

    Returns
    -------
    x : pandas df from groupby
        Single track with filtered proximal states.
        The binary proximal column is convolved with 1 of len(time_thresh).
    """
    masked_filt = np.array(x.loc[:, 'Proximal'].rolling(window=temporal_threshold, min_periods=1, center=True).apply(np.nanmean))
    masked_filt_binary = [0] * len(masked_filt)
    for j in range(len(masked_filt_binary)):
        if (masked_filt[j] >= 0.5) and (j == 0):                             # first iteration
            masked_filt_binary[j] = 1
        if masked_filt[j] >= 0.5:
            masked_filt_binary[j] = 1
        if j > 0:
            if (masked_filt[j] >= 0.5) and (masked_filt[j-1] < 0.5):
                masked_filt_binary[j] = 1
        if (masked_filt[j] > 0.5) and (j == len(masked_filt_binary) - 1):     # last iteration
            masked_filt_binary[j] = 1
        if (masked_filt[j] == 0.5) and (j == len(masked_filt_binary) - 1):     # last iteration
            masked_filt_binary[j] = 0
        if (j != len(masked_filt_binary) - 1) and (j != 0):
            if (masked_filt[j] == 0.5) and (masked_filt[j-1] > 0.5) and (masked_filt[j+1] < 0.5):
                masked_filt_binary[j] = 0
            if (masked_filt[j] == 0.5) and (masked_filt[j-1] >= 0.5) and (masked_filt[j+1] <= 0.5):
                masked_filt_binary[j] = 0
    x.loc[:, 'Proximal_filt'] = masked_filt_binary
    return x


def mean_filter_masked_tracks(masked, time_thresh):
    """
    Parameters
    ----------
    masked : Pandas df
        Contains the position of spots, Frames and distances and raw proximal state segmentation.
    time_thresh : int
        Temporal threshold (in frames) used as the mean averaging window).

    Returns
    -------
    df_filtered : Pandas df
        Contains the position of spots, Frames and distances and proximal state segmentation.
    closed_state_begin : dict of dict
        Contains the index of the proximal state start(s) for each track pair.
    closed_state_end : dict of dict
        Contains the index of the proximal state end(s) for each track pair.
    Censored : dict of dict
        True if the corresponding proximal state was censored, False otherwise.
    """
    df_filtered = masked.groupby('Track_pair', group_keys=False).apply(convolve_proximal, temporal_threshold=time_thresh).reset_index(drop=True)
    df_firstFrame = df_filtered.groupby('Track_pair')['Frame'].first()
    df_firstFrame_trackpair = pd.DataFrame({'Track_pair':df_firstFrame.index, 'Frame':df_firstFrame.values})
    df_lastFrame = df_filtered.groupby('Track_pair')['Frame'].last()
    df_lastFrame_trackpair = pd.DataFrame({'Track_pair':df_lastFrame.index, 'Frame':df_lastFrame.values})

    list_track = list(df_filtered['Track_pair'].unique())
    closed_state_begin = {}
    closed_state_end = {}
    Censored = {}
    for pair in list_track:
        sub_df = df_filtered[df_filtered['Track_pair'] == pair]
        firstframe = df_firstFrame_trackpair[df_firstFrame_trackpair['Track_pair'] == pair]['Frame'].iloc[0]
        lastframe = df_lastFrame_trackpair[df_lastFrame_trackpair['Track_pair'] == pair]['Frame'].iloc[0]
        start_sub, end_sub = get_indices_of_closed_state_starts_and_ends_adapted_for_1_timepoint_closed_state(sub_df)

        closed_state_begin[pair] = {}
        closed_state_end[pair] = {}
        closed_state_begin[pair] = start_sub
        closed_state_end[pair] = end_sub
        
        # Add dictionary of censored closed state (True if censored, False if not censored)
        Censored[pair] = []
        if len(start_sub) != len(end_sub):
            raise ValueError('Not same number of start and end of closed states')
        for i in range(len(start_sub)):
            if (start_sub[i] == firstframe) or (end_sub[i] == lastframe):
                Censored[pair].extend([True])
            else:
                Censored[pair].extend([False])

    return df_filtered, closed_state_begin, closed_state_end, Censored


def get_lifetime_closed_states(Censored, Start, End):
    """
    Parameters
    ----------
    Censored : dict
        Whether each proximal state is censored or not. Keys are the tracks.
    Start : dict
        First timepoint of proximal states (in frames). Keys are the tracks.
    End : dict
        Last timepoint of proximal states (in frames). Keys are the tracks.

    Returns
    -------
    lifetimes : dict
        Proximal state lifetimes.
    lifetimes_list_flatten : list
        Flattened list of lifetimes.
    Censored : dict
        Unmodified from input.
    censored_list_flatten : list
        Flattened list of censoring. Indexes are the same as the flatten list of lifetimes
    """
    lifetimes = {}
    for pair in Start.keys():
        lifetimes[pair] = []
        for i in range(len(Start[pair])):
            lifetimes[pair].extend([End[pair][i] - Start[pair][i] + 1])
    
    lifetimes_list_flat = [lifetimes[k] for k in lifetimes.keys()]
    lifetimes_list_flatten = []
    for p in range(len(lifetimes_list_flat)):
        lifetimes_list_flatten.extend(lifetimes_list_flat[p])

    censored_list_flat = [Censored[k] for k in Censored.keys()]
    censored_list_flatten = []
    for p in range(len(censored_list_flat)):
        censored_list_flatten.extend(censored_list_flat[p])

    return lifetimes, lifetimes_list_flatten, Censored, censored_list_flatten
    

def get_closed_state_frequency(data, frame_s):
    """
    Parameters
    ----------
    data : Pandas df
        Contains the position of spots, Frames and distances and proximal state segmentation.

    Returns
    -------
    mean_frequency : float
        Frequency of appearance of proximal states (if no proximal state, count frequency of 0).
    """
    list_track_pair = list(data['Track_pair'].unique())
    n_occurences = 0
    length_tot = 0
    for pair in list_track_pair:
        sub_df = data[data['Track_pair'] == pair]
        if np.isnan(sub_df['Proximal_filt']).any() == False:
            n_occurences += sum(1 for k, _ in groupby(list(sub_df['Proximal_filt'])) if k)
            length_tot += len(sub_df['Proximal_filt']) * frame_s
    mean_frequency = n_occurences / length_tot
    return mean_frequency


def retrieve_closed_states_without_closed_states_before_wdw_timepoints_fulltracksOnly_WithScores(masked, start_closed, wdw):
    """
    Parameters
    ----------
    masked : Pandas DF
        Contains the position of spots, Frames and distances and proximal state segmentation.
    start_closed : dict of dict
        Contains the index of the proximal state start(s) for each track pair.
    wdw : int
        Indicates the size of the window used to fit closing rates (before the closed state start).

    Returns
    -------
    closed : dict
        For each track, contains the track id and the 'wdw' distances before the inferred proximal state(s).
        Last value of each array is the first proximal state timepoint.
    scores : dict
        For each track, contains the track id and the 'wdw' scores before the inferred proximal state(s).
        Last value of each array is the first proximal state timepoint.
    """
    closed = {}
    scores = {}
    list_track_pair = list(masked['Track_pair'].unique())
    for pair in list_track_pair:
        sub_df = masked[masked['Track_pair'] == pair].reset_index(drop=True)
        if 1 in sub_df['Proximal_filt'].values:
            closed[pair] = {}
            scores[pair] = {}
            for start in start_closed[pair]:
                start_frame = start - wdw - 1
                closed_state = sub_df[sub_df['Frame'].between(start_frame, start)]
                if ((1 not in closed_state.loc[start_frame:start-1, 'Proximal_filt'].values) and (len(closed_state.loc[start_frame: start - 1, 'Proximal_filt']) == wdw + 1)):       # Removes cases where a closed state is preceded by a proximal state in wdw timepoints before proximal state start
                    closed[pair][start] = closed_state['Distance'].values
                    scores[pair][start] = closed_state['Score'].values
    return closed, scores


def make_array_of_distances_scores_from_dict(dict_dist, dict_score, wdw):
    """
    Parameters
    ----------
    dict_dist : dict
        Distances aligned on closed states. Keys are Track_pair ids.
    dict_score : dict
        Scores aligned on closed states. Keys are Track_pair ids.
    wdw : int
        Indicates the size of the window used to fit closing rates (before the closed state start).       

    Returns
    -------
    Closed_arr : array
        Contains distances aligned on proximal state start (last timepoint).
    Score_arr : array
        Contains scores aligned on proximal state start.
    Converts dictionaries of distances and scores into arrays.
    """
    # Count number of full trajectories
    c_full = 0      # Number of full tracks (length equal to window size)
    c_total = 0     # All tracks with a proximal state
    for track in dict_dist.keys():
        c_total += len(dict_dist[track])
        for closed in dict_dist[track].keys():
            if len(dict_dist[track][closed]) == wdw + 2:
                c_full += 1
 
    p = 0
    Closed_arr = np.zeros((wdw + 2, c_total))
    Closed_arr[:] = np.nan
    Score_arr = np.zeros((wdw + 2, c_total))
    Score_arr[:] = np.nan
    for track in dict_dist.keys():
        for closed in dict_dist[track].keys():
            Closed_arr[:len(dict_dist[track][closed]), p] = dict_dist[track][closed]
            Score_arr[:len(dict_score[track][closed]), p] = dict_score[track][closed]
            p += 1
    return Closed_arr, Score_arr


class Score2prec:
    
    def __init__(self,precision):
        self.mean=np.nanmean(precision)
        self.std=np.nanstd(precision)
        
    def __init__(self,mean,std):
        self.mean=mean
        self.std=std

    def getMean(self):
        return self.mean
    
    def getStd(self):
        return self.std
    
    def getGaussianModel(self,xaxis):
        m=np.exp(-0.5*(xaxis-self.mean)*(xaxis-self.mean)/(self.std*self.std))
        m=np.divide(m,np.max(m))
        return m
    
    def computeScoreGaussianRepartionFunctionShifted(self, precision):
        score = 1 - 0.5 * (1 + scipy.special.erf((precision - (self.mean+self.std)) / (np.sqrt(2) * self.std)))
        return score
    
    def GaussianRepartionFunctionShifted2Gaussian(self,score):
        d=scipy.special.erfinv((-(score-1)*2)-1)*(np.sqrt(2)*self.std)+(self.mean+self.std)
        return d





class SpeedFitting:
    def __init__(self,a,b,R2):
        self.a=a
        self.b=b
        self.R2=R2

    def getp012(self,p):
        k=0
        if self.parametersFitted[0]:
            p0=p[k]
            k+=1
        else:
            p0=self.a
        if self.parametersFitted[1]:
            p1=p[k]
            k+=1
        else:
            p1=self.b

        if self.parametersFitted[2]:
            p2=p[k]
            k+=1
        else:
            p2=self.R2
        return p0,p1,p2

    def computeSlope(self,p):
        return (p[2]-p[1])/p[0]

    def _computeSlope_forFit(self,p):
        p0,p1,p2=self.getp012(p)
        return (p2-p1)/p0

    
    def func2slopes(self,x, p=None):

        if p is None:
            p=[self.a,self.b,self.R2]
        y=np.zeros((len(x)))
        y[x<p[0]]=p[2]
        slope=self.computeSlope(p)
        y[x>=p[0]]=slope*x[x>=p[0]]+p[1]
        return y

    #1st order polynomial function
    def _func2slopes_forfit(self,x, p):
        p0,p1,p2=self.getp012(p)
        y=np.zeros((len(x)))
        y[x<p0]=p2
        slope=self._computeSlope_forFit(p)
        y[x>=p0]=slope*x[x>=p0]+p1
        return y


    def fit(self,x_time2D,data,datastd,parametersFitted=[False, False, True], with_multi_init_of_A_and_R2=True):
        """
        Parameters
        ----------
        x_time2D : array
            Timepoints (in seconds).
        data : array
            Mean weighted squared distances.
        datastd : array
            Mean weighted standard error of the mean.
        parametersFitted : Array of booleans, optional
            Whether [t_extr, R2int, R2plateau] should be fitted (True) or fixed (False). The default is [True,True,True].
        with_multi_init_of_textr_and_R2plateau : boolean, optional
            True if perform multi initialization of t_extr and R2int, which minimizes the risk of local minima during fitting. The default is True.

        Returns
        -------
        finalP : float
            Parameters of the piecewise linear model (t_extr, R2int, R2plateau).
        loglik : array
            Log-likelihood.
        """
        if parametersFitted[0] is True:
            print('ERROR: t_extr parameter should not be fitted')
        self.parametersFitted=parametersFitted
        #function to minimize
        def fun_loss(p, x, y):
            func_2_slopes=self._func2slopes_forfit(x, p)
            loss= y-func_2_slopes            
            return loss

        x_time_list=x_time2D.ravel()
        y_data_list=data.ravel()
        datastd_list=datastd.ravel()

        if with_multi_init_of_A_and_R2:
            if parametersFitted[2] is False:
                print('ERROR : R2 parameter should be fitted if with_multi_init_of_A_and_R2 is True')
            step_shift_t_extr=1     # The number of timepoints between different initializations of t_extr
            mini_loss=np.Inf
            mini_lik=np.Inf
            finalP=None
            init_t_extr = np.arange(1, len(x_time2D) - 1, step_shift_t_extr)   # fit a minimum of X timepoints
            best_t_extr = init_t_extr[0]
            for various in init_t_extr:
                self.a=x_time2D[various][0]
                pToFit=[]
                if parametersFitted[0]:
                    pToFit.append(x_time2D[various][0])
                if parametersFitted[1]:
                    pToFit.append(self.b)
                if parametersFitted[2]:
                    pToFit.append(np.nanmean(data[:various]))

                p_2slopes=least_squares(fun_loss, pToFit, args=(x_time_list, y_data_list)).x
                lossinit=np.power(fun_loss(pToFit,x_time_list, y_data_list),2)
                loss=np.power(fun_loss(p_2slopes,x_time_list, y_data_list),2)
                lik=loss
                if np.sum(lik)<mini_loss:
                    mini_loss=np.sum(lik)
                    mini_lik=np.sum(lik)
                    finalP=p_2slopes   
                    best_t_extr = self.a
            # Save the minimal loss over all initialization steps
            loss_mse=mini_loss
            p_2slopes=finalP
            self.a = best_t_extr

        p0,p1,p2=self.getp012(p_2slopes)
        finalP=[p0,p1,p2]

        #update class values
        self.a=p0
        self.b=p1 #b initialized to the last point
        self.R2=p2 #R2 initialized to the mean
        
        log_lik=-0.5*mini_lik
        return finalP,log_lik


def computeBIC(log_lik,k,n):
    bic=np.log(n)*k-2*log_lik 
    return(bic)


def fit_ClosingRate_TakingIntoAccount_LocalizationPrecision(close_states, scores, frame_s, len_loop, mean_R02, mean_Std_prec, textr_is_fitted=False, R2int_is_fitted=False, R2plateau_is_fitted=True, with_multi_init_of_textr_and_R2plateau=True):

    # convert score back to localization precision
    sc_p = Score2prec(mean_Std_prec[0], mean_Std_prec[1])
    loc_prec = sc_p.GaussianRepartionFunctionShifted2Gaussian(scores)
    squared_loc_prec = np.power(loc_prec, 2)

    # Remove localization precision contribution
    closed_states_squared_raw = np.power(close_states, 2)
    closed_states_squared = np.subtract(closed_states_squared_raw, squared_loc_prec)
    # replace infinite values (score of 0)
    closed_states_squared[np.isinf(closed_states_squared)] = closed_states_squared_raw[np.isinf(closed_states_squared)]
    
    # Compute mean weighted squared distance and std
    weighted_dist = np.multiply(closed_states_squared, scores)
    norm_mean_weighted_dist = np.divide(np.sum(weighted_dist, axis=1), np.nansum(scores, axis=1))

    sh = np.shape(closed_states_squared)
    lentrack = sh[0]
    nbtrack = sh[1]
    std_weighted = np.divide(np.sqrt(np.divide(np.sum(np.multiply(np.power(np.transpose(np.subtract(np.transpose(closed_states_squared), norm_mean_weighted_dist)), 2), scores),axis=1), np.sum(scores,axis=1))), np.sqrt(nbtrack))
    x_time = 1 / frame_s * (np.arange(lentrack) - lentrack - 1)

    # get nan values in data and in score
    scoreisnan = np.isnan(scores)  # nan score
    # put score=0 when score is nan (interpolated distance)
    scores[scoreisnan] = 0

    # reshape X from 2D to 1D
    x_time2D = np.transpose([x_time]*nbtrack)
    x_time2DMean = np.transpose([np.mean(x_time2D, axis=1)])

    # parameter initialization
    t_extr_init = np.min(x_time) / 2    # Not taken into account if using multiple initialization ( with_multi_init_of_A_and_R2=True)
    R2int_init = norm_mean_weighted_dist[-1]
    if R2plateau_is_fitted is True:
        R2plateau_init = np.nanmean(norm_mean_weighted_dist[:int(lentrack / 2)])    # Not taken into account if using multiple initialization ( with_multi_init_of_A_and_R2=True)
    else:
        R2plateau_init = R2plateau_init
    
    # Fit
    # 3-parameter model with a plateau and a linear slope
    s_f = SpeedFitting(t_extr_init, R2int_init, R2plateau_init)
    parameters_2slopes, loglik_2slopes = s_f.fit(x_time2DMean, norm_mean_weighted_dist, std_weighted, [textr_is_fitted, R2int_is_fitted, R2plateau_is_fitted], with_multi_init_of_textr_and_R2plateau)
    slope_2slopes = s_f.computeSlope(parameters_2slopes)
    bic_2slopes=computeBIC(loglik_2slopes, 3, lentrack)
    
    # Compute closing rate based on the 3-parameter model
    S0overR02 = len_loop / mean_R02
    t_extr = parameters_2slopes[0]
    R2int = parameters_2slopes[1]
    R2_estimated = parameters_2slopes[2]
    closing_rate_2slopes = - slope_2slopes * S0overR02

    information_criterion = {}
    information_criterion['3_parameter'] = (bic_2slopes, parameters_2slopes, closing_rate_2slopes)
    return closing_rate_2slopes, R2_estimated, slope_2slopes, x_time, t_extr, R2int, information_criterion


def Keep_Only_Desired_Columns(data, columns_to_keep):
    """
    Parameters
    ----------
    data : Pandas df
        any Pandas DataFrame.
    columns_to_keep : list
        List of column names that should be kept in the output pandas df.

    Returns
    -------
    data_desired : Pandas df
        Same as input DataFrame but only with the columns present in 'columns_to_keep'.
    """
    data_desired = data[columns_to_keep]
    return data_desired


def randomize_timepoints_in_one_track(x):
    """
    Parameters
    ----------
    x : Pandas df
        Contains the position of spots, Frames and distances for a single track pair, from pd.groupby.

    Returns
    -------
    sub_df_rdm : Pandas df
        Contains the randomized positions of spots, distances within a single track pair.
    """
    frames = list(x['Frame'])
    sub_df_rdm = x.sample(frac=1, replace=False)      # Randomly shuffle rows in sub df
    sub_df_rdm.loc[:, 'Frame'] = frames              # Keep initial frame indexes
    return sub_df_rdm


def randomize_each_timepoint_in_df(data):
    """
    Parameters
    ----------
    data : Pandas df
        Contains the position of spots, Frames and distances and proximal state segmentation.

    Returns
    -------
    df_rdm : Pandas df
        Pandas df where each timepoint in each single track was randomized (within a track).
    """
    data_to_shuffle = data.copy()
    data_to_shuffle.reset_index(drop=True, inplace=True)
    df_rdm = data_to_shuffle.groupby('Track_pair', group_keys=False).apply(randomize_timepoints_in_one_track).reset_index(drop=True)
    return df_rdm


def MLE_censored_exponential(data, censored, conf=0.95):
    """
    Unmodified from Gabriele et al, Science, 2022
    MLE estimate for exponential distribution, given right-censored data
    Parameters
    ----------
    data : array-like, dtype=float
        the (scalar) values
    censored : array-like, dtype=bool, same shape as data
        whether the corresponding value in `!data` is censored (``True``) or
        completely observed (``False``).
    conf : float in [0, 1], optional
        the confidence level to use in calculating bounds
    Returns
    -------
    m, low, high : float
        point estimate for the mean of the exponential distribution, and
        confidence bounds.
    """
    data = np.asarray(data).flatten()
    censored = np.asarray(censored, dtype=bool).flatten()

    n = np.count_nonzero(~censored)
    alpha = 1-conf

    # Point estimate
    m = np.sum(data) / n

    # Confidence interval
    c = stats.chi2(1).isf(alpha) / (2*n)
    def fitfun(beta): return np.exp(beta) - 1 - beta - c

    res = optimize.root_scalar(fitfun, bracket=(-c-1, 0))
    if not res.flag == 'converged': # pragma: no cover
        raise RuntimeError("Root finding did not converge for upper confidence interval")
    beta_m = res.root

    res = optimize.root_scalar(fitfun, bracket=(0, 2*np.sqrt(c)))
    if not res.flag == 'converged': # pragma: no cover
        raise RuntimeError("Root finding did not converge for lower confidence interval")
    beta_p = res.root
    return m, m*np.exp(-beta_p), m*np.exp(-beta_m)

