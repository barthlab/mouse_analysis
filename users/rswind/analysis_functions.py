# Functions for trace analysis
# Last updated: 09/25/2023 (refactor code)
# Author: Rachel Swindell

import pandas as pd
import numpy as np

################################
#    SHORT HELPER FUNCTIONS    #
################################

# helper functions for lick frequency

def puff_delta(data, index, key_time, key_delay, name):
    '''Caluculates time from sample in a trial to air delivery in that trial and
       returns data in new column of the given name.
    
    Assumes SAT trial structure: IR beam break, then random delay between 200 
    and 800 ms, then 500ms of air puff delivery, then 500ms fixed delay, then 
    water delivery and timeout for 2000 ms before next trial. All samples during
    the random delay period have delay recorded. The delay value is taken 
    from the first sample in a trial, which is guaranteed to be during this
    random delay period. Start of air delivery is calculated as the time of the
    first sample in a trial plus the length of the random delay, which is then 
    subtracted from all samples in the trial.


    Parameters:
    data --- data frame containing trial data
    index --- list of columns to group data by
    key_time --- name of column containg sample time
    key_delay --- name of column containing delay time for that trial
    name --- name of column to write delta to
    '''
    # setting group columns as index and making air_start column required for 
    # row alignment for addition + subtraction
    data = data.set_index(index)
    grouped = data.groupby(index)
    data["air_start"] = grouped[key_time].first() + grouped[key_delay].first()
    data[name] = data[key_time] - data["air_start"]
    data = data.drop(columns="air_start")

    return data.reset_index()

def get_first_sample(data, index, key, name):
    '''Creates a columns with first row of a given column row across all samples
    when grouped by a given index of columns.

       Used to copy timestamp at trial initation across all samples in a trial.

       Parameters:
       data --- data frame containing trial data
       index --- list of columns to group data by
       key --- column to take start time value from
       name --- name of column to write first sample to'''

    data = data.set_index(index)
    data[name] = data.groupby(index)[key].first()
    data = data.reset_index()
    return data

def rolling_time(data, freq_window, index, values_mean, values_first, key):
    '''Calculate licking frequency across index groups as a time-based rolling 
       average and return as a dataframe.

    Mutually exclusive with resample_align. Window is in ms and is open on the 
    left and closed on the right.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns 
             in values, keep, and index
    freq_window --- length of window to calculate rolling window with in ms
    values --- list of columns to calculate rolling average
    keep --- list of columns to keep but not calculate rolling average
    index --- list of columns to index/group by
    key --- column to calculate rolling window with (datetime)
    '''
    # can only do rolling window calculations over index
    group = data.set_index(key).groupby(index)
    # licking frequency in Hz = number of licks/time window 
    licks = group.rolling(window=pd.to_timedelta(freq_window, unit="ms"), closed="right")[values_mean].sum()/(freq_window/1000)
    keep = data.set_index(index)[values_first].set_index(key,append=True)
    data_roll = pd.concat([licks, keep], axis=1).reset_index()

    return data_roll

def bin_by_time(data, bin_size, bin_unit, index, values_mean, values_first, key):
    '''Groups trials by time in given size and unit, and returns data with key 
       column replaced with time bin.

    Takes average across time bin of columns in values_mean; if bin_unit is 
    'ms', calculates frequecny based on given bin_size. Takes first sample of 
    columns in values_first; if index columns group rows in such a way that each
    row is its own unique group, resamples data without combining or dropping
    any rows.
    
    Parameters:
    data --- data frame containing trial data
    bin_size --- size of bin to group by, in bin_unit (int)
    bin_unit --- unit of bin to group by (e.g. ms, min, hr)
    index --- list of columns to group data by
    key --- column to group with, replaced with time bin'''

    gp = pd.Grouper(key=key, freq=pd.to_timedelta(bin_size, unit=bin_unit))
    # append resampler to index group
    index = index.copy()
    index.append(gp)
    group = data.groupby(index)
    
    if bin_unit == "ms":
        licks = group[values_mean].sum()/(bin_size/1000)
    else:
        licks = group[values_mean].mean()

    keep = group[values_first].first()
    data_sample = pd.concat([licks, keep], axis=1).reset_index()

    return data_sample

def delta(data, index, gp, key, key_acc, name):
    '''Calculates time from start of SAT to trial and returns data in given 
    column.

    For the given group, subtracts the values of key_time and key_acc
    in the first row of the group from all values of key_time in the group.
    If used to calculate time from start of SAT, must provide the timebin each
    row is in as the key.

    Parameters:
    data --- data frame containing trial data
    index --- list of columns to index data
    gp --- column or list of columns to
    key --- column with bin times to use for calculation
    name --- name of column to write delta to'''

    data = data.set_index(index)
    # sort to ensure that group takes first sample by ordered key in group
    first_sample_group = data.sort_values(key).groupby(gp)[key].first() 
    # subtract to align delta to start of SAT, rather than to start of train
    first_sample_acc = data.groupby(gp)[key_acc].first()
    delta = data[key] - first_sample_group - first_sample_acc
    data[name] = delta
    return data.reset_index()

# helper functions for timebin aggregation

def drop_group(data, min_trials, min0, min1, index, col, cond0, cond1, key):
    '''Drops bins with fewer than the given number of trials and returns data 
       without those bins.

    Filters bins based on total number of samples in index-based groups, 
    as well as number of samples from given key in a group. This key must be 
    binary (e.g. stimulus or blank; water or no water) for the key specific 
    functionality to work. That is, cond0 and cond1 must be from the same
    column, mutually exclusive, and have no other possible values.

    Parameters:
    data --- data frame containing trial data
    min_trials --- minimum total number of trials to keep a bin (exclusive)
    min0 --- minimum number of trials per time bin for the first trial type
             (exclusive)
    min1 --- minimum number of trials per time bin for the second trial type
             (exclusive)
    index --- set of columns to group data by
    col --- column to take trial number counts from 
            (has no effect on analysis, should not be in index)
    cond0 --- name of first trial type
    cond1 --- name of second trial type
    key --- name of column containing trial type labels
    '''
    gp_index = index.copy()
    gp_index.append(key)
    group = data.groupby(gp_index)
    data = data.set_index(gp_index)

    #filter bins with fewer than min_blank blank or min_water water trials
    counts = group.count().unstack(level=key)
    cond = (counts.loc[:, (col, cond0)] > min0) & (counts.loc[:, (col, cond1)] > min1) 
    data = data[cond]

    #filter bins with fewer than min_trials total trials
    total_group = data.groupby(index) 
    total_counts = total_group.count()
    data = data[total_counts[col] > min_trials]
    data = data.reset_index()

    return data

def mean_bin(data, index, value, keep):
    '''Take mean of data in value and take first row of data in keep, 
       grouped by index, and return aggredate values.'''

    data = data.groupby(index)
    lick = data[value].mean()
    keep = data[keep].first()

    data = pd.concat([lick, keep], axis=1).reset_index()
    return data

def thresh(data, start, end, key):
    '''Drop rows with values of the key column outside of start and end. 
    Requires key to be of type timedelta. Inclusive to start and exclusive of 
    end.'''

    s = pd.to_timedelta(abs(start),unit="ms")
    if start < 0:
        s = -s

    e = pd.to_timedelta(abs(end), unit="ms")
    if end < 0:
        e = -e
    
    data = data[(data[key] >= s) & (data[key] < e)]
    return data

def performance(data, index, keep, key_trial, key_value, cond0, cond1):
    '''Calculate performance (difference in aggregate value between trial types)
    and return as new data frame.

    Groups data by index and given trial type column, then computes 
    cond1 - cond0 for each group to get performance. Takes the first row from 
    cond1 group for kept columns (metadata).
    
    Parameters:
    data --- data frame containing trial data
    index --- list of columns to group data by
    keep --- list of columns to keep
    key_trial --- column containing trial type labels
    key_value --- column containing data to subtract
    cond0 --- name of first trial type
    cond1 --- name of second trial type'''

    group = data.set_index(index).groupby(key_trial)
    blank = group.get_group(cond0)[key_value]
    stim = group.get_group(cond1)[key_value]
    data_perf = stim - blank
    keep = group.get_group(cond1)[keep]
    data_perf = pd.concat([data_perf, keep], axis=1)
    return data_perf.reset_index()

def trial_counts(data, index, keep, key):
    '''Returns total number of rows in key for data grouped by index. 
    
    Keeps first sample from given list of columns in keep. The key column used 
    to select counts does not affect the result as long as it is not in index.'''

    gp = data.groupby(index)
    counts = gp[key].count()
    keep = gp.first()[keep]

    data = pd.concat([counts, keep], axis=1).reset_index()
    return data

def get_bout_start_ind(gp, key, bt_ln):
    bout_start_ind = np.where(gp[key] - gp[key].shift() > bt_ln)[0]
    bout_start_ind = bout_start_ind + gp.index[0]
    bout_start_ind = np.insert(bout_start_ind, 0, gp.index[0])
    bout_start_ind = np.append(bout_start_ind, gp.index[-1])

    return bout_start_ind

def get_bouts(data, group, key, name, bout_len, bout_unit):
    '''Labels trials in bins based on the given bin size.
    
    Two trials are placed in separate bouts if time between the end of one trial
    and the start of the next is more than bout_len.

    Parameters:
    data --- data frame containing trial data
    group --- list of columns to group data by
    key --- column to get sample times from
    name --- name of column to write bout numbers to
    bout_len --- time between two trials to consider them in separate bouts
    bout_unit --- unit of bout_len
    '''

    bt_ln = pd.to_timedelta(bout_len, unit=bout_unit)
    gps = data.groupby(group)
    res = []

    for gp_name, gp in gps: 
        bout_starts = get_bout_start_ind(gp, key, bt_ln)

        gp[name] = pd.NA

        for i in range(1, len(bout_starts)):
            gp.loc[bout_starts[i-1]:bout_starts[i],name] = i - 1

        gp[group] = gp_name
        res.append(gp)

    res = pd.concat(res)
    return res

def get_n_back(data, trial_ind, bout_ind, key, n, cond0, cond1):

    # get first sample for each trial
    trial_start = data.groupby(trial_ind,as_index=False).first()
    bt_gp = trial_start.groupby(bout_ind)[key]
    data = data.set_index(trial_ind)

    # label each trial - need to label each sample
    for i in range(1, n+1):
        cname = f'{i}_back'
        
        # shift to get n-previous trial identity
        trial_start[cname] = bt_gp.shift(i)
        
        # get indicies of trials that had n-previous respective n-back trial identity
        cond0_ind = pd.Index(trial_start.loc[(trial_start[cname] == cond0), trial_ind])
        cond1_ind = pd.Index(trial_start.loc[(trial_start[cname] == cond1), trial_ind])
        
        # label all samples in relevant trials with n-back trial identity
        data.loc[cond0_ind, cname] = cond0
        data.loc[cond1_ind, cname] = cond1

    return data.reset_index()

#################################
#    MAIN ANALYSIS FUNCTIONS    #
#################################


def lickfreq_analysis(data, freq_window, freq_bin, time_bin):
    '''Computes licking frequency at a given sample rate across trial and
    returns raw licking frequency.
    
    Parameters provide flexibility in multiple ways for both cross-trial and 
    cross-training binning.

    Parameters:
    data --- data to be analyzed
    freq_window --- size of window for rolling window analysis across a trial 
                    (milliseconds)
    freq_bin --- size of bin for discrete bin analysis across a trial 
                 (milliseconds)
    time_bin --- size of bin for multi-trial binning over entire training 
                 (minutes)
    '''
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_acc = "acc"

    index_trial = ["condition", "animal","trial no"]

    keep = ["timestamp", "age", "sex", "strain", "acc", "stimulus", "water"]
    values = ["lick", "poke"]


    # do not modify loaded data
    df = data.copy()

    # calculated time from sample to air puff
    data = puff_delta(df, index_trial, key_time, key_delay, key_puff)

    # lick frequency
    if freq_window != 100:
        # rolling window
        keep.append(key_puff)
        data = rolling_time(data, freq_window, index_trial, values, keep, 
                                         key_time)
        keep.pop()
    else:
        # discrete bins
        data = bin_by_time(data, freq_bin, "ms", index_trial, values, keep, key_puff)

    # resampling for cross-trial sample alignment
    data = bin_by_time(data, 100, "ms", index_trial, [], keep + values, key_puff)

    # label start of trial at each sample
    data = get_first_sample(data, index_trial, key_time, key_start)

    # bin by time (minutes)
    index_trial.append(key_puff)
    data = bin_by_time(data, time_bin, "min", index_trial, [], keep + values, key_start)
    index_trial.pop()

    # calulate time bin relative to start of SAT training
    data = delta(data, index_trial, "animal", key_start, key_acc, key_delta)
    
    # convert time bins and trial time to float representations
    data["Time (hr)"] = data[key_delta].to_numpy(dtype="timedelta64[m]").astype("float")/60.
    data["Time (ms)"] = data[key_puff].to_numpy(dtype="timedelta64[ms]").astype("float")

    return data

def aggregate_analysis(data, min_trials, min_blank, min_stimulus):
    '''Computes aggregate values for time bins (mean lick frequency, counts, 
    performance) and return them as separate data frames.
    
    Drops time bins with fewer than the given number of trials. Thresholds
    samples from 200 ms before air delivery (shortest amount of time random 
    delay can be) to 2000 ms after water delivery (shortest amount of time
    until samples for a trial stop).

    Parameters:
    data--- data frame to calculate aggregate values from
    min_trials --- minimum total number of trials to keep a bin
    min_stimulus --- minimum number of stimulus trials to keep a bin
    min_blank --- minimum number of blank trials to keep a bin'''
    
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_acc = "acc"

    index_animal = ["condition", "animal"]
    index_timebin = index_animal + [key_delta, key_puff]

    keep = ["age", "sex", "strain", "acc", "Time (hr)", "Time (ms)"]
    values = ["lick", "poke"]

    # drop bins with fewer than the given number of total, water, or blank trials
    # group by delta for timebins
    # group by key_puff for counting (# samples at any puff delta is number of
    # trials in a bin), if trials are aligned on key_puff
    key = "stimulus"
    cond0 = "blank"
    cond1 = "stimulus"
    data = drop_group(data, min_trials, min_blank, min_stimulus, index_timebin, 
                      "trial no", cond0, cond1, key)
    
    # number of trials for each timebin
    # TODO: delete? (if we have counts_groups we dont need counts) <- check that the groups sum to the total 1st
    # counts = trial_counts(data, index_timebin, keep, "trial no")
    # counts = thresh(counts, 0, 1, key_puff)
    # counts = counts.drop([key_puff, "Time (ms)"], axis=1)

    index_groups = index_timebin + ["stimulus", "water"]
    counts_groups = trial_counts(data, index_groups, keep, "trial no")
    counts_groups = thresh(counts_groups, 0, 1, key_puff)
    counts_groups = counts_groups.drop([key_puff, "Time (ms)"], axis=1)

    # mean licking frequencies

    data_mean = mean_bin(data, index_timebin + ["stimulus"], values, keep + ["water"])

    # threshold trials to 200ms before to 2000ms after air puff
    data_mean = thresh(data_mean, -200, 2000, key_puff)
    
    # performance (Lstim - Lblank)
    perf = performance(data_mean, index_timebin, keep, key, "lick", 
                       cond0, cond1)

    return (data_mean, counts_groups, perf)