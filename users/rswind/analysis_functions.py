# Functions for trace analysis
# Last updated: 09/25/2023 (refactor code)
# Author: Rachel Swindell

import pandas as pd

# TODO: make functions fully generic or specific to this analysis (ie requiring specific column names)
# which is easier to debug long-term? -> generic is much easier to change column names or conditions as necessary

################################
#    SHORT HELPER FUNCTIONS    #
################################

# refactor finished

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
    '''Creates a columns with first row of a given column row across all samples when grouped by a given index of columns.

       Used to copy timestamp at trial initation across all samples in a trial.

       Parameters:
       data --- data frame containing trial data. Must contain all of the columns in index and key
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
    freq_window --- length of window to calculate rolling window with in ms (int)
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
    data --- data frame containing trial data. Must contain all of the columns in index and key
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
    '''Calculates time from start of SAT to trial and returns data in given column.

    For the given group, subtracts the values of key_time and key_acc
    in the first row of the group from all values of key_time in the group.
    If used to calculate time from start of SAT, must provide the timebin each
    row is in as the key.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and key
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

# TODO: refactor and make generic

def bin_prev_identity(data, index, gp, key, cond_pos, cond_neg, name):
    '''Returns tuple of all trials where the previous trial was blank and trials where the previous trial was stimulus.
       
       Parameters:
       data --- data frame containing trial data. Must contain all of the columns in index and key
       index --- list of columns to group data by
       key --- column to take trial labels from'''
    # get first sample of each trial
    gp = data.groupby(index,as_index=False)
    trial_start = gp.first()

    # align previous trial type to trial by animal and condition
    trial_start["one_back"] = trial_start.groupby(["condition", "animal"], as_index=False)[key].shift()
    
    # select trials for which previous trial was blank
    one_back_blank = trial_start[trial_start["one_back"] == "blank"].set_index(index)
    
    # get all samples in trials for which previous trial was blank
    one_back_blank_trials = data[data.set_index(index).index.isin(one_back_blank.index)]

    # same as above for stimulus
    one_back_water = trial_start[trial_start["one_back"] == "water"].set_index(index)
    
    # get all samples in trials for which previous trial was blank
    one_back_water_trials = data[data.set_index(index).index.isin(one_back_blank.index)]
    return (one_back_blank_trials, one_back_water_trials)

# TODO: allow conversion to float before this function is run (keep column/make generic)
def drop_bins(data, min_trials, min_blank, min_water, index, key):
    '''Drops bins with fewer than the given number of trials and returns data without those bins.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and key
    index --- list of columns to group data by

    '''
    gp_index = index.copy()
    gp_index.append(key)
    group = data.groupby(gp_index)
    data = data.set_index(gp_index)

    #filter bins with fewer than min_blank blank or min_water water trials
    counts = group.count().unstack(level=key)
    cond = (counts.loc[:, ("trial no", "blank")] > min_blank) & (counts.loc[:, ("trial no", "water")] > min_water) 
    data = data[cond]

    #filter bins with fewer than min_trials total trials
    total_group = data.groupby(index) 
    total_counts = total_group.count()
    data = data[total_counts["trial no"] > min_trials]
    data = data.reset_index()

    return data

def mean_bin(data, index, value, keep):
    data = data.groupby(index)
    lick = data[value].mean()
    keep = data[keep].first()

    data = pd.concat([lick, keep], axis=1).reset_index()
    return data

def thresh(data, start, end):
    '''threshold licking samples to samples between start and end.'''
    s = pd.to_timedelta(abs(start),unit="ms")
    if start < 0:
        s = -s

    e = pd.to_timedelta(abs(end), unit="ms")
    if end < 0:
        e = -e
    
    data = data[(data["delivery delta"] >= s) & (data["delivery delta"] < e)]
    return data

def performance(data, index, keep):
    group = data.set_index(index).groupby("stimulus")
    water = group.get_group("water")["lick"]
    blank = group.get_group("blank")["lick"]
    data_perf = water - blank
    data_perf = data_perf
    keep = group.get_group("water")[keep]
    data_perf = pd.concat([data_perf, keep], axis=1)
    return data_perf.reset_index()

def trial_counts(data, index, keep, key):
    gp = data.groupby(index)
    counts = gp[key].count()
    keep = gp.first()[keep]

    data = pd.concat([counts, keep], axis=1).reset_index()
    return data



#################################
#    MAIN ANALYSIS FUNCTIONS    #
#################################

def lickfreq_analysis(data, freq_window, freq_bin, time_bin):
    '''Computes licking frequency across trial and returns raw licking frequency.
    
    Parameters provide flexibility in a multiple ways for both across-trial and across-training binning.

    Parameters:
    data --- data to be analyzed
    freq_window --- size of window for rolling window analysis across a trial (in milliseconds)
    freq_bin --- size of bin for discrete bin analysis across a trial (in milliseconds)
    time_bin --- size of bin for multi-trial binning over entire training (in minutes)
    '''
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_acc = "acc"

    index_trial = ["condition", "animal","trial no"]
    index_trial_type = ["condition", "animal", "trial no", "stimulus"]

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

def aggregate_analysis(data, min_trials, min_blank, min_water):
    '''Computes aggregate values for time bins (mean lick frequency, counts, performance) and return them as separate data frames.
    
    Drops time bins with fewer than the given number of trials prior to aggregation.

    Parameters:
    data--- data frame to calculate aggregate values from
    min_trials --- minimum total number of trials to keep a bin
    min_water --- minimum number of water trials to keep a bin
    min_blank --- minimum number of blank trials to keep a bin'''
    # drop bins with fewer than the thresholded number of total, water, or blank trials
    index = ["condition", "animal", "delta","delivery delta"]
    key = "stimulus"
    data = analysis_functions.drop_bins(data, min_trials, min_blank, min_water, index, key)

    

    # number of trials for each timebin
    index = ["condition", "animal", "delta", "delivery delta"]
    keep = ["age", "sex", "strain", "acc", "Time (hr)", "Time (ms)"]
    counts = analysis_functions.trial_counts(data, index, keep, "trial no")
    
    # mean licking frequencies for each animals for each time bin 
    keep = ["age", "sex", "strain", "acc","Time (hr)", "Time (ms)"]
    index = ["condition", "animal", "stimulus", "delta", "delivery delta"]
    value = ["lick"]
    data_mean = analysis_functions.mean_bin(data, index, value, keep)

    # threshold trials to 200ms before to 2000ms after air puff
    data_mean = analysis_functions.thresh(data_mean, -200, 2000)
    
    # calculate performance for each animal for each time bin
    index = ["condition", "animal", "delta", "delivery delta"]
    keep = ["age", "sex", "strain", "Time (hr)", "Time (ms)"]
    perf = analysis_functions.performance(data_mean, index, keep)
    return (data_mean, counts, perf)

def generate_trialcount_plot(counts_data):
    '''Plot average number of trials and return the figure.
    
    Creates an individual plot for each condition. Plots average number of trials per bin vs time bins in training time.'''
    g = sns.catplot(counts_data, x="Time (hr)", y="trial no", col="condition", kind='bar', color="grey")
    
    for ax in g.axes.flat:
        ax.set_ylabel("Number of Trials")
        ax.set_ylim([0, 200])
    return g.fig