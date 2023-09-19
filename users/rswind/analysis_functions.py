# Functions for trace analysis
# Last updated: 09/14/2023 (adding function descriptions)
# Author: Rachel Swindell

import pandas as pd

# TODO: make functions fully generic or specific to this analysis (ie requiring specific column names)
# which is easier to debug long-term? -> generic is much easier to change column names or conditions as necessary

# TODO: make generic
def puff_delta(data, index, name, key_time, key_delay):
    '''Caluculates time from sample in a trial to air delivery in that trial and returns data with new column 'delivery delta'.
    
    Assumes SAT trial structure - random delay between 200ms to 

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and 'timestamp' column
    index --- list of columns to group data by
    key_time --- name of column containg sample time
    key_delay --- name of column containing delay time for that trial
    name --- name of column to write delta to
    '''

    data = data.set_index(index)
    grouped = data.groupby(index)
    air_start = grouped[key_time].first() + grouped[key_delay].first()
    data[name] = data[key_time] - air_start # data["air start"]
    
    return data.reset_index()



def rolling_frequency_average(data, freq_window, values, keep, index, key):
    '''Calculate licking frequency across index groups as a time-based rolling average and return as a dataframe.

    Mutually exclusive with resample_align. Window is in ms and is open on the left and closed on the right.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in values, keep, and index
    freq_window --- length of window to calculate rolling window with in ms (int)
    values --- list of columns to calculate rolling average
    keep --- list of columns to keep but not calculate rolling average
    index --- list of columns to index/group by
    key --- column to calculate rolling window with
    '''
    # can only do rolling window calculations over index
    group = data.set_index(key).groupby(index)
    # licking frequency in Hz = number of licks/time window 
    licks = group.rolling(window=pd.to_timedelta(freq_window, unit="ms"), closed="right")[values].sum()/(freq_window/1000)
    keep = data.set_index(index)[keep].set_index(key,append=True)
    data_roll = pd.concat([licks, keep], axis=1).reset_index()

    return data_roll


# TODO: not currently usable for binning - doesnt ask for set of columns to calculate average on
# add values
def resample_align(data, frq, key, keep, index):
    '''Resample data across index values into the given time-based bins.
    
    Mutually exclusive with rolling_frequency_average. 

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in key, keep, and index
    frq --- bin size in milliseconds. 100ms (base sample frequency) aligns samples across trials without changing values
    key --- name of column to resample by
    keep --- list of columns to keep but not calculate average
    index --- list of columns to index/group by
   
    Resample at the given freq. If index is a set of columns that puts each sample into its own unique group, 
    a 100ms frequency aligns data without changing values'''

    gp = pd.Grouper(key=key, freq=pd.offsets.Milli(frq))
    # append resampler to index group
    index.append(gp)
    group = data.groupby(index)
    keep = group[keep].first()
    data_sample = keep.reset_index()

    return data_sample

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


def bin_by_time(data, bin_size, bin_unit, index, key):
    '''Groups trials by time in given size and unit, and returns data with key column replaced with time bin
    
    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and key
    bin_size --- size of bin to group by, in bin_unit (int)
    bin_unit --- unit of bin to group by (e.g. hr, min) (string compatible with pandas timedelta)
    index --- list of columns to group data by. In order to keep all rows, needs to be specific enough that each row is a unique group
    key --- column with trial starts to group with, replaced with time bin'''

    # TimeGrouper object on key
    gp = pd.Grouper(key=key, freq=pd.to_timedelta(bin_size, unit=bin_unit))
    # list (copy so input list is not destructively modified)
    index = index.copy()
    index.append(gp)
    # group original data by all columns specified and time groups generated with time grouper object
    group = data.groupby(index)
    # keeps first (and usually only) row in each group but writes time group to key
    data = group.first()
    data = data.reset_index()
    return data

# TODO: make fully generic
def delta(data, keep, index, key, name):
    '''Calculates time from start of SAT to trial and returns data in given column.
    
    bin_by_time should be run first to get trial start times

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and key
    keep --- list of columns to keep
    gp --- column or list of columns to group by
    key --- column with bin times to use for calculation
    name --- name of column to write delta to'''

    data = data.set_index(index)
    delta = data[key] - data.sort_values(key).groupby("animal")[key].first() - data.groupby("animal")["acc"].first()
    data[name] = delta
    return data.reset_index()


# TODO: make generic
def bin_prev_identity(data, key, index):
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
    puff_name = "puff delta"

    index_trial = ["condition", "animal","trial no"]
    index_trial_type = ["condition", "animal", "trial no", "stimulus"]

    keep = ["timestamp", "age", "sex", "strain", "acc"]


    # do not modify loaded data
    df = data.copy()

    # calculated time to air puff
    # index = ["condition", "animal","trial no"]
    data = puff_delta(df, index_trial, puff_name, key_time, key_delay)

    # rolling window average licking frequency
    # keep = ["timestamp", "age", "sex", "strain", "acc", "delivery delta"]
    # index = ["condition", "animal", "trial no", "stimulus"]
    keep.append(puff_name)
    values = ["lick", "poke"]
    # TODO: need stimulus ? 
    
    data = rolling_frequency_average(data, freq_window, values, keep, index_trial_type)
    keep.pop()

    # resampling for trial alignment
    # if discrete binning instead of rolling window is used, also bins samples across trial
    # keep = ["timestamp", "lick", "poke", "age", "sex", "strain", "acc"]
    keep_align = keep + values
    # TODO: need stimulus ?
    index = ["condition", "animal", "trial no", "stimulus"]
    data = resample_align(data, freq_bin, puff_name, keep_align, index)

    # label start of trial at each sample
    # index = ["condition", "animal", "trial no"]
    data = get_first_sample(data, index_trial, key_time)

    # bin by time
    index = ["condition", "animal", "trial no", "stimulus", "delivery delta"]
    data = bin_by_time(data, time_bin, "min", index, "trial time")

    # calulate time bin relative to start of SAT training
    # index = ["condition", "animal", "trial no"]
    name = "delta"
    key = "trial time"
    data = delta(data, index, "trial time")
    
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

    # convert time bins and trial time to float representations
    # works better in pre-analysis but cant be moved because of drop_bins functionality
    data["Time (hr)"] = data["delta"].to_numpy(dtype="timedelta64[m]").astype("float")/60.
    data["Time (ms)"] = data["delivery delta"].to_numpy(dtype="timedelta64[ms]").astype("float")

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