# Functions for trace analysis
# Last updated: 09/14/2023 (adding function descriptions)
# Author: Rachel Swindell

import pandas as pd

# TODO: make functions fully generic or specific to this analysis (ie requiring specific column names)
# which is easier to debug long-term?

def rolling_frequency_average(data, freq_window, values, keep, index):
    '''Calculate licking frequency across index groups as a time-based rolling average and return as a dataframe.

    Mutually exclusive with resample_align.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in values, keep, and index
    freq_window --- length of window to calculate rolling window with in ms (int)
    values --- list of columns to calculate rolling average
    keep --- list of columns to keep but not calculate rolling average
    index --- list of columns to index/group by
    '''
    group = data.set_index("timestamp").groupby(index)
    licks = group.rolling(window=pd.to_timedelta(freq_window, unit="ms"), closed="right")[values].sum()/(freq_window/1000)
    keep = data.set_index(index)[keep].set_index("timestamp",append=True)
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
    #append resampler to index group
    index.append(gp)
    group = data.groupby(index)
    keep = group[keep].first()
    data_sample = keep.reset_index()

    return data_sample

def get_trial_start(data, index, key):
    '''Creates a columns with first row of a given column row across all samples when grouped by a given index of columns.

        Used to copy timestamp at trial initation across all samples in a trial.

       Parameters:
       data --- data frame containing trial data. Must contain all of the columns in index and key
       index --- list of columns to group data by
       key --- column to take start time value from'''

    data = data.set_index(index)
    data["trial time"] = data.groupby(index)[key].first()
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
def delta(data, index, key):
    '''Calculates time from start of SAT to trial and returns data with new column 'delta'.
    
    bin_by_time should be run first to get trial start times

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and key
    index --- list of columns to group data by
    key --- column with bin times to use for calculation'''

    data = data.set_index(index)
    delta = data[key] - data.sort_values(key).groupby("animal")[key].first() - data.groupby("animal")["acc"].first()
    data["delta"] = delta
    return data.reset_index()

# TODO: make generic
def puff_delta(data, index):
    '''Caluculates time from sample in a trial to air delivery in that trial and returns data with new column 'delivery delta'.
    
    Data must have 'timestamp' column indicating sample recording time to calculate puff delta.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns in index and 'timestamp' column
    index --- list of columns to group data by
    '''

    data = data.set_index(index)
    grouped = data.groupby(index)
    data["air start"] = grouped["timestamp"].first() + grouped["delay"].first()
    data["delivery delta"]  = data["timestamp"] - data["air start"]
    data = data.drop(columns="air start")
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
    group = data.set_index(index).groupby("trial type")
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