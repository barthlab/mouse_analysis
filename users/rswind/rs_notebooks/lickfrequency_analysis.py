# Analysis of (formatted) sensory association training data
# Last updated: 2024-29-05
# Author: Rachel Swindell

import pandas as pd
import numpy as np
import argparse

################################
#    SHORT HELPER FUNCTIONS    #
################################

def time_to_float(data, name, key, frequency):
    '''Return dataframe with a new column representing `key` column as a float.
    
    Parameters:
    data --- dataframe containing trial data. Must contain `key`.
    name --- name of new column
    key --- name of timedelta64 column to convert
    frequency --- string frequency at which to convert timedelta (eg. ms, s, h)'''
    data[name] = data[key].to_numpy(dtype=f"timedelta64[{frequency}]").astype("float")
    return data

def puff_delta(data, index, key_time, key_delay, name):
    '''Returns dataframe with new column of the given name with caluculated 
    time from each sample in a trial to airpuff delivery in that trial.
    
    Assumes SAT trial structure: IR beam break, then random delay between 200 
    and 800 ms, then 500ms of air puff delivery, then 500ms fixed delay, then 
    water delivery and timeout for 2000 ms or until animal stops licking 
    before next trial. All samples during the random delay period have a non-zero delay. 
    The delay value is taken from the first sample in a trial, which is 
    guaranteed to be during this random delay period. Start of air delivery
    is calculated as the time of the first sample in a trial plus the length of 
    the random delay. This value is subtracted from all samples in the trial to 
    get time from each sample to airpuff delivery.

    Parameters:
    data --- data frame containing trial data. Must contain `key_time`, `key_delay`,
             and all columns in `index`.
    index --- list of columns to group data by
    key_time --- name of column containg sample time
    key_delay --- name of column containing delay time for that trial
    name --- name of column to write delta to
    '''

    data = data.set_index(index)
    grouped = data.groupby(index)
    data["air_start"] = grouped[key_time].first() + grouped[key_delay].first()
    data[name] = data[key_time] - data["air_start"]
    data = data.drop(columns="air_start")

    return data.reset_index()

def get_first_sample(data, index, key, name):
    '''Returns dataframe with new column that has first row of grouped data 
    iterated across all rows in the group. 

    Used to copy timestamp at trial initation across all samples in a trial.

    Parameters:
    data --- data frame containing trial data. Must contain `key` and all 
             columns in `index`.
    index --- list of columns to group data by
    key --- column to take start time value from
    name --- name of column to write first sample to'''
    
    data = data.set_index(index)
    data[name] = data.groupby(index)[key].first()
    data = data.reset_index()
    return data

def fixed_window_lickfreq_helper(data,  index, values_mean, values_first, puff, start, end):
    '''Filter data to a fixed window of trial and return lick frequency in that
    window. Lick frequency (in Hz) is calculated as the total number of licks in
    the fixed window divided by the length of the window.

    Parameters:
    data --- data frame containing trial data. Must contain all columns in
             `values`, `keep`, and `index`
    values_mean --- list of columns to calculate lick frequeny
    values_first --- list of columns to keep
    index --- list of columns to index and group by
    puff --- name of column containing time to airpuff delivery
    freq_window --- length of window to calculate rolling window with in ms
    start --- start of fixed window relative to airpuff delivery in ms
    end --- end of fixed window relative to airpuff delivery in ms
    '''
    freq_window = end - start
    data = data[data[puff] >= pd.to_timedelta(start, unit="ms")]
    data = data[data[puff] <= pd.to_timedelta(end, unit="ms")]
    group = data.groupby(index)
    licks = group[values_mean].sum()/(freq_window/1000)
    keep = group[values_first].first()
    data_freq = pd.concat([licks, keep], axis=1).reset_index()
    return data_freq

def rolling_lickfreq(data, freq_window, index, values_mean, values_first, key):
    '''Calculate licking frequency across index groups as a time-based rolling 
    average.

    Mutually exclusive with resample_align. Window a fixed length in ms and is open on the 
    left and closed on the right. Licking frequency (in Hz) is the total number
    of licks in a window divided by the length of the window.

    Parameters:
    data --- data frame containing trial data. Must contain all columns in
             `values_mean`, `values_first` `keep`, `index`, and `key`
    freq_window --- length of window to calculate rolling window with in ms
    values_mean --- list of columns to calculate rolling average
    values_first --- list of columns to keep but not calculate rolling average
    index --- list of columns to index/group by
    key --- column to calculate rolling window with (datetime)
    '''
    # set timestamp column to index because pandas only allows rolling window 
    # calculations over index
    group = data.set_index(key).groupby(index)
    # licking frequency in Hz = number of licks/time window 
    licks = group.rolling(window=pd.to_timedelta(freq_window, unit="ms"), closed="right")[values_mean].sum()/(freq_window/1000)
    keep = data.set_index(index)[values_first].set_index(key,append=True)
    data_roll = pd.concat([licks, keep], axis=1).reset_index()

    return data_roll

def bin_by_time(data, bin_size, bin_unit, index, values_agg, values_first, key, 
                origin="start_day",offset=None, fn="mean"):
    '''Groups trials by time in given size and unit and apply provided function, and returns data with key 
    column replaced with time bin.

    If `fn` is "lickfreq", calculates lickfreq (number of licks/timebin size) for
    each timebin. If fn is `mean`, `std`, `var`, or `se`, calculates the respective
    measure. All other calculations must be a user-provided function.

    Parameters:
    data --- data frame containing trial data. Must contain all columns in 
             `index`, `values_mean`, `values_first` and `key`.
    bin_size --- size of bin to group by, in bin_unit (int)
    bin_unit --- unit of bin to group by (e.g. ms, min, hr)
    index --- list of columns to group data by
    values_agg --- set of columns to aggregate values for
    values_first --- set of columns to keep
    key --- column to group with, replaced with time bin
    origin --- start point to calculate timebins from
    offset --- amount to offset timebins from origin'''

    gp = pd.Grouper(key=key, freq=pd.to_timedelta(bin_size, unit=bin_unit), 
                    offset=offset, origin=origin)
    # append resampler to index group
    index = index.copy()
    index.append(gp)
    group = data.groupby(index)

    # calculate lick frequency
    if fn == 'lickfreq':
        licks = group[values_agg].sum()/(bin_size/1000)
    # handle some other common aggregations that are faster with builtins
    elif fn == 'mean':
        licks = group[values_agg].mean()
    elif fn == 'std':
        licks = group[values_agg].std()
    elif fn == 'var':
        licks = group[values_agg].var()
    elif fn == 'se':
        licks = group[values_agg].se()
    # all other calcualtions
    else:
        licks = group[values_agg].apply(fn)    

    keep = group[values_first].first()
    data_sample = pd.concat([licks, keep], axis=1).reset_index()

    return data_sample

def delta(data, index, key, key_acc, key_offset, name):
    '''Calculates time from start of SAT to trial and returns data in given 
    column.

    For the given group, subtracts the values of key_time and key_acc
    in the first row of the group from all values of key_time in the group.
    If used to calculate time from start of SAT, must provide the timebin each
    row is in as the key.

    Parameters:
    data --- data frame containing trial data
    index --- list of columns to index data
    gp --- column or list of columns to group by
    key --- column with bin times to use for calculation
    name --- name of column to write delta to'''

    data = data.set_index(index)
    delta = data[key] - (data[key_offset] + data[key_acc])
    data[name] = delta
    return data.reset_index()


def mean_bin(data, index, value, keep):
    '''Return mean of columns in `value` grouped by `index`, with columns in 
    `keep` retained.
    
    Parameters:
    data --- data frame containing trial data. Must contain columns in `index`,
             `value`, and `keep`.
    index --- list of columns to group by
    value --- list of columns to average
    keep --- list of columns to retain
    '''

    data = data.groupby(index)
    lick = data[value].mean()
    keep = data[keep].first()

    data = pd.concat([lick, keep], axis=1).reset_index()
    return data

def thresh(data, start, end, key):
    '''Threshold instaneanous sample data to be within a given timeframe from
    onset of airpuff delivery.
    
    data --- data frame containing trial data. Must contain `key`.
    start --- start of threshold timeframe in ms
    end --- end of threshold timeframe in ms
    key --- name of column containing time to puff delivery
    '''

    s = pd.to_timedelta(abs(start),unit="ms")
    if start < 0:
        s = -s

    e = pd.to_timedelta(abs(end), unit="ms")
    if end < 0:
        e = -e
    
    data = data[(data[key] >= s) & (data[key] < e)]
    return data

def performance(data, index, keep, key_trial, key_value, cond0, cond1):
    '''Calculate performance (difference between trial types)
    and return as new data frame.

    Groups data by index and given trial type column, then computes 
    cond1 - cond0 for each group to get performance. Takes the first row from 
    cond1 group for kept columns.
    
    Parameters:
    data --- data frame containing trial data. Must contain all columns in 
             `index`, `keep`, `key_trial`, and `key_value`.
    index --- list of columns to group data by
    keep --- list of columns to keep
    key_trial --- column containing trial type labels
    key_value --- column containing data to calculate performance
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
    
    Parameters:
    data --- data frame containing trial data. Must contain all columns in 
             `index`, `keep`, and `key`.
    index --- list of columns to group data by
    keep --- list of columns to keep
    key --- column to calculate trial counts from
    '''

    gp = data.groupby(index)
    counts = gp[key].nunique()
    keep = gp.first()[keep]

    data = pd.concat([counts, keep], axis=1).reset_index()
    return data

def get_bout_start_ind(gp, key, bt_ln):
    '''Get indecies representing the start of bouts, as defined by a given bout break length.
    a trial is considered the start of a bout if the end of the previous trial 
    ended longer than the given bout break length before.
    
    Parameters:
    gp ---  data to calculate bout start indecies for
    key --- name of column containing trial time
    bt_ln --- minimum break between two trials to classify them into separate bouts'''

    bout_start_ind = np.where(gp[key] - gp[key].shift() > bt_ln)[0]
    bout_start_ind = bout_start_ind + gp.index[0]
    bout_start_ind = np.insert(bout_start_ind, 0, gp.index[0])
    bout_start_ind = np.append(bout_start_ind, gp.index[-1])

    return bout_start_ind

def get_bouts(data, group, key, name, bout_len, bout_unit):
    '''Labels trials in bouts based on the given break lenth and unit.
    
    Two trials are placed in separate bouts if time between the end of one trial
    and the start of the next is more than `bout_len`.

    Parameters:
    data --- data frame containing trial data. <ust contain the columns in 
             `group` and `key`.
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
        # get the rows that are bout starts
        bout_starts = get_bout_start_ind(gp, key, bt_ln)

        gp[name] = pd.NA
        # label all rows between i-th index and i+1-th index as i
        for i in range(1, len(bout_starts)):
            gp.loc[bout_starts[i-1]:bout_starts[i],name] = i - 1

        gp[group] = gp_name
        res.append(gp)

    res = pd.concat(res)
    return res

def get_nth_part_trials(start, end, gp):
    '''Get set of trials from gp starting at `start` and ending at `end`. For 
    getting set of trials which represent nth fraction of a day.

    Parameters:
    gp --- group to get trials for
    start --- dataframe of starting trial no correpsonding to group names
    end --- datafram of ending trial no corresponding to group names
    '''
    if ~gp.empty:
        return gp[(gp["trial no"] > start.loc[gp.name]) & (gp["trial no"] <= end.loc[gp.name])]
    else:
        return []

def get_nth_part_day(data, frac, nth_part):
    '''Get set of trials representing the nth fraction of all trials an animal
    did in a day.
    
    Parameters:
    data --- data frame containing trial data
    frac --- denomenator of fraction of trials
    nth_part --- numenator of fraction of trials'''
# TODO: doesn't work correctly with string days -- FIXED?
    if nth_part > frac:
        print("Partition to return must be less than the number of partitions")
        return
    
    # total number of trials per day
    total_trials_day= data.sort_values(['condition', 'animal', 'trial no']).groupby(['condition', 'animal', 'Day']).last()['trial no']
    total_trials_day1 = total_trials_day.reset_index().groupby(['condition', 'animal']).first().reset_index().set_index(['condition', 'animal', 'Day'])['trial no']
    total_trials_day = total_trials_day.groupby(['condition', 'animal']).diff().fillna(total_trials_day1)
    
    # number of trials in a partition per day, rounded down
    trials_per_part = np.floor(total_trials_day*(1/frac))
    
    # first and last trial in partition for each day
    start = (trials_per_part*(nth_part-1))
    end = (trials_per_part*nth_part)
    
    # adjust first and last trial numbers for number of trials done on previous 
    # days
    for i in range(1, len(data["Day"].unique()) + 1):
        start = start + total_trials_day.groupby("animal").shift(i).fillna(0)
        end = end + total_trials_day.groupby("animal").shift(i).fillna(0)

    # convert to int to allow indexing
    start = start.astype(int)
    end = end.astype(int)
    
    # get trials in partition
    part = data.groupby(['condition', "animal", "Day"]).apply(lambda x: get_nth_part_trials(start, end, x))
    return part.reset_index(drop=True)

def day_to_label(day):
    '''Convert day as a float time to start of SAT representation to a string representation.'''
    day = int(day)
    if day < 0:
        return "ACC " + str(day + 1)
    return "SAT " + str(day + 1)

def sum_trials(data, group, keep, key):
    '''Get total number of trials grouped by given group. Only use on aggregated
    trial data.
    
    Parameters:
    data --- data frame containing trial data. Must contain all columns in
             `group`, `keep`, and `key`.
    group --- list of columns to group by
    keep --- list of columns to keep
    key --- column containing trial counts'''
    tot = data.groupby(group)[key].sum().rename('trial no')
    keep = data.groupby(group).first()[keep]
    res = pd.concat([tot, keep], axis=1).reset_index()
    return res



#################################
#    MAIN ANALYSIS FUNCTIONS    #
#################################

def align_to_timebin(data, time_bin, key_dict, keep, index, values):
    '''Align and bin trial data to timebin. Labels each sample with a 
    timebin as a timedelta and float, but does not perform any
    aggregation.
    
    Parameters:
    data --- data frame containing trial data.
    time_bin --- length of bin in minutes to group data by
    key_dict --- names of relevant columns in data
    keep --- list of columns to keep
    index --- list of columns to group by
    values --- list of columns that will be aggregated
    '''
    time = key_dict["timestamp"]
    acc = key_dict["acc"]
    millisbin = key_dict['millis bin']
    offset = key_dict['offset']

    key_start = "trial start"
    key_delta = "delta"

    # label samples with timebin as timestamp
    data = get_first_sample(data, index, time, key_start)
    if millisbin != None:
        index.append(millisbin)
        keep.remove(millisbin)
    data = bin_by_time(data, time_bin, "min", index, [], keep + [time, offset] + values, key_start,
                       offset=pd.to_timedelta(12, unit="h"))
    if millisbin != None:
        index.pop()
        keep.append(millisbin)

    # calculate timebin from start of sensory training
    data = delta(data, index, key_start, acc, offset, key_delta)

    data = time_to_float(data, "Time (hr)", key_delta, "m")
    data["Time (hr)"] = data["Time (hr)"]/60.
    return data

def align_to_day(data, key_dict, keep, index, values):
    '''Align and bin trial data to day. Labels each sample with a 
    day as a timedelta and float, but does not perform any
    aggregation.
    
    Parameters:
    data --- data frame containing trial data.
    time_bin --- length of bin in minutes to group data by
    key_dict --- names of relevant columns in data
    keep --- list of columns to keep
    index --- list of columns to group by
    values --- list of columns that will be aggregated
    '''
    time = key_dict["timestamp"]
    acc = key_dict["acc"]
    millisbin = key_dict['millis bin']
    offset = key_dict['offset']

    key_start = "trial start"
    key_delta = "delta"
    key_day = "day_delta"

    # label samples with day as timestamp
    data = get_first_sample(data, index, time, key_start)
    if millisbin != None:
        index.append(millisbin)
        keep.remove(millisbin)
    data = bin_by_time(data, 1, "day", index, [], 
                       keep + values + [time, offset, key_delta, "Time (hr)"], key_start, 
                       offset=pd.to_timedelta(12, unit="h"))
    if millisbin != None:
        index.pop()
        keep.append(millisbin)

    # calculate day from start of sensory training
    data = delta(data, index, key_start, acc, offset, key_day)

    data = time_to_float(data, "Day", "day_delta", "D")
    
    data = data.drop(columns=key_start)
    return data

def fixed_window_lickfreq(data, window, time_bin, key_dict, keep, values):
    '''Get lick frequency for a given window in the trial. Window is open on the
    right and clsoed on the left.
    
    Usually anticipatory lick frequency (700-1000 ms after airpuff delivery).
    
    Parameters:
    data --- data to be analyzed
    window --- tuple of (start, end) representing start and end of analysis window
    time_bin --- time bin in minutes to group trials into
    key_dict --- dictionary of relevant column names
    keep --- list of columns to keep
    values --- list of columns to calculate frequency for'''

    puff = key_dict["millis bin"]
    time = key_dict["timestamp"]
    index_trial = ["condition", "animal","trial no"]
    start, end = window
    key_start = "first sample"

    df = data.copy()
    data = puff_delta(df, index_trial, key_dict["timestamp"], key_dict["delay"], puff)

    # using first sample of each trial to align to timebin and day aligns trial
    # to bin where trial started in edge case where trial spans bins and 
    # anticipatory period happens after bin end
    # matches Matlab v16 analysis behavior
    data = get_first_sample(data, index_trial, time, key_start)
    key_dict['timestamp'] = 'first sample'

    # get lick frequency for each trial in given window
    data = fixed_window_lickfreq_helper(data, index_trial, values, 
                            keep + [time, key_start, key_dict['offset']],puff, start, end)
    
    # align trials to timebin and day
    data = align_to_timebin(data, time_bin, key_dict, keep, index_trial, values)
    data = align_to_day(data, key_dict, keep, index_trial, values)
    key_dict['timestamp'] = time

    return data

def instentanteous_lickfreq(data, freq_window, freq_bin, time_bin, key_dict, keep, values):
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
    key_dict --- dictionay of relevant column names
    keep --- 
    '''
    puff = key_dict["millis bin"]
    time = key_dict["timestamp"]
    index_trial = ["condition", "animal","trial no"]
    offset = key_dict['offset']
    
    df = data.copy()
    # get time to airpuff delivery
    data = puff_delta(df, index_trial, time, key_dict["delay"], puff)
    
    # correct for pandas behavior when resampling with negative timedeltas
    # if the minimum delay (800 ms) is not present in the dataset
    l = len(data.index)
    if (data['puff delta'].min()) != pd.to_timedelta('-800 ms'):
        data.loc[l, 'trial no'] = 0
        data.loc[l, 'puff delta'] = pd.to_timedelta('-800ms')
    
    # calculate instentanteous or rolling window lick frequency for each trial
    if freq_window != None:
        # rolling window
        keep.append(puff)
        data = rolling_time(data, freq_window, index_trial, values, keep + offset, 
                                         time)
        keep.pop()
    else:
        # discrete bins
        data = bin_by_time(data, freq_bin, "ms", index_trial, values, keep + [time, offset], puff, origin=None)
    
    # correct for pandas behavior when resampling with negative timedeltas
    # if the minimum delay (800 ms) is not present in the dataset
    data = data.loc[data['trial no'] > 0]
    
    # align sample times for cross-trial averaging
    data = bin_by_time(data, 100, "ms", index_trial, [], keep + [time, offset] + values, puff, origin=None)
    data = time_to_float(data, "Time (ms)", puff, "ms")
    # align to timebin and day
    data = align_to_timebin(data, time_bin, key_dict, keep + [puff,"Time (ms)"], index_trial, values)
    data = align_to_day(data, key_dict, keep + [puff,"Time (ms)"], index_trial, values)

    return data

def drop_group(data, mins, key_dict, key, index, keep):
    '''Drops bins with fewer than the given number of trials and returns data 
    without those bins.

    Filters bins based on total number of samples in index-based groups, 
    as well as number of samples from given key in a group. This key must be 
    binary (e.g. stimulus or blank; water or no water) for the key specific 
    functionality to work. That is, cond0 and cond1 must be from the same
    column, mutually exclusive, and have no other possible values.

    Parameters:
    data --- data frame containing trial data. Must contain columns in `index`,
             `col`, and `key`.
    mins --- dictionary containing
        `min_trials` --- minimum total number of trials to keep a bin (exclusive)
        `min_stimulus` --- minimum number of stimulus trials
        `min_blank` --- minimum number of blank trials
    key --- name of column containing trial number
    index --- list of columns to group data by
    keep --- list of columns to keep
    '''
    kp = keep.copy()
    # include time bin for aggregate grouping    
    index_timebin = index + [key_dict["time bin"]]
    if key_dict["time bin"] in kp:
        kp.remove(key_dict["time bin"])

    # drop bins with fewer than the given number of total, water, or blank trials
    stimulus = key_dict["stimulus type"]
    cond0, cond1 = np.sort(data[stimulus].unique())
    group = data.groupby(index_timebin + [stimulus])
    data = data.set_index(index_timebin)

    #filter bins with fewer than min_blank blank or min_water water trials
    counts = group[key].nunique().unstack(level=stimulus).fillna(0)
    cond = (counts.loc[:, cond0] >= mins["min_blank"]) & (counts.loc[:, cond1] >= mins["min_stimulus"]) 
    data = data[cond]

    #filter bins with fewer than min_trials total trials
    total_group = data.groupby(index_timebin) 
    total_counts = total_group[key].nunique()
    data = data[total_counts >= mins["min_trials"]]
    data = data.reset_index()
    return data

def agg(data, key_dict, index, keep, values):
    '''Aggregate lick frequency values across timebin and day.
    
    Parameters:
    data --- '''
    kp = keep.copy()
    millisbin = key_dict["millis bin"]
    millis = key_dict["Time (ms)"]
    stimulus = key_dict["stimulus type"]
    water = key_dict["water delivery"]

    index_timebin = index + [key_dict["time bin"]]
    if key_dict["time bin"] in kp:
        kp.remove(key_dict["time bin"])

    # number of trials for each timebin by stimulus type
    index_groups = index_timebin + [stimulus, water]
    cond0, cond1 = np.sort(data[stimulus].unique())
    if stimulus in kp:
        kp.remove(stimulus)
    if water in kp:
        kp.remove(water)
    
    counts_groups = trial_counts(data, index_groups, kp, key_dict["trial no"])
    
    # mean licking frequencies
    if millisbin != None:
        index_timebin = index_timebin + [millisbin]
        if millisbin in kp:
            kp.remove(millisbin)
            
    data_mean = mean_bin(data, index_timebin + [stimulus], values, kp + [water])

    # performance (Lstim - Lblank)
    perf = performance(data_mean, index_timebin, kp, stimulus, key_dict["licks"], 
                       cond0, cond1)
    #poke_perf = performance(data_mean, index_timebin, kp, stimulus, 'pokestart', 
    #                   cond0, cond1)

    return (data_mean, counts_groups, perf)#, poke_perf)

def write_analysis(out_files, new_data, columns, add=False):
    mode = 'a' if add else 'w'
    
    for k in range(len(new_data)):
        if add:
            with open(out_files[k]) as f:
                old_cols = f.readline().strip('\n').split(',')
            out = new_data[k][old_cols]
        else:
            out = new_data[k][columns[k]]
        out.to_csv(out_files[k], index=False, mode=mode, header=(not add), columns=columns[k])
        
def run_analysis(df, kind, keys, mins, keep, values, index, analysis_args, out_files=None, add=False):
    all_cols = index + keep + ["Day", "Time (hr)"]
    if kind =='full':
        all_cols = all_cols + ["Time (ms)"]
        data = instentanteous_lickfreq(df, *analysis_args, keys, keep, values) #freq_window, freq_bin, time_bin
        keys["time bin"] = "delta"
    elif kind == 'bin':
        data = fixed_window_lickfreq(df, *analysis_args, keys, keep, values) #r, time_bin
        keys["time bin"] = "delta"
    elif kind == 'nth_part':
        data = get_nth_part_day(df, *analysis_args) #num_parts, nth_part
        keys["time bin"] = "day_delta"
    else:
        print(f'Unknown analysis type {kind}.')
    
    keep = list(data.columns) 
    for c in index + values + ['trial no']:
        keep.remove(c)
    #keep + ["Time (hr)", "delta", "Day", "day_delta", "puff delta", "Time (ms)"]
    means, counts, perf = agg(data, keys, index, keep, values)
    
    data_cols = all_cols + ['trial no', 'delta', 'day_delta'] + values 
    means_cols = all_cols + values
    perf_cols = means_cols.copy()
    for i in ['stimulus', 'water', 'type']:
        if i in perf_cols:
            perf_cols.remove(i)
    trial_cols = all_cols + ['trial no']

    if kind =='nth_part':
        if out_files:
            columns = [data_cols, means_cols, perf_cols, trial_cols]
            write_analysis(out_files,[data, means, perf, counts], columns, add=add)
        return (means, counts, perf)
    
    daykeep = list(counts.columns)
    for c in index + ['Day', 'trial no', 'timestamp', 'offset']:
        if c in daykeep:
            daykeep.remove(c)

    daytots = sum_trials(counts, index + ["Day"], daykeep, "trial no")
    data_filtered = drop_group(data, mins, keys, "trial no", index, keep)
    filtmeans, filtcounts, filtperf = agg(data_filtered, keys, index, keep, values)
    
    if kind == 'full':
        columns = [data_cols,means_cols, means_cols, perf_cols, perf_cols]        
        output = (data, means, filtmeans, perf,  filtperf)
    elif kind=='bin':
        columns = [data_cols,means_cols, means_cols, perf_cols, perf_cols, trial_cols, trial_cols]        
        output = (data, means, filtmeans, perf,  filtperf, counts, daytots)
    if out_files: 
        write_analysis(out_files,output, columns, add=add)
    return output

if __name__ =='__main__':
    parser = argparse.ArgumentParser(
        description=
        '''Analyze lick frequency data. Must be pre-formatted output from format_arduino_data.py'''
    )
    parser.add_argument('files', nargs='+',help='Input formatted files.')
    parser.add_argument('-t', '--type', choices=['full', 'bin', 'nth_part'], help="Analysis type to run. 'full' analyzes all data, 'bin' selects a given bin, and 'nth_part' selects a given set of trials by day.")
    parser.add_argument('-g', '--args', nargs='+',help="Argument options for analysis. 'full' - freq_window, freq_bin, and time_bin. 'bin' - time_bin and analysis_window. 'nth_part' - num_parts and nth_part.")
    parser.add_argument('-a', '--add', action='store_true',help='If provided, append additional animals to output file. Requires only one input file is provided. Default behavior is to overwrite data in output file.')
    parser.add_argument('-o', metavar='OUT', nargs='+', help='Set of files to write output to. Requires only one input file is provided. Must be provided in this order: (data, means, counts, perf, daytots, data_filtered, filtmeans, filtcounts, filtperf, filtdaytots).')
    parser.add_argument('-m', '--drop_minimums', nargs=3, help='Minimum number of total, stimulus, and blank trials required to not drop a bin, in that order.')
    parser.add_argument('-k', '--columns', nargs='+', help='Set of columns to keep during analysis.')
    parser.add_argument('-v', '--values', nargs='+',help='Set of columns to analyze.')
    parser.add_argument('-i','--index', nargs='+', help='Set of columns to group by for analysis.')
    args = parser.parse_args()
    
    keys = {
        "millis bin":"puff delta",
        "Time (ms)":"Time (ms)",
        "stimulus type":"stimulus",
        "water delivery":"water",
        "trial no":"trial no",
        "time bin":"delta",
        "licks":"lick",
        "timestamp":"timestamp",
        "delay":"delay",
        "acc":"acc",
        'offset':'offset'
    }
    if args.type != 'full':
        keys['millis bin'] = None

    mins = {
        "min_trials": args.drop_minimums[0], 
        "min_blank": args.drop_minimums[2], 
        "min_stimulus": args.drop_minimums[1]
        }

    for fil in args.files:
        df = pd.read_csv(fil)
        run_analysis(df, type, keys, mins, args.columns,args.values, args.index, args.args, args.o, args.add)

