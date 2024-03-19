# Functions for trace analysis
# Last updated: 09/25/2023 (refactor code)
# Author: Rachel Swindell

import pandas as pd
import numpy as np

################################
#    SHORT HELPER FUNCTIONS    #
################################

# helper functions for lick frequency

def time_to_float(data, name, key, frequency):
    data[name] = data[key].to_numpy(dtype=f"timedelta64[{frequency}]").astype("float")
    return data

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

def one_bin_lickfreq(data, freq_window, index, values_mean, values_first):
    '''Calculate licking frequency for one bin. Data must be prefilitered to
    the length of freq_window to get accurate results.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns 
             in values, keep, and index
    freq_window --- length of window to calculate rolling window with in ms
    values_mean --- list of columns to calculate rolling average
    values_first --- list of columns to keep but not calculate rolling average
    index --- list of columns to index/group by
    '''
    group = data.groupby(index)
    licks = group[values_mean].sum()/(freq_window/1000)
    keep = group[values_first].first()
    data_freq = pd.concat([licks, keep], axis=1).reset_index()
    return data_freq

def rolling_time(data, freq_window, index, values_mean, values_first, key):
    '''Calculate licking frequency across index groups as a time-based rolling 
       average and return as a dataframe.

    Mutually exclusive with resample_align. Window is in ms and is open on the 
    left and closed on the right.

    Parameters:
    data --- data frame containing trial data. Must contain all of the columns 
             in values, keep, and index
    freq_window --- length of window to calculate rolling window with in ms
    values_mean --- list of columns to calculate rolling average
    values_first --- list of columns to keep but not calculate rolling average
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

def bin_by_time(data, bin_size, bin_unit, index, values_mean, values_first, key, 
                origin="start_day",offset=None):
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

    gp = pd.Grouper(key=key, freq=pd.to_timedelta(bin_size, unit=bin_unit), 
                    offset=offset, origin=origin)
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

# helper functions for timebin aggregation

def drop_group_helper(data, min_trials, min0, min1, index, col, cond0, cond1, key):
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
    counts = group.count().unstack(level=key).fillna(0)
    cond = (counts.loc[:, (col, cond0)] >= min0) & (counts.loc[:, (col, cond1)] >= min1) 
    data = data[cond]

    #filter bins with fewer than min_trials total trials
    total_group = data.groupby(index) 
    total_counts = total_group.count()
    data = data[total_counts[col] >= min_trials]
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
    counts = gp[key].nunique()
    keep = gp.first()[keep]

    data = pd.concat([counts, keep], axis=1).reset_index()
    return data

def get_bout_start_ind(gp, key, bt_ln):
    bout_start_ind = np.where(gp[key] - gp[key].shift() > bt_ln)[0]
    bout_start_ind = bout_start_ind + gp.index[0]
    bout_start_ind = np.insert(bout_start_ind, 0, gp.index[0])
    bout_start_ind = np.append(bout_start_ind, gp.index[-1])

    return bout_start_ind


#################################################
#    SMALLER FUNCTIONS FOR OPTIONAL ANALYSES    #
#################################################

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

def set_prev_info(data, n_back, back_type):
    mean_statistics, counts, perf = data
    mean_statistics["n_back"] = n_back
    counts["n_back"] = n_back 
    perf["n_back"] = n_back
    mean_statistics["back_type"] = back_type
    counts["back_type"] = back_type
    perf["back_type"] = back_type
    return mean_statistics, counts, perf
    
# untested for n > 1
def agg_by_prev_trial(data, agg_data, n, min_trials, min_blank_trials, min_water_trials, keep):
    mean_statistics, counts, perf = agg_data
    mean_statistics["n_back"] = 0
    counts["n_back"] = 0
    perf["n_back"] = 0

    mean_statistics["back_type"] = "none"
    counts["back_type"] = "none"
    perf["back_type"] = "none"

    for i in range(1, n+1):
        cname = f'{i}_back'
        stimcond = data[cname] == "stimulus"
        bcond = data[cname] == "blank"
        for j in range(1, i):
            cname = f'{j}_back'
            stimcond = stimcond & (data[cname] == "stimulus")
            bcond = bcond & (data[cname] == "blank")
        prev_stim = data[data[cname] == "stimulus"]
        prev_blank = data[data[cname] == "blank"]

        ps_data = aggregate_analysis(prev_stim, min_trials, min_blank_trials, 
                                     min_water_trials, keep)
        ps_data = set_prev_info(ps_data, i, "stimulus")
        ps_mean, ps_counts, ps_perf = ps_data
        

        pb_data = aggregate_analysis(prev_blank, min_trials, min_blank_trials, 
                                     min_water_trials, keep)
        pb_data = set_prev_info(pb_data, i, "blank")
        pb_mean, pb_counts, pb_perf = pb_data
        

        mean_statistics = pd.concat([mean_statistics, ps_mean, pb_mean])
        counts = pd.concat([counts, ps_counts, pb_counts])
        perf = pd.concat([perf, ps_perf, pb_perf])
        return mean_statistics, counts, perf

def get_trials(start, end, gp):
    if ~gp.empty:
        return gp[(gp["trial no"] > start.loc[gp.name]) & (gp["trial no"] <= end.loc[gp.name])]
    else:
        return []

def get_nth_part_day(data, frac, nth_part):
# TODO: doesn't work correctly with string days
    if nth_part > frac:
        print("Partition to return must be less than the number of partitions")
        return
    ind = ["condition", "animal", "Day"]
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
    # convert to int for indexing
    start = start.astype(int)
    end = end.astype(int)
    # get trials in partition
    part = data.groupby(['condition', "animal", "Day"]).apply(lambda x: get_trials(start, end, x))
    return part.reset_index(drop=True)

def day_to_label(min, day):
    day = int(day)
    if day < 0:
        return "ACC " + str(day + 1)
    return "SAT " + str(day + 1)

def sum_trials(data, group, keep, key):
    tot = data.groupby(group)[key].sum().rename('trial no')
    keep = data.groupby(group).first()[keep]
    res = pd.concat([tot, keep], axis=1).reset_index()
    return res



#################################
#    MAIN ANALYSIS FUNCTIONS    #
#################################

def align_to_timebin_and_day(data, time_bin, key_dict, keep, index, values):
    time = key_dict["timestamp"]
    acc = key_dict["acc"]
    millisbin = key_dict['millis bin']
    offset = key_dict['offset']

    key_start = "trial start"
    key_delta = "delta"
    key_day = "day_delta"

    # bin by time (minutes)
    data = get_first_sample(data, index, time, key_start)
    if millisbin != None:
        index.append(millisbin)
        keep.remove(millisbin)

    data = bin_by_time(data, time_bin, "min", index, [], keep + [time, offset] + values, key_start,
                       offset=pd.to_timedelta(12, unit="h"))
    if millisbin != None:
        index.pop()
        keep.append(millisbin)
    data = delta(data, index, key_start, acc, offset, key_delta)
    data = time_to_float(data, "Time (hr)", key_delta, "m")
    data["Time (hr)"] = data["Time (hr)"]/60.

    # bin by day
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
    data = delta(data, index, key_start, acc, offset, key_day)
    data = time_to_float(data, "Day", "day_delta", "D")
    #minday = int(data["Day"].min())
    #data["Day"] = data["Day"].apply(lambda x: day_to_label(minday, x))

    data = data.drop(columns=key_start)
    return data

def bin_lickfreq(data, freq_window, time_bin, key_dict, keep, values):
    puff = key_dict["millis bin"]
    time = key_dict["timestamp"]
    index_trial = ["condition", "animal","trial no"]
    ant_start, ant_end = freq_window
    key_start = "first sample"

    df = data.copy()
    data = puff_delta(df, index_trial, key_dict["timestamp"], key_dict["delay"], puff)

    # for binning trials based on trial initiation
    # important for edge case where trial spans bins and anticipatory 
    # period happens after bin end
    data = get_first_sample(data, index_trial, time, key_start)
    timestamp = key_dict['timestamp']
    key_dict['timestamp'] = 'first sample'

    # filter for anticipatory licking (300 ms before water)
    data = data[data[puff] >= pd.to_timedelta(ant_start, unit="ms")]
    data = data[data[puff] <= pd.to_timedelta(ant_end, unit="ms")]
    
    # number of anticipatory licks divided by 
    # length of anticipatory licking period    
    data = one_bin_lickfreq(data, ant_end - ant_start, index_trial, values, 
                            keep + [time, key_start, key_dict['offset']])
    data = align_to_timebin_and_day(data, time_bin, key_dict, keep, index_trial, values)
    key_dict['timestamp'] = timestamp

    return data

def full_lickfreq(data,freq_window, freq_bin, time_bin, key_dict, keep, values):
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
    puff = key_dict["millis bin"]
    time = key_dict["timestamp"]
    index_trial = ["condition", "animal","trial no"]
    offset = key_dict['offset']

    df = data.copy()
    data = puff_delta(df, index_trial, time, key_dict["delay"], puff)
    # lick frequency
    if freq_window != None:
        # rolling window
        keep.append(puff)
        data = rolling_time(data, freq_window, index_trial, values, keep + offset, 
                                         time)
        keep.pop()
    else:
        # discrete bins
        data = bin_by_time(data, freq_bin, "ms", index_trial, values, keep + [time, offset], puff)

    # resampling for cross-trial sample alignment
    data = bin_by_time(data, 100, "ms", index_trial, [], keep + [time, offset] + values, puff)
    data = time_to_float(data, "Time (ms)", puff, "ms")
    data = align_to_timebin_and_day(data, time_bin, key_dict, keep + [puff,"Time (ms)"], index_trial, values)

    return data

def adj_index(index, kp, key_dict):
    # include time bin for aggregate grouping    
    index_timebin = index + [key_dict["time bin"]]
    if key_dict["time bin"] in kp:
        kp.remove(key_dict["time bin"])

    # include millis bin for counting (number of records at any sample is number
    # of trials in a bin) - only for full trial analysis (i.e. more than one 
    # sample per trial)
    millisbin = key_dict["millis bin"]
    if millisbin != None:
        index_timebin = index_timebin + [millisbin]
        if millisbin in kp:
            kp.remove(millisbin)
    return index_timebin, kp

def drop_group(data, mins, key_dict, key, index, keep):
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

    # threshold trials to 200ms before to 2000ms after air puff
    if millisbin != None:
        data_mean = thresh(data_mean, -200, 2000, millisbin)

    # performance (Lstim - Lblank)
    perf = performance(data_mean, index_timebin, kp, stimulus, key_dict["licks"], 
                       cond0, cond1)
    #poke_perf = performance(data_mean, index_timebin, kp, stimulus, 'pokestart', 
    #                   cond0, cond1)

    return (data_mean, counts_groups, perf)#, poke_perf)