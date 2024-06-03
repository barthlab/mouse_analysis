# various unused and deleted codes


def generate_barplots(perf_data):
    '''Plot various performance measures.
    
    Unused plotting code that could be useful as different ways of visualizing data in the future. Returns void.'''

    # average performance trace by animal
    cond = (perf_data["Time (hr)"] < 24) & (perf_data["Time (hr)"] > -24) 
    g = sns.relplot(data=perf_data[cond],kind="line",x="Time (ms)", 
                    y="lick",col="animal", col_wrap=5,hue="Time (hr)", palette="coolwarm", errorbar=None,legend="full")
    for ax in g.axes.flat:
        ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
        ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
        ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)

    # average performance trace by timebin, 4h timebins, with each animals' performance in the background
    cond = (perf_data["Time (hr)"] < 24) & (perf_data["Time (hr)"] > -24) 
    g = sns.relplot(data=perf_data[cond],kind="line",x="Time (ms)", 
                    y="lick",col="Time (hr)", row="condition", errorbar=None, zorder=5)
    for (train_type, time), ax in g.axes_dict.items():
        an_cond = perf_data["condition"] == train_type
        time_cond = perf_data["Time (hr)"] == time
        sns.lineplot(
                data=perf_data[cond & an_cond & time_cond], x="Time (ms)", y="lick", hue="animal", 
                palette=["grey"], alpha=0.5,linewidth=1, ax=ax, legend=False)
        ax.axhline(y=0, xmin=0, xmax=1, ls="--", lw=0.75,color="black", zorder=0)
        ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
        ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)

    # average performance by timebin across all animals
    index = ["condition", "delta", "delivery delta"]
    keep = ["age", "sex", "strain", "Time (hr)", "Time (ms)"]
    gp = perf_data.groupby(index)
    perf_avg = gp["lick"].mean()
    keep = gp[keep].first()
    perf_avg = pd.concat([perf_avg, keep], axis=1).reset_index()

    # average time until performance becomes positive from start of air for each timebin
    cond = (perf_avg["lick"] > 0) & (perf_avg["Time (ms)"] > 100)
    pt = perf_avg[cond].groupby(["condition", "delta"]).first().reset_index()

    # plot a few bins during acclimation
    cond = (pt["Time (hr)"] < -4) & (pt["Time (hr)"] > -24)
    g = sns.catplot(data=pt[cond], kind="bar", x="condition", y="Time (ms)", col="Time (hr)", col_wrap=5)
    g.fig.suptitle("Time to Positive Performance", x=0.4, y=1.02)
    g.set_xlabels("Condition")
    g.set_ylabels("Time (ms)")    

    # plot a few bins during SAT1
    cond = (pt["Time (hr)"] < 20) & (pt["Time (hr)"] > 0)
    g = sns.catplot(data=pt[cond], kind="bar", x="condition", y="Time (ms)", col="Time (hr)", col_wrap=5)
    g.fig.suptitle("Time to Positive Performance", x=0.4, y=1.02)
    g.set_xlabels("Condition")
    g.set_ylabels("Time (ms)")    
    cond = (perf_avg["Time (hr)"] < 20) & (perf_avg["Time (hr)"] >= 0) & (perf_avg["Time (ms)"] > 0) & (perf_avg["Time (ms)"] < 1000)
    pt_puff = perf_avg[cond].sort_values("lick").groupby(["condition", "delta"]).first().reset_index()

    # average time from puff to minimum performance (quantifing dip in performance after air)
    g = sns.catplot(data=pt_puff, kind="bar", x="condition", y="Time (ms)", col="Time (hr)", col_wrap=5)
    g.fig.suptitle("Time to Minimum Performance", x=0.5, y=1.02)
    g.set_xlabels("Condition")
    g.set_ylabels("Time (ms)")    

    # average magnitude of minimum performance
    g = sns.catplot(data=pt_puff, kind="bar", x="condition", y="lick", col="Time (hr)", col_wrap=5)
    g.fig.suptitle("Magnitude of Minimum Performance", x=0.5, y=1.02)
    g.set_xlabels("Condition")
    g.set_ylabels("Performance")

#number of blank trials
pd.Series(NJO1["trial no"][NJO1["water"]==9].unique()).count()
#number of water trials
pd.Series(NJO1["trial no"][NJO1["water"]==5].unique()).count()
#plot number of trials over time
NJO1.set_index("timestamp")["trial no"].resample('4H').count().plot()

#gets anticipatory licking window with same method as matlab code
#works but is SLOW (2-4s on 2d of training)
def fltr (g):
  idx = SAT_animals["trial no"][g.index[0]]
  return g[(g["timestamp"] > antwindleft[idx]) & (g["timestamp"] < antwindright[idx])]

puffstart = NJO1_grouped.first()["timestamp"] + NJO1_grouped.first()["delay"]
antwindleft = (puffstart + pd.to_timedelta(700, unit="ms"))
antwindright = (puffstart + pd.to_timedelta(1, unit="s"))
antwindow = NJO1_grouped.apply(fltr) #SLOW
antwindow = antwindow.reset_index(drop=True)
antwindow

#split anticipatory licking frequency by trial type and plot by 4 hour bins
antwindow_group = antwindow[antwindow["lick"] > 0].groupby("trial no")
antlickfreq = pd.concat([antwindow_group["timestamp"].first(), antwindow_group["trial type"].first(), antwindow_group["lick"].count()/.3], axis=1)
water = antlickfreq[antlickfreq["trial type"] == "water"].set_index("timestamp").resample("4H").mean() 
blank = antlickfreq[antlickfreq["trial type"] == "blank"].set_index("timestamp").resample("4H").mean()
pd.concat([water, blank], axis=1).fillna(0).plot()


#transform with user-defined function is SLOW

NJO1_grouped = NJO1

def delt (x):
    x = x.reset_index()    
   # pos = x.columns.get_loc("timestamp")
    res = x.iloc[1:, 1] - x.iat[0, 1]
    res = pd.concat([pd.Series(pd.to_timedelta(0)), res])
    res.index = x["index"]
    return res

#NJO1_grouped["delta"] = NJO1_grouped.groupby("trial no")["timestamp"].transform(delt)

NJO1_grouped.loc[NJO1_grouped["lick"] == 2,"lick"] = 1

NJO1_grouped = NJO1_grouped.set_index("timestamp")
lickfreq = NJO1_grouped.groupby("trial no")["lick"].resample("300 ms").sum()/.3
#lck = lickfreq.reset_index()
#lck["delta"] = lck.groupby("trial no")["timestamp"].transform(delt)
#lck = lck.set_index(["trial no", "delta"])
lck["lick"].unstack(level=0).plot()
#NJO1_grouped
#water trials

#NJO1_go = NJO1_go.set_index(["timestamp", "trial no", "trial type"])

#grouped = NJO1_go[NJO1_go["water"]==3].groupby(["trial no"]).tail(10).groupby(["trial no"])

#pd.concat([grouped["timestamp"].tail(1).to_frame().reset_index(),grouped["lick"].mean().to_frame().reset_index()], axis=1).set_index("timestamp").resample("4H").mean()["lick"]






def freqCurve(data, freq_bin, trial_bin):
    '''
    Returns lick frequency resampled by given frequency and binned into given trial bin size.
    Arguments:
    data - data to resample
    freq_bin - size of bin to calculate lick and poke frequency with, in ms
    trial_bin - size of bin to group trials by, in H
    '''
    #find start of water delivery (trial start + delay to puff + 1s)
    grouped = data.groupby(["condition", "animal", "trial no"])
    data["water start"] = grouped["timestamp"].first() + grouped["delay"].first() + pd.to_timedelta("1s")
    data["delivery delta"]  = data["timestamp"] - data["water start"]
    data = data.reset_index().set_index("timestamp")
    #freq window for a frequency across entire trial
    #licks are 1 if lick an 0 if no lick, sum licks over freq window and divide by window size
    licks = data.groupby(["condition", "animal", "trial no", "trial type", pd.Grouper(key="delivery delta", freq=pd.offsets.Milli(freq_bin))])[["lick", "poke"]].sum()/(freq_bin/1000)

    #start of trial across all samples in trial, for by-trial resampling
    licks["trial time"] = data.reset_index().groupby(["condition", "animal", "trial no"])["timestamp"].first()

    #resamples trial time into bins, leaving other columns/index alone (easiest to do this way to maintain alignment of columns)
    #in particular does not effect lick frequency or delivery delta calculations
    licks = licks.reset_index().groupby(["condition", "animal", "trial no", "trial type", "delivery delta", pd.Grouper(key="trial time", freq=pd.offsets.Hour(trial_bin))]).first()
    licks = licks.reset_index().set_index("animal")
    return licks
l_2_4 = freqCurve(animals, 200, 4)

#get performance in each bin across entire trial
#drop bins with < 10 trials for either group (parameterize in function)
licks = l_2_4
#mean frequency across a trial time group for water and blank
lk_group = licks.groupby(["condition", "animal", "delivery delta", "trial time", "trial type"])[["lick", "poke"]]
#mean lick/poke frequency and trial counts for each bin
licks = lk_group.mean().unstack(level="trial type")
counts = lk_group.count().unstack(level="trial type")
#at least 10 water trials and 10 blank trials
cond = (counts.loc[:, ("trial no", "blank")] > 10) & (counts.loc[:, ("trial no", "water")] > 10) 
licks = licks[cond]
#calculate performance (water - blank) for licks
licks.loc[:, ("lick", "perf")] = licks.loc[:, ("lick", "water")] - licks.loc[:,("lick", "blank")]
#for pokes
licks.loc[:, ("poke", "perf")] = licks.loc[:, ("poke", "water")] - licks.loc[:,("poke", "blank")]
licks = licks.stack(level="trial type").reset_index().set_index("animal")
licks


# from trace_notebook
freq = analysis_helper.freqCurve(aa, 300)
freq
#resamples trial time into bins, leaving other columns/index alone (easiest to do this way to maintain alignment of columns)
#in particular does not effect lick frequency or delivery delta calculations
group = freq.groupby(["condition", "animal", "trial no", "trial type", "delivery delta", pd.Grouper(key="timestamp", freq="4h")])
freqbinned = group.first()
freqbinned= freqbinned.reset_index()
freqbinned["delta"] = freqbinned["timestamp"]
freqbinned["timestamp"] = freq["timestamp"]

freqbinned
perf = analysis_helper.thresh(perf, 2000, 3000, False, True)

perf["delta"] = perf["delta"] - perf.sort_values("delta").groupby("animal")["delta"].first() - perf.groupby("animal")["acc"].first()
perf["delta"] = perf["delta"].astype("timedelta64[h]")
perf["delivery delta"] = perf["delivery delta"].astype("timedelta64[ms]")
cond = (perf["delta"] < 24) & (perf["delta"] > -24)
g = sns.relplot(data=perf[(perf["trial type"] == "perf") & cond ],kind="line",x="delivery delta", 
                y="lick",col="condition", col_wrap=2,hue="delta", palette="coolwarm", errorbar=None,legend="full")
for ax in g.axes.flat:
    ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
    ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
    ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)

data=freqbinned.reset_index().set_index(["condition", "animal", "delivery delta", "delta"]) 
lk_group = data.groupby(["condition", "animal", "delivery delta", "delta", "trial type"])
counts = lk_group.count().unstack(level="trial type")
cond = (counts.loc[:, ("trial no", "blank")] > 10) & (counts.loc[:, ("trial no", "water")] > 10) 
filtered = data[cond].reset_index()


perf = analysis_helper.performance(filtered)
perf = analysis_helper.thresh(perf, 2000, 3000, False, True)

#perf = perf.set_index(["condition", "animal"])
perf["delta"] = perf["delta"] - perf.sort_values("delta").groupby("animal")["delta"].first() - perf.groupby("animal")["acc"].first()
perf["delta"] = perf["delta"].astype("timedelta64[h]")
perf["delivery delta"] = perf["delivery delta"].astype("timedelta64[ms]")
perf = perf.reset_index()

#filtered["delta"] = filtered["delta"] - filtered.sort_values("delta").groupby("animal")["delta"].first() - filtered.groupby("animal")["acc"].first()
#filtered["delta"] = filtered["delta"].astype("timedelta64[h]")
#filtered["delivery delta"] = filtered["delivery delta"].astype("timedelta64[ms]")
#filtered = filtered.reset_index()
perf


cond = (perf["delta"] < 24) & (perf["delta"] > -24)
g = sns.relplot(data=perf[(perf["trial type"] == "perf") & cond],kind="line",x="delivery delta", 
                y="lick",col="condition", col_wrap=2,hue="delta", palette="coolwarm", errorbar=None,legend="full")
for ax in g.axes.flat:
    ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
    ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
    ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)
#end = perf.groupby(["condition", "animal", "delta"]).max()["delivery delta"]
start = perf.groupby(["condition", "animal", "delta"]).min()#.reset_index()#.set_index(["condition", "animal"])
#perf = perf.set_index(["condition", "animal"])
flt = perf[(perf["delivery delta"] > 0)& (perf["delivery delta"] < 1000)]
flt[flt["trial type"] == "perf"].sort_values(["animal", "delta", "delivery delta"])
means = flt.groupby(["condition", "delta", "delivery delta", "trial type"]).mean().reset_index()
delta = flt.groupby(["condition", "delta", "trial type"])[["delivery delta", "lick"]].diff()
flt[["deltax", "deltay"]] = delta
flt["slope"] = flt["deltay"]/(flt["deltax"]/1000)
slope = flt.groupby(["condition", "trial type", "delta"])["slope"].mean().reset_index()
#means[means["delta"] == 16]
slope
cond = (slope["delta"] < 24) & (slope["delta"] > -24)
g = sns.catplot(data=slope[(slope["trial type"] == "perf") & cond],x="delta", 
                y="slope",col="condition", col_wrap=2,hue="delta", palette="coolwarm", kind="bar")
cond = (flt["delta"] < 24) & (flt["delta"] > -24)
g = sns.catplot(data=flt[(flt["trial type"] == "perf") & cond],x="delta", 
                y="slope",col="condition", col_wrap=2,hue="delta", palette="coolwarm", kind="box")


data = data_sample
trial_bin="1h"

data = data.set_index(["condition", "animal", "trial no"])
data["trial time"] = data.groupby(["condition", "animal", "trial no"])["timestamp"].first()
data = data.reset_index()

group = data.groupby(["condition", "animal", "trial no", "trial type", "delivery delta", pd.Grouper(key="trial time", freq=trial_bin)])
data = group.first()
data = data.reset_index()

data = analysis_helper.delta(data, key="trial time")

min_trials = 10
min_blank = 2
min_water = 2

data = data.set_index(["condition", "animal", "delta","delivery delta", "trial type"])
group = data.groupby(["condition", "animal", "delta","delivery delta", "trial type"])
counts = group.count().unstack(level="trial type")
cond = (counts.loc[:, ("trial no", "blank")] > min_blank) & (counts.loc[:, ("trial no", "water")] > min_water) 
data = data[cond]

data = data.reset_index().set_index(["condition", "animal", "delta","delivery delta"])
total_group = data.groupby(["condition", "animal", "delta","delivery delta"]) 
total_counts = total_group.count()

data = data[total_counts["trial no"] > min_trials]
data = data.reset_index()

keep=["age", "sex", "strain", "SAT_start", "acc"]

data_bin = data.groupby(["condition", "animal", "trial type", "delta", "delivery delta"])
lick = data_bin["lick"].mean()
keep = data_bin[keep].first()

data_bin = pd.concat([lick, keep], axis=1).reset_index()

data_thresh = analysis_helper.thresh(data_bin, 200, 3000, False, True)

group = data_thresh.set_index(["condition", "animal", "delta", "delivery delta"]).groupby("trial type")
water = group.get_group("water")["lick"]
blank = group.get_group("blank")["lick"]
data_perf = water - blank
data_perf = data_perf.reset_index()
data_perf["delta_float"] = data_perf["delta"].astype("timedelta64[h]")
data_perf["delivery delta_float"] = data_perf["delivery delta"].astype("timedelta64[ms]")
data_perf






def aggregate_analysis(data, min_trials, min_blank, min_stimulus, keep):
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
    min_blank --- minimum number of blank trials to keep a bin
    keep --- list of columns to add to the columns kept after analysis'''
    
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_acc = "acc"

    index_animal = ["condition", "animal"]
    index_timebin = index_animal + [key_delta, key_puff]

    keep = ["age", "sex", "strain", "acc", "cage", "Time (hr)", "Time (ms)"] + keep
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

def bin_aggregate_analysis(data, min_trials, min_blank, min_stimulus, keep):
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
    min_blank --- minimum number of blank trials to keep a bin
    keep --- list of columns to add to the columns kept after analysis'''
    
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_acc = "acc"

    index_animal = ["condition", "animal"]
    index_timebin = index_animal + [key_delta]

    keep = ["age", "sex", "strain", "acc","cage", "Time (hr)"] + keep
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
    
    index_groups = index_timebin + ["stimulus", "water"]
    counts_groups = trial_counts(data, index_groups, keep, "trial no")

    # mean licking frequencies
    data_mean = mean_bin(data, index_timebin + ["stimulus"], values, keep + ["water"])
    
    # performance (Lstim - Lblank)
    perf = performance(data_mean, index_timebin, keep, key, "lick", 
                       cond0, cond1)

    return (data_mean, counts_groups, perf)


def lickfreq_analysis(data, freq_window, freq_bin, time_bin):
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_day = "day_delta"
    key_acc = "acc"

    index_trial = ["condition", "animal","trial no"]

    keep = ["timestamp", "age", "sex", "strain", "acc", "cage", "stimulus", "water", "type"]
    values = ["lick", "poke"]


    # do not modify loaded data
    df = data.copy()

    # calculated time from sample to air puff
    data = puff_delta(df, index_trial, key_time, key_delay, key_puff)

    # lick frequency
    if freq_window != None:
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
    data = time_to_float(data, "Time (hr)", key_delta, "m")
    data["Time (hr)"] = data["Time (hr)"]/60.
    
    # bin by day
    data = get_first_sample(data, index_trial, key_time, key_start)
    index_trial.append(key_puff)

    data = bin_by_time(data, 1, "day", index_trial, [], 
                       keep + values + ["delta"], key_start, 
                       offset=pd.to_timedelta(12, unit="h"))
    index_trial.pop()
    data = delta(data, index_trial, "animal", key_start, key_acc, key_day)

    # convert time bins and trial time to float representations
    data = time_to_float(data, "Time (ms)", key_puff, "ms")



    return data

def bin_lickfreq_analysis(data, freq_window, time_bin):
    '''Computes licking frequency in the given licking window and
    returns raw licking frequency.
    
    Parameters provide flexibility in multiple ways for both cross-trial and 
    cross-training binning.

    Parameters:
    data --- data to be analyzed
    freq_window --- tuple indicating start and 
    time_bin --- size of bin for multi-trial binning over entire training 
                 (minutes)
    '''
    # parameters
    key_time = "timestamp"
    key_delay = "delay"
    key_puff = "puff delta"
    key_start = "trial start"
    key_delta = "delta"
    key_day = "day_delta"
    key_acc = "acc"

    index_trial = ["condition", "animal","trial no"]

    keep = ["timestamp", "age", "sex", "strain", "acc", "cage", "stimulus", "water", "type"]
    values = ["lick", "poke"]

    ant_start, ant_end = freq_window

    # do not modify loaded data
    df = data.copy()

    # calculated time from sample to air puff
    data = puff_delta(data, index_trial, key_time, key_delay, key_puff)

    # filter for anticipatory licking (300 ms before water)
    data = data[data[key_puff] >= pd.to_timedelta(ant_start, unit="ms")]
    data = data[data[key_puff] <= pd.to_timedelta(ant_end, unit="ms")]
    
    # number of anticipatory licks divided by 
    # length of anticipatory licking period
    data = one_bin_lickfreq(data, ant_end - ant_start, index_trial, values, keep)
    
    # label start of trial at each sample
    data = get_first_sample(data, index_trial, key_time, key_start)

    # bin by time (minutes)
    data = bin_by_time(data, time_bin, "min", index_trial, [], keep + values, key_start)

    # calulate time bin relative to start of SAT training
    data = delta(data, index_trial, "animal", key_start, key_acc, key_delta)

    # bin by day
    data = get_first_sample(data, index_trial, key_time, key_start)

    # label days
    data = bin_by_time(data, 1, "day", index_trial, [], 
                       keep + values + ["delta"], key_start, 
                       offset=pd.to_timedelta(12, unit="h"))
    data = delta(data, index_trial, "animal", key_start, key_acc, key_day)
    data = time_to_float(data, "Day", "day_delta", "D")
    minday = int(data["Day"].min())
    data["Day"] = data["Day"].apply(lambda x: day_to_label(minday, x))

    # convert time bins and trial time to float representations
    data = time_to_float(data, "Time (hr)", key_delta, "m")
    data["Time (hr)"] = data["Time (hr)"]/60.
    
    return data

