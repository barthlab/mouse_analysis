

import numpy as np
import pandas as pd

def delta(data, key):
    '''calculates time since start of SAT. key is column to use to calculate'''
    data = data.set_index(["condition", "animal", "trial no"])
    delta = data[key] - data.sort_values(key).groupby("animal")[key].first() - data.groupby("animal")["acc"].first()
    data["delta"] = delta#.astype("timedelta64[m]")
    return data.reset_index()

def deliverydelta(df):
    '''find start of air delivery (trial start + delay to puff)'''
    df = df.set_index(["condition", "animal", "trial no"])
    grouped = df.groupby(["condition", "animal", "trial no"])
    df["air start"] = grouped["timestamp"].first() + grouped["delay"].first()
    df["delivery delta"]  = df["timestamp"] - df["air start"]
    df = df.drop(columns="air start")
    return df.reset_index()

def trialdelta(df):
    '''finds number of trials since air was turned on'''
    trials = df.groupby(["animal", "trial no"]).first().reset_index().set_index("animal")["trial no"]
    start = pd.to_timedelta(0, unit="ms")
    td = (trials - df[(df["delta"] >= start)].groupby("animal")["trial no"].first()).rename("trial delta")
    td = pd.concat([trials, td], axis=1).set_index("trial no", append=True)
    df = df.set_index(["animal", "trial no"])
    df["trial delta"] = td
    return df.reset_index()

def freqCurve(data, freq_bin, keep):
    '''
    Returns licks resampled by given frequency, in milliseconds.
    '''
    #freq window for a frequency across entire trial
    #licks are 1 if lick an 0 if no lick, sum licks over freq window and divide by window size
    group = data.groupby(["condition", "animal", "trial no", "trial type", pd.Grouper(key="delivery delta", freq=pd.offsets.Milli(freq_bin))])
    timestamp = group["timestamp"].first()
    keep = group[keep].first()
    licks = group[["lick", "poke"]].sum()/(freq_bin/1000)

    return pd.concat([licks, timestamp, keep], axis=1).reset_index()

def trialBin(data, trial_bin):
    '''Returns data resampled by given trial bin (as a string)'''
    #start of trial across all samples in trial, for by-trial resampling
    data["trial time"] = data.groupby(["condition", "animal", "trial no"])["timestamp"].first()
    data["acc"] = data.groupby(["condition", "animal", "trial no"])["acc"].first()

    #resamples trial time into bins, leaving other columns/index alone (easiest to do this way to maintain alignment of columns)
    #in particular does not effect lick frequency or delivery delta calculations
    group = data.reset_index().groupby(["condition", "animal", "trial no", "trial type", "delivery delta", pd.Grouper(key="trial time", freq=trial_bin)])
    data = group.first()
    data = data.reset_index()
    return data

def drop_bins(data, min_blank, min_water):
    '''Drops trial bins with fewer than min_trials blank or water trials'''
    lk_group = data.groupby(["condition", "animal", "delivery delta", "trial time", "trial type"])
    counts = lk_group.count().unstack(level="trial type")
    cond = (counts.loc[:, ("trial no", "blank")] > min_blank) & (counts.loc[:, ("trial no", "water")] > min_water) 
    filtered = data[cond]
    return filtered.stack(level="trial type").reset_index()

def performance(licks):
    '''Calculates performance (water - blank) for each time bin (assumes binned data).'''
    acc = licks.groupby(["condition", "animal"]).first()["acc"]
    lk_group = licks.groupby(["condition", "animal", "delivery delta", "delta", "trial type"])
    #mean lick/poke frequency and trial counts for each bin
    licks = lk_group[["lick", "poke"]].mean().unstack(level="trial type")
    #calculate performance (water - blank) for licks
    licks.loc[:, ("lick", "perf")] = licks.loc[:, ("lick", "water")] - licks.loc[:,("lick", "blank")]
    #for pokes
    licks.loc[:, ("poke", "perf")] = licks.loc[:, ("poke", "water")] - licks.loc[:,("poke", "blank")]
    licks = licks.stack(level="trial type").reset_index().set_index(["condition", "animal"])
    licks["acc"] = acc
    return licks

def thresh(data, start, end, start_pos, end_pos):
    '''threshold licking samples to samples between start and end. start_pos and end_pos indicate whether the start and end should be taken as positive or negative'''
    s = pd.to_timedelta(start,unit="ms")
    e = pd.to_timedelta(end, unit="ms")
    
    if not start_pos: s = -s
    if not end_pos: e = -e

    data = data[(data["delivery delta"] >= s) & (data["delivery delta"] < e)]
    
    return data

def intervals(animals, keepcols):
    '''Calculates the time between the end of a trial and the start of the next. keepcols is a list of columns to keep'''
    grouped = animals.groupby(["condition", "animal", "trial no"])
    interval = (grouped["timestamp"].first().shift(-1) - grouped["timestamp"].last()).rename("interval")
    inter_intervals = pd.concat([grouped[keepcols].first(),interval], axis=1)
    
    #the last trial has no next trial, and so difference should be NaT
    mask = inter_intervals.reset_index().groupby("animal").last().loc[:, "condition":"trial no"].reset_index()
    mask = mask[["condition", "animal", "trial no"]]
    mask_idx = list(mask.itertuples(index=False, name=None))
    inter_intervals.loc[mask_idx, "interval"] = pd.NaT
    
    return inter_intervals.reset_index()

def tweakplot(g, ymin, ymax, xmin, xmax, xlab="", ylab="", legend="", title=""):
    g.set_xlabels(xlab)
    g.set_ylabels(ylab)
    g._legend.set_title(legend)
    g.fig.suptitle(title, y=1.02)
    for ax in g.axes.flat:
        ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
        ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)
        ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(xmin, xmax)