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
