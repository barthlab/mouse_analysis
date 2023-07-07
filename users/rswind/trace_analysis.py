#!python3
# last edit 7/6/2023 Rachel Swindell

import re
import os
import argparse

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

pd.options.mode.chained_assignment = None  # default='warn'

# to run:
# python3 trace_analysis.py --metadata metadata-file --trial_bin_size bin-size --bin_size bin-size [--group] path-to-data

# read + interpret console commands
# adapted from code in mouse_analysis.py written by Alex Kiroff
def parse_args():
    #default_acclimation_time = 48  # hours
    default_bin_time = 4  # hours
    default_trial_bin = 300 # ms
    default_last_x_percent = 20 # percent

    parser = argparse.ArgumentParser(description="Advanced mouse behavior analysis")
    #parser.add_argument("--animal_name", required=True, type=str, help="Name of the animal, to name the output files")
    parser.add_argument("--trial_bin", required=True, type=float, default=default_trial_bin, help="Licking frequency bin size in ms (default: {default_trial_bin})")
    #parser.add_argument("--csv_files", required=True, type=str, nargs="+", help="In order paths to the CSV files")
    parser.add_argument("--metadata", required=True, type=str, help="Path to file with animal metadata")
    #parser.add_argument("--acclimation_time", type=float, default=default_acclimation_time, help=f"Acclimation time in hours (default: {default_acclimation_time})")
    parser.add_argument("--bin_time", type=float, default=default_bin_time, help=f"Bin time in hours (default: {default_bin_time})")
    parser.add_argument("--last_x_percent", type=float, default=default_last_x_percent, help=f"Percentage of data to consider (0 < last_x_percent < 100) (default: {default_last_x_percent})")

    parsed = parser.parse_args()

    # Check if all csv_files are valid files
    for csv_file in parsed.csv_files:
        if not os.path.isfile(csv_file):
            parser.error(f"The provided CSV file '{csv_file}' does not exist.")

    # Check if acclimation_time is greater than 0
    #if parsed.acclimation_time < 0:
    #    parser.error("Acclimation time must be greater than 0.")

    # Check if trial_bin is greater than 0
    if parsed.trial_bin < 0:
        parser.error("Frequency sample size must be greater than 0.")

    # Check if bin_time is greater than 0
    if parsed.bin_time < 0:
        parser.error("Bin time must be greater than 0.")

    # Check if last_x_percent is in the range (0, 100)
    if not 0 < parsed.last_x_percent < 100:
        parser.error("last_x_percent must be in the range (0, 100).")

    return parsed


# open and format animal data
def time_from_file(filename):
    "Convert filename to timestamp for start of trials"
    datetime = re.findall("\d\d_\d\d_\d\d_T_\d\d_\d\d_\d\d", filename)[0]
    datetime = datetime.replace("_", "-", 2).replace("_T_", " ").replace("_", ":")
    return pd.Timestamp(datetime)


def load_file(filename):
    '''Given a string filename, loads in the data, extracts the start time from the filename, and formats the timestamps based on the start time.'''
    animal = pd.read_csv(filename, header=None)
    animal = animal.rename(columns={0:"timestamp",1:"poke", 2:"lick", 3:"water", 4:"delay"})
    datetime = time_from_file(filename)
    animal["timestamp"] = pd.to_datetime(animal["timestamp"], unit="s", origin=pd.Timestamp(datetime))
    animal["delay"] = pd.to_timedelta(animal["delay"], unit="ms")
    return animal.iloc[:, 0:5]

def get_trials (trial):
    # generates a list of idx of length x
    (trial_no, len) = trial
    return pd.Series([trial_no for i in range(len)], dtype="int")

def enumerate_trials(animal):
    '''Given the set of data for an animal, labels each sample in a trial with the trial type and the trial number, numbering trials consecutively from 1.
    The number of trials, and the number of each trial type, matches v16 of the matlab analysis.
    '''
    #first sample, all samples where current water code is not 7 (timeout) and previous is 7,  last sample
    trial_boundary = pd.concat([animal.iloc[[0]], animal[(animal["water"]!=7) & (animal["water"].shift() == 7)], animal.iloc[[animal.shape[0]-1]]])
    #indexes of first sample in each trial
    trial_boundary_indicies = pd.Series(trial_boundary.index)

    #get number of samples per trial (difference in sample number between previous index and current index)
    num_samples = trial_boundary_indicies.diff().fillna(0).astype('int').tolist()
    #label with index
    trial_count = pd.Series(enumerate(num_samples))

    #label all samples in a trial with its trial number
    trial_count = trial_count.map(get_trials).explode()[1:]
    trial_count.index = range(0, trial_count.shape[0])
    #add trial labels to animal data
    animal.insert(1, "trial no", trial_count)

    # set trial no for last sample
    animal.loc[animal.shape[0] - 1, "trial no",] = trial_boundary.shape[0] - 1

    return animal

def label_trials(animal):
    # get all trials labeled with a 3 (water trials)
    go = animal.groupby(["trial no"]).filter(lambda x: (x["water"]==3).any())
    # label all these trials as water
    go["trial type"] = ["water" for i in range(go.shape[0])]
    animal["trial type"] = go["trial type"]
    # label remaning trials as blank
    animal["trial type"] = animal["trial type"].fillna(value="blank")
    return animal

def format_data(data):
    '''Set animal and condition, update lick value to binary 1/0, format acclimation time, and drop last sample of each trial'''
    data.loc[data["lick"] == 2,"lick"] = 1
    data["acc"] = pd.to_timedelta(data["acc"], unit="day")
    # drop last sample since it is duplicate
    data = data[data.groupby(["condition", "animal", "trial no"]).cumcount(ascending=False) > 0]
    return data

def make_animal_df(andir, metadata, animal_name):
    fs = os.listdir(andir)
    #ensure files concatenated in time order
    fs.sort()
    animal = []        
    #this will load all files in a given directory - should only have the relevant training files
    for f in fs:
        f_path = andir + "\\" + f
        if (os.path.isfile(f_path)):      
            animal.append(load_file(f_path))

    animal = pd.concat(animal, ignore_index=True)
    animal = enumerate_trials(animal)
    animal = label_trials(animal)

    #load metadata if the animal has it
    if not (metadata[metadata["Animal ID"] == animal_name].empty):     
        an_meta =  metadata[metadata["Animal ID"] == animal_name].reset_index()
        animal["acc"] = an_meta.loc[0, "ACC days"]
        animal["age"] = an_meta.loc[0, "Age"]
        animal["sex"] = an_meta.loc[0, "Sex"]
        animal["strain"] = an_meta.loc[0, "Strain"]
        animal["animal"] = animal_name
    
    animal = format_data(animal) 
    return animal

def make_condition_df(condir, condition, metadata):
    andirs = os.listdir(condir)
    animals = []
    for animal_name in andirs:
        animal_path = condir + "\\" + animal_name
        animals.append(make_animal_df(animal_path, metadata, animal_name))
    animals = pd.concat(animals, ignore_index=True)
    animals["condition"] = condition
    return animals

# analysis functions

def rolling_frequency_average(data, freq_window, values, keep, index):
    '''Calculate licking frequency as a rolling average of 
    Parameters:
    values - set of columns to calculate rolling average
    keep - set of columns to keep but not calculate rolling average
    index - set of columns to index/group by
    '''
    group = data.set_index("timestamp").groupby(index)
    licks = group.rolling(window=pd.to_timedelta(freq_window, unit="ms"))[values].sum()/(freq_window/1000)
    keep = data.set_index(index)[keep].set_index("timestamp",append=True)
    data_roll = pd.concat([licks, keep], axis=1).reset_index()

    return data_roll

def resample_align(data, frq, key, keep, index):
    '''Resample at the given freq. 100ms aligns data without changing values'''
    gp = pd.Grouper(key=key, freq=pd.offsets.Milli(frq))
    #append resampler to index group
    index.append(gp)
    group = data.groupby(index)
    keep = group[keep].first()
    data_sample = keep.reset_index()

    return data_sample

def get_trial_start(data, index, key):
    data = data.set_index(index)
    data["trial time"] = data.groupby(index)[key].first()
    data = data.reset_index()
    return data

def bin_by_time(data, trial_bin, index, key):
    gp = pd.Grouper(key=key, freq=trial_bin)
    index.append(gp)
    group = data.groupby(index)
    data = group.first()
    data = data.reset_index()

def delta(data, index, key):
    '''calculates time since start of SAT. key is column to use to calculate'''
    data = data.set_index(index)
    delta = data[key] - data.sort_values(key).groupby("animal")[key].first() - data.groupby("animal")["acc"].first()
    data["delta"] = delta
    return data.reset_index()

def deliverydelta(data, index):
    '''find start of air delivery (trial start + delay to puff)'''
    data = data.set_index(index)
    grouped = data.groupby(index)
    data["air start"] = grouped["timestamp"].first() + grouped["delay"].first()
    data["delivery delta"]  = data["timestamp"] - data["air start"]
    data = data.drop(columns="air start")
    return data.reset_index()

def drop_bins(data, min_trials, min_blank, min_water, index, key):
    gp_index = index.copy()
    gp_index.append(key)
    group = data.groupby(gp_index)
    #filter bins with fewer than min_blank blank or min_water water trials
    counts = group.count().unstack(level=key)
    cond = (counts.loc[:, ("trial no", "blank")] > min_blank) & (counts.loc[:, ("trial no", "water")] > min_water) 
    data = data[cond]
    #filter bins with fewer than min_trials total trials
    total_group = data.groupby(index) 
    total_counts = total_group.count()
    data = data[total_counts > min_trials]
    data = data.reset_index()
    return data


# analysis: traces over entire trial -> rolling window with bin size flexibiliy 
#           bin over training time -> bin size flexibility
#           drop trials -> number of trials flexibility
# plots: all bins average trace plot
def calculate_statistics(data, freq_window, freq_bin, time_bin, min_trials, min_water, min_blank):
    

# plot all traces and average trace for a bin
if name == "__main__":
    
