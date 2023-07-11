# Functions to load and format raw data from cages
# Last updated: 7/7/2023
# Author: Rachel Swindell

# TODO: slow - base off of mouse_analysis.py for speedups?

import re
import os

import pandas as pd

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
    data.loc[:, "acc"] = pd.to_timedelta(data["acc"], unit="day")
    # drop last sample since it is duplicate
    data = data[data.groupby(["animal", "trial no"]).cumcount(ascending=False) > 0]
    return data

def make_animal_df(andir, metadata, animal_name, acc_col_name, default_acc):
    
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
        animal["acc"] = an_meta.loc[0, acc_col_name]
        animal["age"] = an_meta.loc[0, "Age"]
        animal["sex"] = an_meta.loc[0, "Sex"]
        animal["strain"] = an_meta.loc[0, "Strain"]   
    else:
        print(f"{animal_name}: no metadata")
        animal["acc"] = default_acc

    animal["animal"] = animal_name
    animal = format_data(animal) 
    return animal

def make_condition_df(condir, condition, metadata, acc_col_name, default_acc):
    andirs = os.listdir(condir)
    animals = []
    for animal_name in andirs:
        animal_path = condir + "\\" + animal_name
        animals.append(make_animal_df(animal_path, metadata, animal_name, acc_col_name, default_acc))
    animals = pd.concat(animals, ignore_index=True)
    animals["condition"] = condition
    return animals
