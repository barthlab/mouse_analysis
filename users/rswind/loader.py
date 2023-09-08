# Functions to load and format raw data from cages
# Last updated: 7/7/2023
# Author: Rachel Swindell

# TODO: slow - base off of mouse_analysis.py for speedups?

import re
import os
import sys

import pandas as pd

# open and format animal data
def time_from_file(filename):
    "Convert filename to timestamp for start of trials and return timestamp."

    datetime = re.findall("\d\d_\d\d_\d\d_?~?T_?~?\d\d_\d\d_\d\d", filename)[0]
    datetime = datetime.replace("_", "-", 2).replace("~", "-", 2).replace("_T_", " ").replace("_", ":")
    return pd.Timestamp(datetime)

# TODO: change  file loading to handle pseudo data (include 6th column to id stimulus vs blank trials)
def load_file(filename):
    '''Loads, formats, and returns lick frequecny data as a data frame.
    
    Loads data from csv format. Extracts the start time from the filename and uses it to format timestamps from milliseconds since start of file to datetimes.
    Converts delays at beginning of trial to timedeltas.

    Parameters:
    filename --- path of file to load (string) '''

    animal = pd.read_csv(filename, header=None)
    animal = animal.rename(columns={0:"timestamp",1:"poke", 2:"lick", 3:"water", 4:"delay"})
    datetime = time_from_file(filename)
    animal["timestamp"] = pd.to_datetime(animal["timestamp"], unit="s", origin=pd.Timestamp(datetime))
    animal["delay"] = pd.to_timedelta(animal["delay"], unit="ms")
    return animal.iloc[:, 0:5]

def get_trials(trial):
    '''Generates a pandas Series of a given length populated by a given value.
    
    Used to label all samples in a trial with the trial number.

    Parameters
    trial --- tuple of (value, length)'''

    # generates a list of idx of length x
    (trial_no, len) = trial
    return pd.Series([trial_no for i in range(len)], dtype="int")

def enumerate_trials(animal):
    '''Numbers each sample in a trial with the trial number and returns the input dataframe with labeled trials.

    Given the set of data for an animal, labels each sample in a trial with the trial type and the trial number, numbering trials consecutively from 1.
    The number of trials, and the number of each trial type, matches v16 of the matlab analysis.
    '''
    # first sample, all samples where current water code is not 7 (timeout) and previous is 7,  last sample
    trial_boundary = pd.concat([animal.iloc[[0]], animal[(animal["water"]!=7) & (animal["water"].shift() == 7)], animal.iloc[[animal.shape[0]-1]]])
    # indexes of first sample in each trial
    trial_boundary_indicies = pd.Series(trial_boundary.index)

    # get number of samples per trial (difference in sample number between previous index and current index)
    num_samples = trial_boundary_indicies.diff().fillna(0).astype('int').tolist()
    # label with index
    trial_count = pd.Series(enumerate(num_samples))

    # label all samples in a trial with its trial number
    trial_count = trial_count.map(get_trials).explode()[1:]
    trial_count.index = range(0, trial_count.shape[0])
    # add trial labels to animal data
    animal.insert(1, "trial no", trial_count)

    # set trial no for last sample
    animal.loc[animal.shape[0] - 1, "trial no",] = trial_boundary.shape[0] - 1

    return animal

# TODO: change labeling to handle pseudo data (label as stimulus vs blank and as water vs no water -  check how mouse_analysis is doing it)
# might need to do it differently since it uses lists vs series/dataframes
def label_trials(animal):
    '''Decides what type a trial is (water or no water) and return animal with all trials labeled.
    
    Trial types are decided by the code in the 3rd column (if the trial contains a 3, it is a water trial). 
    Requires trials to already be grouped individually and labeled with trial number.
    '''
    # get all trials labeled with a 3 (water trials)
    go = animal.groupby(["trial no"]).filter(lambda x: (x["water"]==3).any())
    # label all these trials as water
    go["trial type"] = ["water" for i in range(go.shape[0])]
    animal["trial type"] = go["trial type"]
    # label remaning trials as blank
    animal["trial type"] = animal["trial type"].fillna(value="blank")
    return animal

# TODO: handle acclimation as an hourly input
def format_data(data):
    '''Make small formatting changes and return input dataframe.

    Formatting chagnes made: 
     - update lick column to binary values
     - format acclimation time to day-based timedelta
     - drop last sample of each trial, as it is a duplicate sample that indicates the end of a trial and holds no information'''
    
    data.loc[data["lick"] == 2,"lick"] = 1
    data.loc[:, "acc"] = pd.to_timedelta(data["acc"], unit="day")
    data = data[data.groupby(["animal", "trial no"]).cumcount(ascending=False) > 0]

    return data

# TODO: handle missing metadata file more gracefully 
# TODO: handle acclimation in hours or days
def make_animal_df(andir, metadata, animal_name, acc_col_name, default_acc):
    '''Load and format all data files for one animal and return as a data frame.

    Will fail if there are non-data files in the given directory or if animal id column in metadata file is not named correctly.
    Prints failure warning messages if metadata is missing. Specifically looks for number of days (or hours) of acclimation, and
    sex, age, and strain of animal.

    Parameters:
    andir --- directory containing the data files for that animal (and no other files)
    metadata --- dataframe containing metadata about an animal associated with animal id
                 if animal id column is not named exactly 'Animal ID', will fail
    animal_name --- ID code of animal
    acc_col_name --- name of column containing length of acclimation in metadata file
    default_acc --- default length of acclimation to assume if none is listed'''

    fs = os.listdir(andir)
    # ensure files concatenated in time order
    fs.sort()
    animal = []        
    # this will load all files in a given directory - should only have the relevant training files
    for f in fs:
        f_path = andir + "\\" + f
        if (os.path.isfile(f_path)): 
            try:     
                animal.append(load_file(f_path))
            except UnicodeDecodeError:
                print(f"Check that only behavior files are in {andir}")
                sys.exit(1)
    if animal == []:
        print(f"{animal_name} has no behavior files")
        return pd.DataFrame()
    else:
        animal = pd.concat(animal, ignore_index=True)
        animal = enumerate_trials(animal)
        animal = label_trials(animal)

        # load metadata if the animal has it
        try:
            metadata[metadata["Animal ID"] == animal_name].empty
        except KeyError:
            print("Animal ID column must be named 'Animal ID'")
            sys.exit(1)

        if not (metadata[metadata["Animal ID"] == animal_name].empty):
            an_meta =  metadata[metadata["Animal ID"] == animal_name].reset_index()
            try:
                animal["acc"] = an_meta.loc[0, acc_col_name]
            except KeyError:
                print(f"Number of acclimation days not named {acc_col_name}")
            try:
                animal["age"] = an_meta.loc[0, "Age"]
            except:
                animal["age"] = pd.NA
                print(f"{animal_name} Age not in metadata file")
            try:
                animal["sex"] = an_meta.loc[0, "Sex"]
            except KeyError:
                animal["sex"] = ''
                print(f"{animal_name} Sex not in metadata file")
            try:
                animal["strain"] = an_meta.loc[0, "Strain"]   
            except KeyError:
                animal["strain"] = ''
                print(f"{animal_name} Strain not in metadata file")
        else:
            print(f"{animal_name}: no metadata")
            animal["acc"] = default_acc
    
            
        animal["animal"] = animal_name
        animal = format_data(animal) 
        return animal

def make_condition_df(condir, condition, metadata, acc_col_name, default_acc):
    '''Load and format all data files for all animals in a condition, and return as a data frame.
    
    Requires all animals to be included to be in the same directory, which should have no other folders or files in it.
    
    Parameters:
    condit --- path of directory containing animal subdirectories which contain data files
    condition --- name of condition (e.g. Acc, SAT1)
    metadata --- data frame containing relevant metadata about animals in condition (can contain other information as well)
    acc_col_name --- name of column containing length of acclimation in metadata file
    default_acc --- default length of acclimation to assume if none is listed'''
    
    andirs = os.listdir(condir)
    animals = []
    for animal_name in andirs:
        animal_path = condir + "\\" + animal_name
        animals.append(make_animal_df(animal_path, metadata, animal_name, acc_col_name, default_acc))
    animals = pd.concat(animals, ignore_index=True)
    animals["condition"] = condition
    return animals
