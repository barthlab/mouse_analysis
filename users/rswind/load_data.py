import re
import os
import numpy as np
import pandas as pd

def timefromfile(filename):
    "Convert filename to timestamp for start of trials"
    datetime = re.findall("\d\d_\d\d_\d\d_T_\d\d_\d\d_\d\d", filename)[0]
    datetime = datetime.replace("_", "-", 2).replace("_T_", " ").replace("_", ":")
    return pd.Timestamp(datetime)


def loadfile(filename):
    '''Given a string filename, loads in the data, extracts the start time from the filename, and formats the timestamps based on the start time.'''
    animal = pd.read_csv(filename, header=None)
    animal = animal.rename(columns={0:"timestamp",1:"poke", 2:"lick", 3:"water", 4:"delay"})
    datetime = timefromfile(filename)
    animal["timestamp"] = pd.to_datetime(animal["timestamp"], unit="s", origin=pd.Timestamp(datetime))
    animal["delay"] = pd.to_timedelta(animal["delay"], unit="ms")
    return animal.iloc[:, 0:5]

def getTrials (tp):
    (idx, x) = tp
    return pd.Series([idx for i in range(x)], dtype="int")

def labelTrials (animal):
    '''Given the set of data for an animal, labels each sample in a trial with the trial type and the trial number, numbering trials consecutively from 1.
    The number of trials, and the number of each trial type, matches v16 of the matlab analysis.
    '''
    #first sample for each new trial - current water code is not 7 (timeout) and previous is 7 plus first and last sample (for diff alignment)
    trial_boundaries = pd.concat([animal.iloc[[0]], animal[(animal["water"]!=7) & (animal["water"].shift() == 7)], animal.iloc[[animal.shape[0]-1]]])

    trial_count = pd.Series(trial_boundaries.index)    
    #get number of samples per trial (difference in sample number between previous boundary and current boundary)
    trial_count = pd.Series(enumerate(trial_count.diff().fillna(0).astype('int').tolist()))
    #enumerate trial number from first sample to all samples in trial
    trial_count = trial_count.map(getTrials).explode()[1:]
    trial_count.index = range(0, trial_count.shape[0])
    animal.insert(1, "trial no", trial_count)
    #set trial no for last sample
    animal.loc[animal.shape[0] - 1, "trial no",] = trial_boundaries.shape[0] -1

    #label trial types
    go = animal.groupby(["trial no"]).filter(lambda x: (x["water"]==3).any())
    go["trial type"] = ["water" for i in range(go.shape[0])]
    animal["trial type"] = go["trial type"]
    animal["trial type"] = animal["trial type"].fillna(value="blank")
    return animal

def init_data(data, animal, cond):
    '''Set animal and condition and update lick value to binary 1/0'''
    data.loc[data["lick"] == 2,"lick"] = 1
    data["animal"] = animal
    data["condition"] = cond
    return data

def load_sat(pSAT, animals, name, hm4di_metadata, metadata):
    for d in os.listdir(pSAT):
        print(d)

        p = pSAT + "\\" + d
        #guarantees files opened in time order (path list is arbitrary)
        fs = os.listdir(p)
        fs.sort()
        curranimal = [] 
        
        #all files concatenated in order
        time = 0
        SAT_start = 0
        for f in fs:
            if (os.path.isfile(p + "\\" + f)):      
                curranimal.append(loadfile(p + "\\" + f))        
                if time == 0:
                    time = timefromfile(p + "\\" + f)             
        
        curranimal = labelTrials(pd.concat(curranimal, ignore_index=True))

        #align animal metadata
        if not (metadata[metadata["Animal ID"] == d].empty): 
            airon = metadata[metadata["Animal ID"] == d].reset_index().loc[0, "air on"]
            if not pd.isna(airon):
                t = re.findall("\d?\d",airon)

            an_meta =  metadata[metadata["Animal ID"] == d].reset_index()
            acc = an_meta.loc[0, "ACC days"]
            curranimal["age"] = an_meta.loc[0, "Age"]
            curranimal["sex"] = an_meta.loc[0, "Sex"]
            curranimal["strain"] = an_meta.loc[0, "Strain"]
            SAT_start = pd.Timestamp(year=time.year, month=time.month, day=time.day, hour=int(t[0]), minute=int(t[1])) + pd.to_timedelta(acc, unit="D")  
        elif not (hm4di_metadata[hm4di_metadata["Animal ID"] == d].empty):
            acc = hm4di_metadata[hm4di_metadata["Animal ID"] == d].reset_index()["Acclimation (days)"][0]
            SAT_start  = pd.Timestamp(year=time.year, month=time.month, day=time.day, hour=12) + pd.to_timedelta(acc, unit="D")  
        else:
            print("no metadata")
            
        curranimal = init_data(curranimal, d, name)
        curranimal["SAT_start"] = SAT_start
        curranimal["acc"] = acc
        animals = pd.concat([animals, curranimal])
    return animals

def load_100(animals, path_100_80, metadata):
    for d in os.listdir(path_100_80):
        print(d)
        
        p = path_100_80 + "\\" + d
        #animal directory (path list is arbitrary)
        andirs = os.listdir(p)
        andirs = sorted(andirs, key=str.lower)
        curranimal = []
        for d2 in andirs:
            p2 = p + "\\" + d2
            #acc/SAT - sort if multiple files
            fs = os.listdir(p2)
            fs.sort()
            SAT_start = 0
            for f in fs:            
                if os.path.isfile(p2 + "\\" + f):
                    curranimal.append(loadfile(p2 + "\\" + f))
                    if (d2 == "SAT") & (SAT_start == 0):
                        SAT_start = timefromfile(p2 + "\\" + f)
        
        #all files concatenated in order   
        curranimal = labelTrials(pd.concat(curranimal, ignore_index=True))
        curranimal["SAT_start"] = SAT_start
        
        curranimal = init_data(curranimal, d, "SAT_100_80")
        if not (metadata[metadata["Animal ID"] == d].empty): 
            airon = metadata[metadata["Animal ID"] == d].reset_index().loc[0, "air on"]
            if not pd.isna(airon):
                t = re.findall("\d?\d",airon)
                
            an_meta =  metadata[metadata["Animal ID"] == d].reset_index()
            curranimal["age"] = an_meta.loc[0, "Age"]
            curranimal["sex"] = an_meta.loc[0, "Sex"]
            curranimal["strain"] = an_meta.loc[0, "Strain"]
            curranimal["acc"] = an_meta.loc[0, "ACC days"]
        
        animals = pd.concat([animals, curranimal])
        animals["acc"] = animals["acc"].fillna(2)
    return animals

def load_all_data(SATpaths, path_100_80, hm4di_metadata, metadata):
    animals = pd.DataFrame()
    
    for (pSAT, name) in SATpaths:
        print(name)
        animals = pd.concat([animals, load_sat(pSAT, pd.DataFrame(), name, hm4di_metadata, metadata)])
    animals = pd.concat([animals, load_100(pd.DataFrame(), path_100_80, metadata)])
    
    animals = animals.set_index(["condition", "animal", "trial no"])

    return animals