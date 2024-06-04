# TODO: add scrollbar when window gets too large for screen

# TODO: updates during analysis and plotting
# TODO: add separate frame for error messages to prevent stupid resizing


# -*- coding: utf-8 -*-
"""
Created 2024-05-29
Last Edited 2024-06-03

@author: Rachel Swindell
"""

import os
from pathlib import Path
import re
import datetime
import math
from collections import Counter

import tkinter as tk
from tkinter import ttk, filedialog, StringVar, BooleanVar

import pandas as pd
import numpy as np

import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection

import seaborn as sns #v0.13.0

# Create a Tkinter window
root = tk.Tk() #Makes the window
root.title("Behavior Analysis") #Makes the title that will appear in the top left
root.bind_all('<Escape>', lambda e: e.widget.focus_set())

class WrappingLabel(tk.Label):
    def __init__(self, master=None, **kwargs):
        tk.Label.__init__(self, master, **kwargs)
        self.bind('<Configure>', lambda e: self.config(wraplength=self.winfo_width()))


class ResizeEqualLabelFrame(ttk.LabelFrame):
    def __init__(self, master=None,cols=0, rows=0, **kwargs):
        ttk.LabelFrame.__init__(self, master, **kwargs)
        for i in range(cols):
            self.grid_columnconfigure(i, weight=1)
        for i in range(rows):
            self.grid_rowconfigure(i, weight=1)

class EnterButton(ttk.Button):
    def __init__(self, master=None,**kwargs):
        ttk.Button.__init__(self, master, **kwargs)
        self.bind('<Return>', lambda e: self.invoke())

class EnterCheckbutton(ttk.Checkbutton):
    def __init__(self, master=None,**kwargs):
        ttk.Checkbutton.__init__(self, master, **kwargs)
        self.bind('<Return>', lambda e: self.invoke())


class mplEmbedPlot():
    def __init__(self,fig=None, master=None,**kwargs):
        self.frame = ttk.Frame(master, **kwargs)
        self.canvas = FigureCanvasTkAgg(fig, master=self.frame)
        self.canvas.draw()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.frame, pack_toolbar=False)
        self.toolbar.update()
        self.fontlabel = tk.Label(text='Font size')
        #self.fontentry = ttk.Entry()

        #self.fontentry.bind('<FocusOut>', self.update_font)
        #self.fontentry.bind('<Return>', self.update_font)
        self.frame.bind_all('<Escape>', lambda e: e.widget.focus_set())

        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        #self.fontlabel.pack(side=tk.BOTTOM)
        #self.fontentry.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        

    def update_font(self, e):
        size = int(self.fontentry.get())
        for ax in self.canvas.figure.axes:
            ax.tick_params(axis='both', which='major', labelsize=size)
            ax.xaxis.label.set_size(size)
            ax.yaxis.label.set_size(size)
            ax.title.label.set_size(size)
        self.canvas.figure.suptitle.set_size(size)
        

#####################################
# Functions to Load and Format Data #
#####################################

def time_from_file(filename):
    "Convert filename to timestamp."
    dt = re.findall("\d\d_\d\d_\d\d_?~?T_?~?\d\d_\d\d_\d\d", filename)[0]
    dt = dt.replace("_", "-", 2).replace("~", "-", 2).replace("_T_", " ").replace("-T-", " ").replace("_", ":")
    return datetime.datetime.strptime(dt,'%m-%d-%y %H:%M:%S')

def set_noon(t):
    "Get noon on start of training as timestamp"
    return datetime.datetime(t.year, t.month, t.day, 12, 0, 0)

def load_file(filename):
    '''Loads, formats, and returns lick frequecny data as a data frame.
    
    Loads data from csv format. Extracts the start time from the filename and uses it to format timestamps from milliseconds since start of file to datetimes.
    Converts delays at beginning of trial to timedeltas.

    Parameters:
    filename --- path of file to load (string) '''

    data = pd.read_csv(filename, header=None)

    # rename columns
    mp = {0:"timestamp",1:"poke", 2:"lick", 3:"condition code", 4:"delay"}
    if len(data.columns) == 6:
        mp = {0:"timestamp",1:"poke", 2:"lick", 3:"condition code", 4:"delay", 5:"stimulus"}
    data = data.rename(columns=mp)
    
    # convert time columns to correct type
    dt = time_from_file(filename)
    # convert to int in ms to avoid floating point rounding errors in datetime conversion
    data['timestamp'] = (data['timestamp']*1000).astype(int)
    data["timestamp"] = pd.to_datetime(data["timestamp"], unit="ms", origin=dt)
    data["delay"] = pd.to_timedelta(data["delay"], unit="ms")

    # get noon for timebin alignment
    data['offset'] = set_noon(dt)

    return data

def get_trials(trial):
    '''Generates a pandas Series of x of length len.
    
    Used to label all samples in a trial with the trial number.

    Parameters
    trial --- (value, length)'''
    (x, len) = trial
    # generates a list of x of length len
    return pd.Series([x for i in range(len)], dtype="int")

def enumerate_trials(data, start=0):
    '''Labels each sample in a trial with the trial type and the trial number, numbering trials consecutively from `start`.
    The number of trials, and the number of each trial type, matches v16 of the matlab analysis. Returns the input dataframe with labeled trials.
    '''
    # get first sample in each trial
    # first sample, all samples where current water code is not 7 (timeout) and previous is 7,  last sample
    trial_boundary = pd.concat([data.iloc[[0]], data[(data["condition code"]!=7) & (data["condition code"].shift() == 7)], data.iloc[[data.shape[0]-1]]])
    # indexes of first sample in each trial
    trial_boundary_indicies = pd.Series(trial_boundary.index)

    # get number of samples per trial 
    # (difference in index number between previous sample and current sample)
    num_samples = trial_boundary_indicies.diff().fillna(0).astype('int').tolist()
    # label with number of samples with index
    trial_count = pd.Series(enumerate(num_samples))

    # label all samples in a trial with its trial number
    trial_count = trial_count.map(get_trials).explode()[1:]
    trial_count.index = range(0, trial_count.shape[0])
    trial_count = trial_count + start
    
    # add trial labels to data
    data["trial no"] = trial_count
    # set trial no for last sample
    data.loc[data.shape[0] - 1, "trial no",] = trial_boundary.shape[0] - 1 + start

    return data

def label_trials_stimulus(data):
    '''Label samples as stimulus or blank.'''
    data["stimulus"] = data["stimulus"].replace({0:"blank", 1:"stimulus"})
    return data

def label_trials_water(data):
    '''Return data with all trials labeled as water or no water.
    
    Trial types are decided by the code in the 3rd column (if the trial contains a 3, it is a water trial). 
    Requires trials to already be grouped individually and labeled with trial number.
    '''
    # get all trials labeled with a 3 (water trials)
    go = data.groupby(["trial no"]).filter(lambda x: (x["condition code"]==3).any())
    # label all these trials as water
    go["water"] = ["water" for i in range(go.shape[0])]
    data["water"] = go["water"]
    # label remaning trials as blank
    data["water"] = data["water"].fillna(value="no water")
    return data

def label_trials(data):
    '''Label trials with both stimulus or blank and water or no water distinction.'''
    data = label_trials_water(data)
    if "stimulus" in data.columns:
        data = label_trials_stimulus(data)
    else:
        data["stimulus"] = data["water"]
        data["stimulus"] = data["stimulus"].replace({"water":"stimulus", "no water":"blank"})
    data["type"] = data["water"] + " & " + data["stimulus"]
    return data

def drop_last_sample(data):
    '''Drops last sample of each trial. This has the same timestamp as 
    the second-to-last sample and only indicates the trial is over.'''
    data = data[data.groupby(["animal", "trial no"]).cumcount(ascending=False) > 0]
    return data

def make_animal_df(andir, metadata, animal_name, meta_vals,mssg, root):
    '''Load and format all data files for one animal and return as a data frame.

    Requires:
     --- `andir` only contains data files (may contain directories but not other files)
     --- metadata file is provided
     --- metadata contains animal ID in column named "Animal ID"
     --- metadata contains acclimation time in column named "acc"

    If name of animal in metadata file does not match `animal_name` no data will
    be loaded. Column names in addition to "acc" can be provided to load metadata 
    from those columns in the order they are provided. Prints a warning message
    for all additional columns that were provided but had no data, and returns
    those columns as NaN. For animals that have training split into multiple files,
    the entire filename except for the timestamp _must_ be identical, or they
    will not be loaded in the correct order.

    Parameters:
    andir --- directory containing the data files for that animal (and no other files)
    metadata --- dataframe containing metadata about an animal associated with animal id
                 if animal id column is not named exactly 'Animal ID', will fail
    animal_name --- ID code of animal
    meta_vals --- list of columns to include metadata from. Must include 'acc' for length of acclimation (in days).'''
    
    if metadata.empty:
        print("Must include metadata file including at least acclimation time")
        mssg.set("Must include metadata file including at least acclimation time")
        return
    if 'acc' not in metadata.columns:
        print("Column indicating acclimation time must be named 'acc'")
        mssg.set("Column indicating acclimation time must be named 'acc'")
        return
    if 'acc' not in meta_vals:
        print("Must include acc in list of columns to be selected from metadata")
        mssg.set("Must include acc in list of columns to be selected from metadata")
        return

    fs = os.listdir(andir)
    # ensure files concatenated in time order
    # only works if rest of filename is identical
    fs.sort()
    animal = []        
    # load all files in an animal's folder
    start = 0
    for f in fs:
        f_path = andir + "\\" + f
        if (os.path.isfile(f_path)): 
            try:     
                data = load_file(f_path)
            except UnicodeDecodeError:
                print(f"Check that only behavior files are in {andir}")
                mssg.set(f"Check that only behavior files are in {andir}")
                return
            # labeling trials
            data = enumerate_trials(data, start)
            start = data.tail(1)["trial no"].reset_index(drop=True)[0]
            data = label_trials(data)
            animal.append(data)
    # load metadata
    if animal == []:
        print(f"{animal_name} has no behavior files")
        mssg.set(mssg.get() + f"\n{animal_name} has no behavior files. Check if this folder should be run as an animal or as a condition.")
        return pd.DataFrame()
    else:
        animal = pd.concat(animal, ignore_index=True)
        try:
            metadata[metadata["Animal ID"] == animal_name].empty
        except KeyError:
            print("Animal ID column must be named 'Animal ID'")
            mssg.set("Animal ID column must be named 'Animal ID'")
            return
        if not (metadata[metadata["Animal ID"] == animal_name].empty):
            an_meta =  metadata[metadata["Animal ID"] == animal_name].reset_index()
            fail_keys = []
            for key in meta_vals:
                try:
                    animal[key] = an_meta.loc[0, key]
                    if an_meta.loc[:,key].isna().any():
                        if key == 'acc':
                            print(f'Acclimation time (required) is missing for {animal_name}. Skipping animal.')
                            mssg.set(mssg.get() + f'\nAcclimation time (required) is missing for {animal_name}. Skipping animal.')
                            return pd.DataFrame()
                        else:
                            fail_keys += [key]
                except KeyError:
                    fail_keys += [key]
                    if key == 'acc':
                        print("Metadata file must include acclimation time in column named 'acc'.")
                        mssg.set("Metadata file must include acclimation time in column named 'acc'.")
                        return
            if len(fail_keys) > 0:
                print("%s %s not in metadata file" % (animal_name, ', '.join(fail_keys)))     
                mssg.set(mssg.get() + ("\n%s %s not in metadata file" % (animal_name, ', '.join(fail_keys)))) 
        else:
            print(f"{animal_name}: no metadata - must have at least acclimation time. Skipping animal.")
            mssg.set(mssg.get() + f"\n{animal_name}: no metadata - must have at least acclimation time. Skipping animal.")
            return pd.DataFrame()
        
        animal["animal"] = animal_name
        animal['offset'] = animal.loc[0, 'offset']
        animal[fail_keys] = pd.NA
        # one-hot encode lick frequency
        animal.loc[animal["lick"] == 2,"lick"] = 1
        # convert acclimation time to timedelta
        animal.loc[:, "acc"] = pd.to_timedelta(animal["acc"], unit="day")
        animal = drop_last_sample(animal)
        return animal

##########################
# Functions for Plotting #
##########################
def style_axes_helper(ax, title, xlabel, ylabel, xlim, ylim, xmajlocator, xminlocator,
               ymajlocator, yminlocator):
    ax.set_title(f"{title}", fontsize=12)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if xlabel != None:
        ax.set_xlabel(xlabel, fontsize=10)
    if ylabel != None:
        ax.set_ylabel(ylabel, fontsize=10)    
    if xmajlocator != None: ax.xaxis.set_major_locator(xmajlocator)
    if xminlocator != None: ax.xaxis.set_minor_locator(xminlocator)
    ax.tick_params(axis='both', which='major', labelsize=10)
    if ymajlocator != None: ax.yaxis.set_major_locator(ymajlocator)
    if yminlocator != None: ax.yaxis.set_minor_locator(yminlocator)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, xminmult, ymajmult,
               yminmult, hspace=0.45, wspace=0.25, suptitle=None,title=None):
    xmaj = ticker.MultipleLocator(xmajmult) if xmajmult != None else xmajmult
    xmin = ticker.MultipleLocator(xminmult) if xminmult != None else xminmult
    ymaj = ticker.MultipleLocator(ymajmult) if ymajmult != None else ymajmult
    ymin = ticker.MultipleLocator(yminmult) if yminmult != None else yminmult

    t = title
    if g.axes_dict: 
        for (name,ax) in g.axes_dict.items():
            if title == None:
                t=name
            style_axes_helper(ax, t, xlabel, ylabel, xlim, ylim, xmaj, xmin, ymaj, ymin)
    else:
        style_axes_helper(g.figure.axes[0], '', xlabel, ylabel, xlim, ylim, xmaj, xmin, ymaj, ymin)

    g.figure.suptitle(suptitle, fontsize=12, y=0.98, x=0.5)
    g.figure.subplots_adjust(hspace=hspace, wspace=wspace)

def style_hr(g, xlabel, ylabel, xlim, ylim, hspace=0.45, wspace=0.25, 
             suptitle=None, title=None, ymajmult=2, yminmult=1, xmajmult=12):
    style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, 4, ymajmult, yminmult, 
               hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)  

def style_ms(g, xlabel, ylabel, xlim, ylim, xmajmult=500, xminmult=100, 
             ymajmult=2, yminmult=1,hspace=0.45, wspace=0.25, 
             suptitle=None,style=True):
    style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, xminmult, ymajmult, yminmult, hspace=hspace,
               wspace=wspace, suptitle=suptitle)
    if style:
        if g.axes_dict:
            for (name,ax) in g.axes_dict.items():
                ax.add_patch(Rectangle((0, ylim[0]), 500, ylim[1] - ylim[0], color="lightgrey", alpha=0.4, 
                                    zorder=0, fill=True))        
                ax.add_patch(Rectangle((1005, ylim[0]), 100, ylim[1] - ylim[0], color="#BBE9FD", alpha=0.8, 
                                        zorder=0, fill=True))
        else:
            ax = g.figure.axes[0]
            ax.add_patch(Rectangle((0, ylim[0]), 500, ylim[1] - ylim[0], color="lightgrey", alpha=0.4, 
                            zorder=0, fill=True))        
            ax.add_patch(Rectangle((1005, ylim[0]), 100, ylim[1] - ylim[0], color="#BBE9FD", alpha=0.8, 
                                    zorder=0, fill=True))
  
def style_performance(g):
    if g.axes_dict:
        for (name,ax) in g.axes_dict.items(): 
            ax.axhline(y=0, xmin=0, xmax=1, ls="-", color="black", zorder=10, lw=1)
    else:
        g.figure.axes[0].axhline(y=0, xmin=0, xmax=1, ls="-", color="black", zorder=10, lw=1)

def style_trial(g, color='k', alpha=0.3):
    if g.axes_dict: 
        for (name,ax) in g.axes_dict.items():
            ax.fill_between(ax.lines[0].get_data()[0], ax.lines[0].get_data()[1], 
                        color=color, alpha=alpha)
    else:
        ax = g.figure.axes[0]
        ax.fill_between(ax.lines[0].get_data()[0], ax.lines[0].get_data()[1], 
                        color=color, alpha=alpha)
  
def lineplot(data, x=None, y=None, errorbar=None, aspect=1.5, marker='o', 
             mec=None, ms=10, **kwargs):
    
    g = sns.relplot(data, x=x, y=y, kind="line",  aspect=aspect, 
                    errorbar=errorbar, marker=marker, mec=mec,ms=ms,
                    facet_kws={"sharey":False, "sharex":False}, **kwargs)
    return g

def barplot(data, x=None, y=None, hue=None, col=None, color=None, aspect=1,
            palette=None, hue_order=None, legend=None, col_wrap=None, 
            errorbar="se", dodge=False, capsize=0.2, sharex=False, sharey=False,
            native_scale=True, err_kws={"zorder":0.5, "lw":1}, **kwargs):
     
    g = sns.catplot(data, x=x, y=y, palette=palette, col_wrap=col_wrap,
                    hue_order=hue_order, hue=hue, color=color, col=col, 
                    kind="bar", errorbar=errorbar, dodge=dodge, capsize=capsize, 
                    err_kws=err_kws, sharex=sharex, sharey=sharey,aspect=aspect,
                    native_scale=native_scale, legend=legend, **kwargs)
    return g

def plot_hr(data, x="Time (hr)", y="lick", color='k',
            xlim=[-48, 24], ylim=[-10, 10], ymajmult=2, yminmult=1,
            ylabel="", xlabel="Time (hr)", hspace=0.45, 
            wspace=0.25, suptitle='', title=None, **kwargs):
    
    g = lineplot(data, x=x, y=y,color=color, **kwargs)
    
    style_hr(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
             suptitle=suptitle, title=title,ymajmult=ymajmult, yminmult=yminmult,)  
    
    return g 

def plot_ant_lickfreq(data, x="Time (hr)", y="lick", hue="stimulus", 
                      col="condition", palette=["green", "red"], 
                      hue_order=["stimulus", "blank"],
                      ylim=[0, 12], ylabel="Lick Frequency (Hz)", 
                      **kwargs):
    g = plot_hr(data, x=x, y=y, hue=hue, col=col, palette=palette, 
                hue_order=hue_order, ylim=ylim, ylabel=ylabel, **kwargs)
    
    return g
 
def plot_ant_perf(data, x="Time (hr)", y="lick", col="condition",ylabel=f"Performance\nL_s-L_b", **kwargs):
    g = plot_hr(data, x=x, y=y, col=col,ylabel=ylabel, **kwargs) 
    style_performance(g)
    return g

def plot_ant_perf_bar(data, x="Time (hr)", y="lick", hue=None, col="condition", 
                      xlim=[-50, 26], ylim=[-10, 10], palette=None, 
                      hue_order=None, color='k', ylabel="Performance", 
                      xlabel="Time (hr)", hspace=0.45, wspace=0.35, legend=True,
                      col_wrap=None, suptitle='', **kwargs):
    
    g = barplot(data, x=x, y=y, col=col, col_wrap=col_wrap, hue=hue, 
                palette=palette, hue_order=hue_order, color=color, 
                legend=legend, **kwargs)  
    style_hr(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
             suptitle=suptitle)
    style_performance(g)

    return g

def plot_lickfreq(data, x="Time (ms)", y="lick", hue="stimulus", 
                  col="Time (hr)",color=None,legend=True, col_wrap=6,
                  palette=["green", "red"], hue_order=["stimulus", "blank"], 
                  xlim=[-200, 2000], ylim=[0, 12], ylabel="Lick Frequency",
                  xlabel="Time (ms)", hspace=0.45, wspace=0.25, suptitle=None,
                  xmajmult=500, ymajmult=2, xminmult=100, 
                  yminmult=1, style=True,**kwargs):
    
    g = lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                 palette=palette, color=color, hue_order=hue_order, 
                 legend=legend, **kwargs)
    style_ms(g, xlabel, ylabel, xlim, ylim, xmajmult=xmajmult, ymajmult=ymajmult,
             xminmult=xminmult, yminmult=yminmult,hspace=hspace, wspace=wspace, 
             suptitle=suptitle,style=style)
    
    return g

def plot_perf(data, x="Time (ms)", y="lick", hue="Time (hr)", 
              col="condition", color='k', legend=True, col_wrap=None,
              palette=None, hue_order=None, xlim=[-200, 2000], ylim=[-10, 10],
              ylabel="Performance", xlabel="Time (ms)", hspace=0.45, 
              wspace=0.25, suptitle=None, ymajmult=2, yminmult=1,xmajmult=500, xminmult=100, **kwargs):
    g = lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                 palette=palette, color=color, hue_order=hue_order, 
                 legend=legend, **kwargs)
    style_ms(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, ymajmult=ymajmult,
             suptitle=suptitle, yminmult=yminmult, xmajmult=xmajmult, xminmult=xminmult)
    style_performance(g)
    
    return g

def plot_bar(data, x=None, y=None, aspect=1, col=None, 
                   color="#3B3838", palette=None, col_wrap=None, hue_order=None,
                   hue=None, xlim=[-0.5, 2.5], ylim=[0, 600], xlabel=None, hspace=0.45, wspace=0.25,
                   ylabel="", width=0.4, ymajmult=100, yminmult=50,title='', suptitle='',
                   **kwargs):
    g = barplot(data, x=x, y=y, aspect=aspect, col=col, color=color, 
                palette=palette, hue_order=hue_order, hue=hue, 
                col_wrap=col_wrap, width=width, **kwargs)
    style_axes(g, xlabel, ylabel, xlim, ylim, None, None, ymajmult=ymajmult, 
               yminmult=yminmult, hspace=hspace, wspace=wspace, title=title, suptitle=suptitle)

    return g

def plot_trial_hr(data, x="Time (hr)", y="trial no", col="condition", color='k', 
                  hue=None, palette=None, col_wrap=None, hue_order=None, 
                  errorbar='se', err_style='bars', xlabel="Time (hr)", 
                  ylabel="Number of trials", xlim=[-48, 48], ylim=[0, 300], 
                  hspace=0.45, wspace=0.25, suptitle=None, fill=True, title=None, 
                  ymajmult=100, yminmult=50,**kwargs):
    g = lineplot(data, x=x, y=y, hue=hue, col=col, color=color, palette=palette, 
                 hue_order=hue_order, errorbar=errorbar, err_style=err_style, 
                 **kwargs)
    style_hr(g, xlabel, ylabel, xlim, ylim, ymajmult=ymajmult, yminmult=yminmult,
             hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)
    # for name, ax in g.axes_dict.items():
    #     ax.spines[['right']].set_visible(True)
    #     ax.spines[['left']].set_visible(False)
    #     ax.xaxis.set_visible(False)
    #     ax.yaxis.tick_right()

    if fill: style_trial(g)
    return g

def connect_lines(ax,data, x=None,y=None, hue=None, order=None, hue_order=None):
    segs = []
    if hue and hue != x:
        for i in range(len(order)):
            hfilt = data[data[x] ==  order[i]]
            for a in hfilt["animal"].unique():
                anfilt = hfilt[hfilt["animal"] == a]
                pts = []
                for j in range(len(hue_order)):
                    op = list(anfilt[(anfilt[hue] == hue_order[j])][y])
                    if (len(op) != 0):
                        pts.append([i + 1/(len(hue_order)+1)*(j-(len(hue_order)-1)/2), op[0]])
                if(len(pts) != 0):
                
                    segs.append(pts)
    else:
        for a in data["animal"].unique():
            anfilt = data[data["animal"] == a]
            pts = []
            for j in range(len(order)):
                op = list(anfilt[(anfilt[x] == order[j])][y])
                if (len(op) != 0):
                        pts.append([j, op[0]])
            if(len(pts) != 0):
                segs.append(pts)
    linesegs = LineCollection(segs,color='k', lw=0.5, zorder=1)
    ax.add_collection(linesegs)

def plot_bar_strip(data, x=None, y=None, row=None, col=None, col_order=None, connect=False,
                   errorbar='se', color=None, pt_color=None, err_kws={'lw':1, 'zorder':0}, 
                   capsize=0.2, row_order=None,ylim=[0,10], xlim=[-0.5, 3.5], order=None, dodge=False, jitter=False,
                   hue=None, palette=None, hue_order=None, ylabel='', xlabel='',xlab_kws={}, **kwargs):
    g = sns.catplot(data, x=x, y=y, row=row,col=col,col_order=col_order,order=order,
                    errorbar=errorbar, color=color, err_kws=err_kws, capsize=capsize,kind='bar', 
                    row_order=row_order, hue=hue, palette=palette, hue_order=hue_order, **kwargs)
    if g.axes_dict:
        for name, ax in g.axes_dict.items():
            if type(name) == tuple:
                i, j = name
                cdfilt = data[(data[col] == j) & (data[row]==i)] 
                ax.set_title(f'{i} {j}')
            else:
                if col:
                    cdfilt = data[data[col] == name]
                elif row:
                    cdfilt = data[data[row] == name]
                ax.set_title(name)
            sns.stripplot(cdfilt, ax=ax, x=x, color=pt_color, y=y, 
                          dodge=dodge, jitter=jitter, hue=hue, legend=False,
                          hue_order=hue_order, palette=palette)
            if connect:
                connect_lines(ax,cdfilt, hue=hue, order=order, x=x, y=y,hue_order=hue_order)
                            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.tick_params(labelbottom=True, labeltop=False)
            ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)
    else:
        ax = g.axes[0][0]
        sns.stripplot(data, ax=ax, x=x, color=pt_color, y=y, 
                        dodge=dodge, jitter=jitter, hue=hue, legend=False,
                        hue_order=hue_order, palette=palette)
        if connect:
            connect_lines(ax,data, hue=hue, order=order, x=x,y=y, hue_order=hue_order)
                        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(labelbottom=True, labeltop=False)
        ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)

    g.tight_layout()
    return g

def draw_heatmap(*args, **kwargs):
    data = kwargs.pop('data')
    yind = kwargs.pop('yind')
    d = data.pivot(index=args[1], columns=args[0], values=args[2])
    for y in yind:
        if y not in d.index:
            d.loc[y] = np.NaN
    d = d.sort_index()
    sns.heatmap(d, **kwargs)

def row_col_facet(data, func_args, func, row, col, row_order=None, col_order=None, **kwargs):
    if col_order: 
            cols = col_order
    else:
        cols = data[col].unique()
    if row_order:
        rows = row_order
    else:
        rows = data[row].unique()
    fg, axs = plt.subplots(len(rows), len(cols), figsize=(5*(len(cols) + 1), 5*len(rows)))
    for i in range(len(rows)):
        cax=None
        for j in range(len(cols)):
            plt_data = data[(data[row] == rows[i]) & (data[col] == cols[j])]
            cbar = j == (len(cols)- 1)
            if cbar:
                cax = axs[i][len(cols) - 1].inset_axes([1.05, 0, 0.05,1])
            func(*func_args, data=plt_data,ax=axs[i][j],
                                cbar=cbar,cbar_ax=cax,**kwargs)                            
            axs[i][j].set_title(f'{rows[i]}, {cols[j]}')
    return fg

def col_facet(data, func_args, func, col, col_wrap=None,col_order=None, **kwargs):
    if col_order:
        cols = col_order
    else:
        cols = data[col].unique()
    if col_wrap:
        if len(cols) % col_wrap == 0:
            nrows = (len(cols) // col_wrap)
        else:
            nrows = (len(cols) // col_wrap) + 1
        fg, axs = plt.subplots(nrows,col_wrap, figsize=(5*col_wrap, 5*nrows))
    else:
        fg, axs = plt.subplots(1, len(cols), figsize=(5*len(cols), 5))
    axs = axs.flatten()
    cax=None
    for j in range(len(cols)):
        plt_data = data[(data[col] == cols[j])]
        cbar = j == (len(cols)-1)
        if cbar:
            cax = axs[len(cols) - 1].inset_axes([1.05, 0, 0.05,1])
        func(*func_args,data=plt_data,ax=axs[j],
                        cbar=cbar,cbar_ax=cax,**kwargs)                            
        axs[j].set_title(f'{cols[j]}')
    return fg

def row_facet(data, func_args, func, row,row_order=None, **kwargs):
    if row_order:
        rows = row_order
    else:
        rows = data[row].unique()    
    fg, axs = plt.subplots(len(rows), 1, figsize=(5, 5*len(rows)))
    for i in range(len(rows)):
        plt_data = data[(data[row] == rows[i])]
        cbar = True
        func(*func_args, data=plt_data,ax=axs[i],
                        cbar=cbar,**kwargs)                            
        axs[i].set_title(f'{rows[i]}')
    return fg

def single_facet(data, func_args, func,**kwargs):
    fg, ax = plt.subplots(1, 1, figsize=(5,5))
    cax = ax.inset_axes([1.05, 0, 0.05,1])
    func(*func_args, data=data,ax=ax,
                    cbar=True,cbar_ax=cax,**kwargs) 
    return fg

def draw_facet(data, func_args,  func, row=None, col=None, col_wrap=None, col_order=None, row_order=None,**kwargs):
    yind = data[func_args[1]].unique()
    if row and col:
        cols = col_order if col_order else data[col].unique()
        rows = row_order if row_order else data[row].unique()
        if len(rows) == 1 and len(cols) == 1:
            fg = single_facet(data, func_args, func,yind=yind, **kwargs)
            fg.axes[0][0].set_title(f'{rows[0]}, {cols[0]}')
        elif len(rows) == 1:
            fg = col_facet(data[data[row] == rows[0]], func_args, func, col, col_wrap=col_wrap, 
                           col_order=col_order,yind=yind,**kwargs)
        elif len(cols) == 1:
            fg = row_facet(data[data[col] == cols[0]], func_args, func, row, row_order,yind=yind, **kwargs)
        else:
            fg = row_col_facet(data, func_args, func, row, col,
                               row_order=row_order, col_order=col_order, yind=yind,**kwargs)
    elif row:
        rows = data[row].unique()
        if len(rows) == 1:
            fg = single_facet(data, func_args, func,yind=yind, **kwargs)
            fg.axes[0].set_title(f'{rows[0]}')
        else:
            fg = row_facet(data, func_args, func, row, row_order=row_order,yind=yind,**kwargs)
    elif col:
        cols = data[col].unique()
        if len(cols) == 1:
            fg = single_facet(data, func_args, func,yind=yind, **kwargs)
            fg.axes[0].set_title(f'{cols[0]}')
        else:
            fg = col_facet(data, func_args, func, col, col_wrap=col_wrap, 
                           col_order=col_order,yind=yind,**kwargs)
    else:
        fg = single_facet(data, func_args, func,yind=yind, **kwargs)
    fg.tight_layout()
    return fg

def plot_pair_strip(data, x=None, y=None, row=None, col=None, col_order=None, connect=False,
                   color=None, row_order=None,ylim=[0,10], xlim=[-0.5, 1.5], order=None, dodge=False, jitter=False,
                   hue=None, palette=None, hue_order=None, ylabel='', xlabel='',xlab_kws={}, aspect=0.5,legend=False,**kwargs):
    g = sns.catplot(data, y=y,x=x, order=order,kind="strip", hue=hue, 
                    palette=palette, hue_order=hue_order, col=col, row=row,
                    jitter=False, aspect=aspect, legend=legend,col_order=col_order, 
                    row_order=row_order, color=color, **kwargs)
    if g.axes_dict:
        for name, ax in g.axes_dict.items():
            if type(name) == tuple:
                cdfilt = data[(data[row] == name[0]) & (data[col] == name[1])]
                ax.set_title(f'{name[0]}, {name[1]}')
            else:
                if col:
                    cdfilt = data[data[col] == name]
                elif row:
                    cdfilt = data[data[row] == name]
            if connect:
                connect_lines(ax,cdfilt, hue=hue, order=order, x=x, y=y,hue_order=hue_order)                    
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.tick_params(labelbottom=True, labeltop=False)
            ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)
            ax.set_title(name)
    else:
        ax = g.axes[0][0]
        if connect:
                connect_lines(ax,cdfilt, hue=hue, order=order, x=x, hue_order=hue_order)                    
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(labelbottom=True, labeltop=False)
        ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)
        name  = col if col else ''
        name  =  name + row if row else name
        ax.set_title(name)
    return g

def style_retrain(g, xlim, ylim, xticklab, xlab, ylab):
    g.axes[0][0].set_xticklabels(xticklab)
    xticks = g.axes[0][0].xaxis.get_major_ticks()
    if xlim[1] > 48:
        xticks[8].set_visible(False)
    if xlim[1] > 96:
        xticks[14].set_visible(False)
    xticks = g.axes[0][0].xaxis.get_minor_ticks()
    vis = []
    if xlim[1] > 48:
        vis = [13, 14, 15, 16]
        g.axes[0][0].text(60, ylim[0]-0.06*(ylim[1]-ylim[0]), '//', ha='center', va='bottom', fontsize=24)
        g.axes[0][0].text(96, ylim[0]-.2*(ylim[1]-ylim[0]), '5-8d test', ha='center')
    if xlim[1] > 96:
        vis = vis + [25, 26, 27, 28]
        g.axes[0][0].text(144-12, ylim[0]-0.06*(ylim[1]-ylim[0]), '//', ha='center', va='bottom', fontsize=24)
        g.axes[0][0].text(168, ylim[0]-.2*(ylim[1]-ylim[0]), '21-24d test', ha='center')
    for i in vis:
            xticks[i].set_visible(False)
    g.axes[0][0].text(12, ylim[0]-.2*(ylim[1]-ylim[0]), 'Train', ha='center')
    g.axes[0][0].set_xlabel(xlab, labelpad=28)
    g.axes[0][0].set_ylabel(ylab)

def plot_retrain(c, means, perf, tot, xlim, xticklab,xlab='Time (hr)', lckylab='Lick Frequency (Hz)', 
                 perfylab="Performance\n$L_s - L_b$", totylab='Number of trials',
                 lckylim=[-0.1, 10], perfylim=[-6,6], totylim=[0, 400], aspect=2, **kwargs):
    g1 = None
    if not means.empty:
        g1 = plot_ant_lickfreq(means[(means['cond'] == c) & (means['test_time'] == 'train')],
                            y='lick',aspect=aspect,lw=1, ms=7,col='cond', title=c, 
                            legend=False, xlim=xlim, ylim=lckylim, errorbar='se', 
                            err_style='bars', err_kws={"lw":0.5},wspace=0.15,style=False, **kwargs)
        plt_data = means[(means['cond'] == c) & (means['test_time'] == '5d test')]
        if not plt_data.empty:
            sns.lineplot(plt_data,
                        y='lick',lw=1, ms=7,ax=g1.axes[0][0],x='Time (hr)',
                        errorbar='se', err_style='bars', err_kws={"lw":0.5},hue='stimulus', 
                        hue_order=['blank', 'stimulus'], palette=['red', 'green'], 
                        mec=None, marker='o', legend=False)
        plt_data = means[(means['cond'] == c) & (means['test_time'] == '21d test')]
        if not plt_data.empty:
            sns.lineplot(plt_data,
                        y='lick',lw=1, ms=7,ax=g1.axes[0][0],x='Time (hr)',
                        errorbar='se', err_style='bars', err_kws={"lw":0.5},hue='stimulus', 
                        hue_order=['blank', 'stimulus'], palette=['red', 'green'], 
                        mec=None, marker='o',legend=False)
        style_retrain(g1, xlim, lckylim, xticklab, xlab, lckylab)
    g2 = None
    if not perf.empty:
        g2 = plot_ant_perf(perf[(perf['cond'] == c) & (perf['test_time'] == 'train')],
                        ylabel=perfylab, aspect=aspect,lw=1, ms=7,
                        title='',xlim=xlim,col='condition',legend=True, ylim=perfylim, 
                        errorbar='se', err_style='bars', err_kws={"lw":0.5},
                        wspace=0.25,**kwargs)
        sns.lineplot(perf[(perf['cond'] == c) & (perf['test_time'] == '5d test')],
                    y='lick',lw=1, ms=7,ax=g2.axes[0][0],x='Time (hr)', errorbar='se', 
                    err_style='bars', err_kws={"lw":0.5}, mec=None, marker='o', 
                    color='k', legend=False)
        sns.lineplot(perf[(perf['cond'] == c) & (perf['test_time'] == '21d test')],
                    y='lick',lw=1, ms=7,ax=g2.axes[0][0],x='Time (hr)', errorbar='se', 
                    err_style='bars', err_kws={"lw":0.5},color='k', mec=None, 
                    marker='o', legend=False)
        style_retrain(g2, xlim, perfylim, xticklab,xlab, perfylab)
    g3 = None
    if not tot.empty:
        g3 = plot_trial_hr(tot[(tot['cond'] == c) & (tot['test_time'] == 'train')],
                        xlim=xlim, col='condition',fill=False,ylim=totylim, 
                        aspect=aspect*1.5, title='', ms=7, lw=1,err_kws={"lw":0.5}, 
                        ylabel=totylab,wspace=0.04,**kwargs)
        sns.lineplot(tot[(tot['cond'] == c) & (tot['test_time'] == '5d test')],
                    y='trial no',lw=1, ms=7,ax=g3.axes[0][0],x='Time (hr)', 
                    errorbar='se', err_style='bars', err_kws={"lw":0.5},color='k', 
                    mec=None, marker='o', legend=False)
        sns.lineplot(tot[(tot['cond'] == c) & (tot['test_time'] == '21d test')],
                    y='trial no',lw=1, ms=7,ax=g3.axes[0][0],x='Time (hr)',
                    errorbar='se', err_style='bars', err_kws={"lw":0.5},color='k', 
                    mec=None, marker='o', legend=False)
        style_retrain(g3, xlim, totylim, xticklab, xlab, totylab)
        for i in range(len(g3.axes[0][0].lines)):
            if i % 2 == 0:
                    g3.axes[0][0].fill_between(g3.axes[0][0].lines[i].get_data()[0], g3.axes[0][0].lines[i].get_data()[1], 
                                            color='k', alpha=0.3)
    return (g1, g2, g3)

#######################################
# Short Helper Functions for Analysis #
#######################################

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

####################################
# Longer Analysis Helper Functions #
####################################

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
    data = bin_by_time(data, time_bin, "hr", index, [], keep + [time, offset] + values, key_start,
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
    

    return (data_mean, counts_groups, perf)#, poke_perf)

############################
# Main Format and Analysis #
############################
def make_condition_df_helper(andirs, condir, condition, metadata, meta_vals, mssg, root, prog, proglab,name, i=0, animals=[]):
    if i == len(andirs):
        for a in animals:
            if not isinstance(a, pd.DataFrame):
                prog.after(1, (lambda : evaluate_format(None, mssg, proglab, prog,name)))
        animals = pd.concat(animals, ignore_index=True)
        animals["condition"] = condition
        prog.after(1, (lambda : evaluate_format(animals, mssg, proglab, prog,name)))
        return
    if i < len(andirs):
        animal_name = andirs[i]
        animal_path = condir + "\\" + animal_name
        animal = make_animal_df(animal_path, metadata, animal_name, meta_vals, mssg, root)
        animals.append(animal)
        root.after(1, lambda : make_condition_df_helper(andirs, condir, condition, metadata, meta_vals, mssg, root,prog, proglab, name,i + 1, animals))

def make_condition_df(condir, condition, metadata, meta_vals, mssg, root, prog, proglab, name):
    '''Load and format all data files for all animals in a condition, and return as a data frame.
    
    Requires all animals to be included to be in the same directory, which should have no other folders or files in it.
    
    Parameters:
    condir --- path of directory containing animal subdirectories which contain data files
    condition --- name of condition
    metadata --- data frame containing relevant metadata about animals in condition (can contain other information as well)
    meta_vals --- list of colums to include metadata from. Must include 'acc' for length of acclimation (in days).'''
    andirs = os.listdir(condir)
    if len(andirs) == 0:
        print('No animals in condition.')
        mssg.set('No animals in condition.')
        return pd.DataFrame()
    animals = []
    root.after(1, lambda : make_condition_df_helper(andirs, condir, condition, pd.read_excel(metadata), meta_vals, mssg, root,prog, proglab,name, i=0, animals=animals))

def run_animals(metadata, meta_df, files, root, cols, mssg, names, proglab, prog, name,i=0, animals=[]):
    if i == len(files):
        for a in animals:
            if not isinstance(a, pd.DataFrame):
                prog.after(1, (lambda : evaluate_format(None, mssg, proglab, prog,name)))
        animals = pd.concat(animals, ignore_index=True)
        if len(names) == 1:  
            animals['condition'] = names[0]
        prog.after(1, (lambda : evaluate_format(animals, mssg, proglab, prog,name)))
        return
        

    if len(metadata) > 1:
        meta_df = pd.read_excel(metadata[i])
    # name animal with directory name
    an_name = files[i].split('\\')[-1].split('/')[-1]
    animal = []
    animal = make_animal_df(files[i], meta_df, an_name, cols, mssg, root)
    print(animal)
    if len(names) != 1:
        animal['condition'] = names[i]
    animals = animals + [animal]
    root.after(1, lambda : run_animals(metadata,  meta_df, files, root, cols, mssg, names, proglab, prog, name, i+1, animals=animals))

def run_formatter(files, as_conditions, metadata, names, columns, mssg, root, prog, proglab, name):
    '''
    Returns formatted training data as a data frame, and writes that data to the
    output file if provided. Can format one or more condition directiories, or one or more animals.
    If multiple animals are being loaded at the same time and `names` is length 1,
    all animals are given the same condition name. If adding to a file of previously
    formatted data, will take the set of columns in the previously formatted data.

    Parameters:
    files --- list of animal or condition director(ies) to format
    as_conditions --- boolean indicating whether `files` is animal or condition directories
    output --- path to write formatted output to
    add --- boolean indicating whether to overwrite or add to `output` file
    metadata --- path to metadata file
    names --- list of names corresponding to provided files
    columns --- set of column names to load from metadata file'''
    
    mssg.set("Formatting")
    
    cols = ['acc']
    if columns:
        cols = cols + columns
    
    # treat provided list of files as condition directories containing animal subdirectories
    if as_conditions:
        if not len(files) == len(metadata):
            print(f'Number of condition directories ({len(files)}) and corresponding metadata files ({len(metadata)}) must match.')
            mssg.set(mssg.get() + f'\nNumber of condition directories ({len(files)}) and corresponding metadata files ({len(metadata)}) must match.')
            return
        if names:
            if not len(files) == len(names):
                print('Number of conditions and provided names must match.')
                mssg.set(mssg.get() + '\nNumber of conditions and provided names must match.')
                return
        conds = []
        root.after(1, lambda : make_condition_df(files[0], names[0], metadata[0], cols, mssg, root,prog, proglab, name))

    # treat provided list of files as animal directories containing training data
    else:
        if not names or len(names) != 1:
            print('Must provide exactly one name to load one condition')
            mssg.set(mssg.get() + '\nMust provide exactly one name to load one condition')
            return
        if not len(metadata) == 1 or not (len(metadata) == len(files)):
            print('Provide one metadata file for all animals, or number of animals and metadata files must match.')
            mssg.set(mssg.get() + '\nProvide one metadata file for all animals, or number of animals and metadata files must match.')
            return

        animals = []
        meta_df = pd.read_excel(metadata[0])
        root.after(1, lambda : run_animals(metadata,meta_df, files, root, cols, mssg, names,proglab, prog,name, i=0, animals=animals))

def format_wrapper(as_cond, mssg, prog, proglab, name, fname):
    run_formatter([files_var.get()], as_cond, 
                        [metadata_var.get()], 
                        [name], keep_var.get().split(', '),
                        mssg, root, prog, proglab, fname)
    
def evaluate_format(formatted_res, mssg, proglab, prog, name):
    if not isinstance(formatted_res, pd.DataFrame):
        mssg.set(mssg.get() + "\nFormatting failed!")
        proglab.config(foreground='red')
        prog.after(1, (lambda : root.after(2000, (lambda: prog.destroy()))))
    else:
        prog.after(1, (lambda : after_format_wrapper(formatted_res, prog, name)))

def after_format_wrapper(formatted_res, prog, name):
     
    mins = {
        "min_trials": int(min_var.get()), 
        "min_blank": int(min_blank_var.get()),   
        "min_stimulus": int(min_stimulus_var.get()) # minimum number of water trials 
    }

    index = ["condition", "animal"]
    full_keys = {
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
    
    keys = full_keys.copy()
    keys['millis bin'] = None
    index = ["condition", "animal"]
    keep = keep_var.get().split(', ') + ['acc', 'stimulus', 'water', 'type'] 
    
    full_suffix = ['data', 'means', 'means_filtered', 
                   'performance', 'performance_filtered', 
                   ]
    ant_suffix = full_suffix + ['trial counts', 'trials summed day'] 

    last20_suffix = ['data', 'means', 'performance', 'trial counts']

    analysis_output_path = add_file_label_var.get()
    if add_val.get():
        analysis_output = os.listdir(analysis_output_path)
        analysis_output.sort()
        full_out = [analysis_output_path + '\\' +  analysis_output[i] for i in range(7,12)]
        ant_out = [analysis_output_path + '\\' +  analysis_output[i] for i in range(0,7)]
        last20_out = [analysis_output_path + '\\' +  analysis_output[i] for i in range(12,16)]
    else:
        num, denom = num_val.get().split('/')
        full_out = [analysis_output_path + '\\inst_'  + name + '_' + suff + '.txt' for suff in full_suffix]
        ant_out = [analysis_output_path + f'\\fixed {start_val.get()}-{end_val.get()}_'  + name + '_' + suff + '.txt' for suff in ant_suffix]
        last20_out = [analysis_output_path + f'\\{num}-{denom} part _'  + name + '_' + suff  + '.txt' for suff in last20_suffix]
    out_files = full_out + ant_out + last20_out
    # ignore columns provided in 'keep' that are not present in data
    for c in keep:
        if c not in formatted_res.columns:
            keep.remove(c)  

    if len(values_var.get().split(', ')) == 1:
        vals = [values_var.get()]
    else:
        vals = values_var.get().split(', ')             

    tb = int(time_bin.get())

    root.after(1, lambda : inst_analyze_wrapper(formatted_res, full_keys, mins, keep, vals, index, tb, prog, keys, out_files))

def fixed_window_step1(data, time_bin, key_dict, keep, values, 
                      start, end, agg_arg):
    data = puff_delta(data, ['condition', 'animal', 'trial no'], "timestamp", "delay", "puff delta")
    root.after(1, lambda: fixed_window_step2(data, time_bin, key_dict, keep, values, 
                                     start, end, agg_arg))

def fixed_window_step2(data, time_bin, key_dict, keep, values, 
                    start, end, agg_arg):
    # using first sample of each trial to align to timebin and day aligns trial
    # to bin where trial started in edge case where trial spans bins and 
    # anticipatory period happens after bin end
    # matches Matlab v16 analysis behavior
    data = get_first_sample(data, ['condition', 'animal', 'trial no'], 'timestamp', "first sample")
    key_dict['timestamp'] = 'first sample'
    root.after(1, lambda: fixed_window_step3(data, time_bin, key_dict, keep, values,
                                     start, end, agg_arg))

def fixed_window_step3(data, time_bin, key_dict, keep,
                        values, start, end, agg_arg):
    # get lick frequency for each trial in given window
    data = fixed_window_lickfreq_helper(data, ['condition', 'animal', 'trial no'], values, 
                            keep + ['timestamp', 'first sample','offset'],'puff delta', start, end)
    root.after(1, lambda: fixed_window_step4(data, time_bin, key_dict, keep, values, agg_arg))

def fixed_window_step4(data, time_bin, key_dict, keep,
                        values, agg_arg):
    index_trial = ['condition', 'animal', 'trial no']
    # align trials to timebin and day
    data = align_to_timebin(data, time_bin, key_dict, keep, index_trial, values)
    data = align_to_day(data, key_dict, keep, index_trial, values)
    key_dict['timestamp'] = "timestamp"
    root.after(1, lambda: aggregate_values(data,time_bin, key_dict,['condition', 'animal'], values, keep, "bin", *agg_arg))

def fixed_window_lickfreq(data, window, time_bin, key_dict, keep, values,agg_arg):
    puff = key_dict["millis bin"]
    time = key_dict["timestamp"]
    index_trial = ["condition", "animal","trial no"]
    start, end = window
    key_start = "first sample"
    df = data.copy()
    root.after(1, lambda: fixed_window_step1(df, time_bin, key_dict, keep,
                                      values, start, end, agg_arg))

def inst_step1(df, step2_args, agg_args):
    # get time to airpuff delivery
    data = puff_delta(df, ['condition', 'animal','trial no'], "timestamp", "delay", "puff delta")
    
    # correct for pandas behavior when resampling with negative timedeltas
    # if the minimum delay (800 ms) is not present in the dataset
    l = len(data.index)
    if (data['puff delta'].min()) != pd.to_timedelta('-800 ms'):
        data.loc[l, 'trial no'] = 0
        data.loc[l, 'puff delta'] = pd.to_timedelta('-800ms')
    root.after(1, lambda: inst_step2(data, *step2_args, agg_args))

def inst_step2(data, time_bin, key_dict, values, keep, freq_bin, agg_args):
    # calculate instentanteous or rolling window lick frequency for each trial
    keep1 = keep.copy()
    index_trial = ['condition', 'animal', 'trial no']
    data = bin_by_time(data, freq_bin, "ms", index_trial, values, keep + ["timestamp", "offset"], "puff delta", origin=None, fn="lickfreq")
    # correct for pandas behavior when resampling with negative timedeltas
    # if the minimum delay (800 ms) is not present in the dataset
    data = data.loc[data['trial no'] > 0]
    root.after(1, lambda: inst_step3(data,time_bin, keep, values, key_dict, agg_args))

def inst_step3(data, time_bin, keep, values, key_dict, agg_args):
     # align sample times for cross-trial averaging
    index_trial = ['condition', 'animal', 'trial no']
    data = bin_by_time(data, 100, "ms", index_trial, [], keep + ["timestamp", "offset"] + values, "puff delta", origin=None)
    data = time_to_float(data, "Time (ms)", "puff delta", "ms")
    root.after(1, lambda: inst_step4(data, time_bin, key_dict, keep, values, agg_args))
    
def inst_step4(data,time_bin, key_dict, keep, values, agg_args):
    # align to timebin and day
    index_trial = ['condition', 'animal', 'trial no']
    data = align_to_timebin(data, time_bin, key_dict, keep + ["puff delta","Time (ms)"], index_trial, values)
    data = align_to_day(data, key_dict, keep + ["puff delta","Time (ms)"], index_trial, values)
    root.after(1, lambda: aggregate_values(data, time_bin, key_dict, ['condition', 'animal'],values, keep, "full", *agg_args))

def instentanteous_lickfreq(data, freq_bin, time_bin, key_dict, keep, values, agg_args):
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
    df = data.copy()
    step2_args = [time_bin, key_dict, values, keep, freq_bin]
    root.after(1, lambda: inst_step1(df, step2_args, agg_args))

    return data

def write_analysis(out_files, new_data, columns, add=False):
    mode = 'a' if add else 'w'
    
    for k in range(len(new_data)):
        if add:
            with open(out_files[k]) as f:
                old_cols = f.readline().strip('\n').split(',')
            out = new_data[k][old_cols]
            cols = old_cols
        else:
            out = new_data[k][columns[k]]
            cols = columns[k]
        out.to_csv(out_files[k], index=False, mode=mode, header=(not add), columns=cols)
        
def aggregate_values(data, tb, keys, index, values,keep, kind, out_files, mins, prog, add, df):
    all_cols = index + keep + ["Day", "Time (hr)"]
    if kind =='full':
        fs = [n for n in out_files if "inst" in n]
        all_cols = all_cols + ["Time (ms)"]
        keys["time bin"] = "delta"
    elif kind == 'bin':
        fs = [n for n in out_files if "fixed" in n]
        keys["time bin"] = "delta"
    elif kind == 'nth_part':
        fs = [n for n in out_files if "part" in n]
        keys["time bin"] = "day_delta"
    fs.sort()
    
    keep1 = list(data.columns) 
    for c in index + values + ['trial no']:
        keep1.remove(c)
    means, counts, perf = agg(data, keys, index, keep1, values)
    
    data_cols = all_cols + ['trial no', 'delta', 'day_delta'] + values 
    means_cols = all_cols + values
    perf_cols = means_cols.copy()
    for i in ['stimulus', 'water', 'type']:
        if i in perf_cols:
            perf_cols.remove(i)
    trial_cols = all_cols + ['trial no']
    if kind =='nth_part':
        if fs:
            columns = [data_cols, means_cols, perf_cols, trial_cols]
            write_analysis(fs,[data, means, perf, counts], columns, add=add)
        return (means, counts, perf)
    
    daykeep = list(counts.columns)
    for c in index + ['Day', 'trial no', 'timestamp', 'offset']:
        if c in daykeep:
            daykeep.remove(c)

    daytots = sum_trials(counts, index + ["Day"], daykeep, "trial no")
    data_filtered = drop_group(data, mins, keys, "trial no", index, keep1)
    filtmeans, filtcounts, filtperf = agg(data_filtered, keys, index, keep1, values)
    
    if kind == 'full':
        columns = [data_cols,means_cols, means_cols, perf_cols, perf_cols]        
        output = (data, means, filtmeans, perf,  filtperf)
    elif kind=='bin':
        columns = [data_cols,means_cols, means_cols, perf_cols, perf_cols, trial_cols, trial_cols]        
        output = (data, means, filtmeans, perf,  filtperf, counts, daytots)
    if fs: 
        write_analysis(fs,output, columns, add=add)
    if kind == "full":
        root.after(1, lambda : fixed_window_wrapper(df, keys, mins, keep, values, index, tb, prog, out_files))
    elif kind=='bin':
        root.after(1, lambda: last20_wrapper(data, keys, mins, keep, values, index, out_files, prog))
    elif kind=="nth_part":
        print('here')
        root.after(1, lambda: finish_analysis(prog))

def run_analysis(df, kind, keys, mins, keep, values, index, analysis_args, prog, out_files=None, add=False):
    agg_args = [out_files, mins, prog, add, df]
    if kind =='full':
        data = instentanteous_lickfreq(df, *analysis_args, keys, keep, values, agg_args) #freq_window, freq_bin, time_bin
    elif kind == 'bin':
        data = fixed_window_lickfreq(df, *analysis_args, keys, keep, values, agg_args) #r, time_bin
    elif kind == 'nth_part':
        data = get_nth_part_day(df, *analysis_args) #num_parts, nth_part
        aggregate_values(data, 0, keys, index, values, keep, "nth_part", *agg_args)
    else:
        print(f'Unknown analysis type {kind}.')
    
def inst_analyze_wrapper(formatted_res, full_keys, mins, keep, vals, index, tb, prog, keys, out_files):
    if inst_anal.get():
        mssg.set("Instentaneous analysis")
        run_analysis(formatted_res, 'full', 
                            full_keys, mins, keep, vals, index, 
                            [100, tb], prog, out_files, add_val.get())
    else:
        root.after(1, lambda : fixed_window_wrapper(formatted_res, keys, mins, keep, vals, index, tb, prog, out_files))

def fixed_window_wrapper(formatted_res, keys, mins, keep, vals, index, tb, prog, out_files):
    if fixed_anal.get():
        mssg.set("Fixed window analysis")
        analysis_window = (int(start_val.get()), int(end_val.get()))
        keys['millis bin'] = None
        ant_analysis = run_analysis(formatted_res, 'bin', 
                        keys, mins, keep, vals, index, [analysis_window, tb],prog, out_files, add_val.get())
    else:
        print("here2")
        root.after(1, lambda: finish_analysis(prog))

def last20_wrapper(ant_analysis, keys, mins, keep, vals, index, out_files, prog):
    num, denom = num_val.get().split('/')
    if nth_frac_val.get():
        mssg.set("Last 20% analysis")
        last20 = run_analysis(ant_analysis, 'nth_part', 
                        keys, mins, keep, vals, index, [int(num), int(denom)], prog, out_files, add_val.get())
    else:
        root.after(1, lambda: finish_analysis(prog))


def close_and_plot(prog):
    if analyplot.get():
        root.after(2000, prog.destroy)
        generate_plots()

def finish_analysis(prog):
    mssg.set("Done!")
    root.after(1, lambda: close_and_plot(prog))


mssg = StringVar(value="Analysis")

def analysis():
    if files_var.get() == '':
        files_errmsg.set("Please select data to analyze.")
        files_err.grid(row=1, column=0, padx=defaultpad, pady=defaultpad, sticky='nsw')
        return
    if metadata_var.get() == '':
        metadata_errmsg.set("Please select metadata file.")
        metadata_err.grid(row=3, column=0, padx=defaultpad, pady=defaultpad, sticky='nsw')
        return
    if add_file_label_var.get() == '':
        add_file_errmsg.set('Please select output file.')
        add_file_err.grid(row=1, column=0, padx=defaultpad, pady=defaultpad, sticky='nsw')
        return
    if add_val.get():
        if not validate_previous_analysis():
            return
    prog = tk.Toplevel(root)
    prog.title('Analysis Progress')
    prog.geometry('500x500+100+100')
    proglab = ttk.Label(prog, textvariable=mssg)

    proglab.grid(row=0, column=0, padx=2, pady=2)

    as_cond = (type_var.get()=="Condition")
    name = Path(files_var.get()).parts[-1]
    
    fname = name_var.get()
    if fname == '':
        fname=name    
       
    prog.after(1, lambda : format_wrapper(as_cond, mssg, prog, proglab, name, fname))    


def validate_previous_analysis():
    analermsg.set('')
    flermsg.set('')
    analerlab.grid_forget()
    flerlab.grid_forget()
    prev_anal_files = os.listdir(add_file_label_var.get())
    prev_anal_files.sort()
        
    if len(prev_anal_files) == 0:
        add_file_errmsg.set("No previous analysis.")
        add_file_err.grid(row=1, column=0, padx=defaultpad, pady=defaultpad, sticky='nsw')
        return
    
    fixed = [n for n in prev_anal_files if 'fixed' in n]
    nth = [n for n in prev_anal_files if 'part' in n]
    full = [n for n in prev_anal_files if 'inst' in n]

    fixed_str = 'fixed, ' if fixed_anal.get() else ''
    inst_str = 'inst, ' if inst_anal.get() else ''
    prev_inst_str = 'inst, ' if len(full) > 0 else ''
    prev_fixed_str = 'fixed, ' if len(fixed) > 0 else ''
    nth_str = 'part' if len(fixed) > 0 else ''
    prev_nth_str = 'part' if len(nth) > 0 else ''

    if fixed_str != prev_fixed_str or inst_str != prev_inst_str or nth_str != prev_nth_str:
        analermsg.set('Current analysis selection (%s %s %s) does not match previous analysis (%s %s %s)' 
                    % (inst_str,fixed_str, nth_str, prev_inst_str, prev_fixed_str, prev_nth_str))
    
    fixed_start, fixed_end = fixed[0].split('_')[0].split(' ')[-1].split('-')

    if fixed_start != start_val.get():
        analermsg.set(analermsg.get() + f'\nCurrent fixed start ({start_val.get()}) does not match previous fixed start ({fixed_start})')
    if fixed_end != end_val.get():
        analermsg.set(analermsg.get() + f'\nCurrent fixed end ({end_val.get()}) does not match previous fixed end ({fixed_end})')
    
    num, denom = nth[0].split(' ')[0].split('-')
    cnum, cdenom = num_val.get().split('/')
    if num != cnum or denom != cdenom:
        analermsg.set(analermsg.get() +f'\nCurrent frac day ({cnum}/{cdenom}) does not match previous frac day ({num}/{denom})')
    
    testfile = pd.read_csv(add_file_label_var.get() + '\\' + prev_anal_files[-1])

    prevcols = list(testfile.columns.drop(['condition', 'animal', 'acc', 'stimulus', 'water', 'type', 'Day', 'Time (hr)', 'Time (ms)', 'trial no', 'lick', 'poke','delta', 'day_delta'], errors='ignore'))
    currcols = keep_var.get().split(', ')
    if Counter(prevcols) != Counter(currcols):
        flermsg.set(f"Current columns to keep ({', '.join(currcols)}) do not match previous columns to keep ({', '.join(prevcols)})")

    prevtim = testfile.groupby(['condition', 'animal', 'Time (hr)']).first().reset_index()['Time (hr)'].diff(-1).abs()[0]
    currtim = float(time_bin.get())
    if not math.isclose(prevtim, currtim):
        flermsg.set(flermsg.get() + f'\nCurrent time bin ({currtim}) does not match previous timebin ({prevtim})')

    prevval = list(testfile.columns.drop(['condition', 'animal', 'acc', 'stimulus', 'water', 'type', 'Day', 'Time (hr)', 'trial no', 'Time (ms)', 'delta', 'day_delta'] + prevcols, errors='ignore'))
    currval = values_var.get().split(', ')
    if Counter(prevval) != Counter(currval):
        flermsg.set(flermsg.get() + f"\nCurrent values ({', '.join(currval)}) do not match previous values ({', '.join(prevval)})")

    if flermsg.get() != '':
        flerlab.grid(row=0, column=0, padx=defaultpad, pady=defaultpad, sticky='nsew')
    if analermsg.get() != '':
        analerlab.grid(row=0, column=0, padx=defaultpad, pady=defaultpad, sticky='nsew')
        print(analermsg.get())
    if flermsg.get() == '' and analermsg.get() == '':
        return True
    else:
        return False

# Function to open a file dialog
def open_file(file_path_var, file_err):
    file_path = filedialog.askopenfilename()
    file_path_var.set(file_path)
    file_err.grid_forget()

def open_directory(file_path_var, file_err):
    file_path = filedialog.askdirectory()
    file_path_var.set(file_path)
    file_err.grid_forget()
    
def enableFixedWindow():
    if fixed_anal.get():
        fixed_anal_start_en.state(['!disabled'])
        fixed_anal_end_en.state(['!disabled'])
        nthfrac_cbutt.state(['!disabled'])
        en_num.state(['!disabled'])
    else:
        fixed_anal_start_en.state(['disabled'])
        fixed_anal_end_en.state(['disabled'])
        nthfrac_cbutt.state(['disabled'])
        en_num.state(['disabled'])

def enableAddFolder():
    if add_val.get():
        add_file_butt.state(['!disabled'])
    else:
        add_file_butt.state(['disabled'])

def enableNthFrac():
    if nth_frac_val.get():
        en_num.state(['!disabled'])
    else:
        en_num.state(['disabled'])

# adapted from tkinter docs
errmsg = StringVar()
formatmsg = 'Start and end should be numbers between -800 and 2000'
def check_num(newval):
    valid = re.match('^-?[0-9]*$', newval) is not None and int(newval) > -800 and int(newval) < 2000
    if not valid:
        fixed_anal_err.grid(row=5, column=0,padx=defaultpad, pady=defaultpad,sticky='nsw')
        errmsg.set(formatmsg)
    else:
        fixed_anal_err.grid_forget()
    return valid

check_num_wrapper = (root.register(check_num), '%P')
# example tooltip
# tooltip_qEPSCConvert = ToolTip(ConvertButton, "Convert .txt files to .csv")

# Widgets for analysis
analysis_lf = ResizeEqualLabelFrame(root, rows=5, cols=1,text="1. Analysis")

files_lf = ResizeEqualLabelFrame(analysis_lf, rows=2, cols=8,text='Select files for analysis')
flermsg = StringVar()
flerlab = WrappingLabel(files_lf, textvariable=flermsg, foreground='red')
files_var = StringVar(value='')
files_butt = EnterButton(files_lf, text='Select folder(s)', command=(lambda : open_directory(files_var, files_err)))
files_lab = WrappingLabel(files_lf, textvariable=files_var)
#files_lab.bind('<Configure>', lambda e: files_lab.config(width=files_lf.winfo_width()))
files_errmsg = StringVar()
files_err = WrappingLabel(files_lf, textvariable=files_errmsg, foreground='red')


metadata_var = StringVar(value='')
metadata_butt = EnterButton(files_lf, text='Select metadata file', command=(lambda : open_file(metadata_var, metadata_err)))
metadata_lab = WrappingLabel(files_lf, textvariable=metadata_var, width=50)
metadata_errmsg = StringVar()
metadata_err = WrappingLabel(files_lf, textvariable=metadata_errmsg, foreground='red')

formatimemsg = 'Time bin must be a whole number (min 1)'
errtimemsg = StringVar()
def check_whole_time(newval):
    valid = re.match('^[0-9]*$', newval) is not None and int(newval) > 0
    if not valid:
        time_err.grid(row=6, column=0,padx=defoutpad, pady=defoutpad,sticky='nsw')
        errtimemsg.set(formatimemsg)
    else:
        time_err.grid_forget()
    return valid

check_whole_time_wrapper = (root.register(check_whole_time), '%P')

time_bin = StringVar(value='4')
time_en = ttk.Entry(files_lf, textvariable=time_bin,validate='focusout', validatecommand=check_whole_time_wrapper)
time_lab = WrappingLabel(files_lf,text='Time bin (hours):')
time_err = WrappingLabel(files_lf, textvariable=errtimemsg, foreground='red')

keep_var = StringVar(value='Age, Cage, Sex')
keep_en = ttk.Entry(files_lf, textvariable=keep_var)
keep_lab = WrappingLabel(files_lf, text='Metadata columns to keep:')

# TODO: validate values in 'lick', 'poke'
values_var = StringVar(value='lick')
values_en = ttk.Entry(files_lf, textvariable=values_var)
values_lab = WrappingLabel(files_lf,text='Values to analyze:')

analy_lf = ResizeEqualLabelFrame(analysis_lf,rows=2, cols=10, text='Analysis Types:')

analermsg = StringVar()
analerlab = WrappingLabel(analy_lf, textvariable=analermsg, foreground='red')
inst_anal = BooleanVar(value=False)
inst_cbutt = EnterCheckbutton(analy_lf, text='Instentaneous analysis',variable=inst_anal, onvalue=True, offvalue=False)

fixed_anal = BooleanVar(value=True)
fixed_cbutt = EnterCheckbutton(analy_lf, text='Fixed window (anticipatory) analysis',variable=fixed_anal, onvalue=True, offvalue=False, command=enableFixedWindow)
start_val = StringVar(value='700')
startlabel = WrappingLabel(analy_lf, text='Start of fixed window:')
fixed_anal_start_en = ttk.Entry(analy_lf,textvariable=start_val, validate='focusout', validatecommand=check_num_wrapper)

end_val = StringVar(value='1000')
endlabel = WrappingLabel(analy_lf, text='End of fixed window:')

fixed_anal_end_en = ttk.Entry(analy_lf,textvariable=end_val, validate='focusout', validatecommand=check_num_wrapper)
fixed_anal_err = WrappingLabel(analy_lf, textvariable=errmsg, foreground='red')

nth_frac_val = BooleanVar(value=True)
nthfrac_cbutt = EnterCheckbutton(analy_lf, text='Fraction of day analysis (default last 20%)', variable=nth_frac_val, onvalue=True, offvalue=False, command=enableNthFrac)

nthfrac_errmsg = StringVar()
def check_frac(newval):
    valid = re.match('^[0-9]/[0-9]$', newval) is not None
    if not valid:
        nth_frac_err.grid(row=9, column=0,padx=defaultpad, pady=defaultpad,sticky='nsw')
        nthfrac_errmsg.set('Numerator and denomenator must be digits.')
        return valid
    num = newval.split('/')[0]
    denom = newval.split('/')[1]
    valid = valid and int(num) <= int(denom)
    if not valid:
        nth_frac_err.grid(row=9, column=0,padx=defaultpad, pady=defaultpad,sticky='nsw')
        nthfrac_errmsg.set('Numerator must be less than denomenator.')
    else:
        nth_frac_err.grid_forget()
    return valid
check_frac_wrapper = (root.register(check_frac), '%P')

num_val = StringVar(value='5/5')
en_num = ttk.Entry(analy_lf, textvariable=num_val, validate='focusout', validatecommand=check_frac_wrapper)
nth_frac_lab = WrappingLabel(analy_lf, text='Fraction of day:')
nth_frac_err = WrappingLabel(analy_lf, textvariable=nthfrac_errmsg, foreground='red')

formattrialmsg = 'Trial thresholds must be a whole number.'
errtrialmesg = StringVar(value='')
def check_whole(newval):
    valid = re.match('^[0-9]*$', newval) is not None and int(newval) > 0
    if not valid:
        errtrialmesg.set(formattrialmsg)
        min_err.grid(row=3, column=0,padx=defaultpad, pady=defaultpad,sticky='nsw')
    else: 
        min_err.grid_forget()
    return valid

check_whole_wrapper = (root.register(check_whole), '%P')

trial_thresh_lf = ResizeEqualLabelFrame(analysis_lf, rows=2, cols=6,text='Minimum trial thresholds')
min_var = StringVar(value='10')
min_en = ttk.Entry(trial_thresh_lf, textvariable=min_var, validate='focusout', validatecommand=check_whole_wrapper)
min_lab = WrappingLabel(trial_thresh_lf, text='Minimum number of trials:')

min_blank_var = StringVar(value='1')
min_bl_en = ttk.Entry(trial_thresh_lf, textvariable=min_blank_var, validate='focusout', validatecommand=check_whole_wrapper)
min_bl_lab = WrappingLabel(trial_thresh_lf, text='Minimum number of blank trials:')

min_stimulus_var = StringVar(value='1')
min_st_en = ttk.Entry(trial_thresh_lf, textvariable=min_stimulus_var, validate='focusout', validatecommand=check_whole_wrapper)
min_st_lab = WrappingLabel(trial_thresh_lf, text='Minimum number of stimulus trials:')
min_err = WrappingLabel(trial_thresh_lf, textvariable=errtrialmesg, foreground='red')


add_lf = ResizeEqualLabelFrame(analysis_lf, cols=2, rows=3,text='Output folder')

add_file_label_var = StringVar(value='')
add_file_errmsg = StringVar()
add_file_err = WrappingLabel(add_lf, textvariable=add_file_errmsg, foreground='red')
add_file_butt = EnterButton(add_lf, text='Select output folder',command=(lambda : open_directory(add_file_label_var, add_file_err)))
add_file_label = WrappingLabel(add_lf, textvariable=add_file_label_var, width=50)

add_val = BooleanVar(value=False)
add_cbutt = EnterCheckbutton(add_lf, text='Add to previous analysis',variable=add_val, onvalue=True, offvalue=False)


type_lf = ResizeEqualLabelFrame(analysis_lf,rows=2, cols=4, text='File type')

def updateName():
    name_lab_var.set(f'{type_var.get()} name (optional):')

type_var = StringVar(value='Animal')
animal_rad = ttk.Radiobutton(type_lf, text='Animal', variable=type_var, value='Animal', command=updateName)
condition_rad = ttk.Radiobutton(type_lf, text='Condition', variable=type_var, value='Condition',command=updateName)

name_lab_var = StringVar()
updateName()
name_lab = WrappingLabel(type_lf, textvariable=name_lab_var)
name_var = StringVar()
name_en = ttk.Entry(type_lf, textvariable=name_var)

run_anal_butt = ttk.Button(analysis_lf, text='Analyze', command=analysis)

dtype_map = {'condition':str,'animal':str, 'stimulus':str, 'lick':float, 'Strain':str, 'Cage':str, 'type':str, 'Time (ms)':float, 'Time (hr)':float, 'Day':float, 'water':str}
# Widgets for Plotting
def get_means(df, gpcols):
    means = df.groupby(gpcols)['lick'].mean()
    kpcols = df.columns.drop(gpcols + ['lick'])
    kp  = df.groupby(gpcols)[kpcols].first()
    means = pd.concat([means, kp], axis=1).reset_index()
    return means

def generate_ant_perf(data, pltwindow):
    if ant_perf_var.get():
        try:
            perf = data['fixed performance filtered']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='Fixed Performance: No fixed analysis')
            errlab.grid()
            return
        perf['Time (hr)'] = perf['Time (hr)'] + 2
        xlim = [perf['Time (hr)'].min() - 2, perf['Time (hr)'].max() + 2]
        g = plot_ant_perf(perf, ylabel="Performance", aspect=2, xlim=xlim,ms=7, ylim=[-10, 10])
        g.figure.subplots_adjust(left=0.15, bottom=0.15)
        if dispvar.get():
            perfcan = mplEmbedPlot(g.figure, master=pltwindow)
            pltwindow.add(perfcan.frame, text='Fixed Performance')
        if savevar.get():
            g.savefig(savefilevar.get()+"\\fixed_perf.png", transparent=True)

def generate_ant_lf(data, pltwindow):
    if ant_lf_var.get():
        try:
            ant_means = data['fixed means filtered']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No fixed analysis')
            errlab.grid()
            return
        ant_means['Time (hr)'] = ant_means['Time (hr)'] + 2
        xlim = [ant_means['Time (hr)'].min() - 2, ant_means['Time (hr)'].max() + 2]
        g = plot_ant_lickfreq(ant_means, xlim=xlim, aspect=2, ms=7)
        g.figure.subplots_adjust(left=0.1, bottom=0.15)
        if dispvar.get():
            lfcan = mplEmbedPlot(g.figure, master=pltwindow)
            pltwindow.add(lfcan.frame, text='Fixed Lick Frequency')
        if savevar.get():
            g.savefig(savefilevar.get()+"\\fixed_lickfreq.png", transparent=True)

def generate_fullpfhr(data, pltwindow):
    if fullpfhr.get():
        try:
            perf = data['inst performance filtered']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No instentaneous analysis')
            errlab.grid()
            return
            return
        for c in perf['condition'].unique():
            g = plot_perf(perf[perf['condition'] == c], xlim=[-700, 2000], 
                          col="Time (hr)", col_wrap=6, hue=None, color='k',suptitle=c,
                          ms=2, lw=0.7, ymajmult=5, yminmult=2.5, xminmult=250, wspace=0.2,hspace=0.2 )
            g.tick_params(labelsize=8, pad=0.5)
            i=0
            nrow = math.ceil(len(g.axes)/6)
            for n, ax in g.axes_dict.items():
                ax.set_xlabel(ax.get_xlabel(), fontsize=8)
                ax.set_ylabel(ax.get_ylabel(), fontsize=8)
                ax.set_title(f'Time (hr):{ax.get_title()}', fontsize=8)
                if i % 6 != 0:
                    ax.set_ylabel('')
                if i < len(g.axes)-6:
                    ax.set_xlabel('')
                ax.tick_params(which='major', axis='x', labelrotation=45)
                i +=1
            g.figure.set_size_inches(6*1.5, nrow*1.5)
            g.tight_layout()
            g.figure.subplots_adjust(left=0.05, bottom=0.15)
            if dispvar.get():
                flperfcan = mplEmbedPlot(g.figure, master=pltwindow)
                pltwindow.add(flperfcan.frame, text='Instentaneous Performance')
            if savevar.get():
                g.savefig(savefilevar.get()+f"\\perffull_{c}.png", transparent=True)

def genereate_fullckhr(data, pltwindow):
     if fulllckhr.get():
        try:
            means_full = data['inst means filtered']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No instentaneous analysis')
            errlab.grid()
            return
        for c in means_full['condition'].unique():
            g = plot_lickfreq(means_full[means_full['condition'] == c], suptitle=c,
                              xlim=[-700, 2000], col="Time (hr)",col_wrap=6, ms=2, lw=0.7, 
                              ymajmult=5, yminmult=2.5, xminmult=250, wspace=0.2, hspace=0.2)
            g.tick_params(labelsize=8, pad=0.5)
            i=0
            for n, ax in g.axes_dict.items():
                ax.set_xlabel(ax.get_xlabel(), fontsize=8)
                ax.set_ylabel(ax.get_ylabel(), fontsize=8)
                ax.set_title(f'Time (hr):{ax.get_title()}', fontsize=8)
                if i % 6 != 0:
                    ax.set_ylabel('')
                if i < len(g.axes)-6:
                    ax.set_xlabel('')
                ax.tick_params(axis='x', which='major', labelrotation=45)
                i +=1
            nrow = math.ceil(len(g.axes)/6)
            g.figure.set_size_inches(6*1.5+1, nrow*1.5)
            g.legend.set_in_layout(True)
            g.figure.tight_layout()

            g.legend.set_title("Stimulus Type", prop=mpl.font_manager.FontProperties(size=8))
            for t in g.legend.get_texts():
                t.set_size(8)
            g.legend.set(loc='center right',bbox_to_anchor=(1.01, 0.5))
            g.figure.subplots_adjust(left=0.05, bottom=0.15,right=0.9)
            
            if dispvar.get():
                lfdisp = mplEmbedPlot(g.figure, master=pltwindow)
                pltwindow.add(lfdisp.frame, text='Instentaneous Lick Freq')
            if savevar.get():
                g.savefig(savefilevar.get()+f"\\lckfreqfull_{c}.png", transparent=True)

def generate_htmppf(data, pltwindow):
    if htmppf.get():
        try:
            data['inst performance filtered']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No instentaneous analysis')
            errlab.grid()
            return
        gpcols = ['condition', 'Time (hr)', 'Time (ms)']
        perf_full_means = get_means(data['inst performance filtered'], gpcols)
        htmp_data = perf_full_means.sort_values(['Time (hr)', 'Time (ms)'])
        htmp_data = htmp_data[(htmp_data['Time (ms)'] < 2000) & (htmp_data['Time (ms)'] > -700)]
        htmp_args=['Time (ms)', 'Time (hr)', 'lick']
        fg = draw_facet(htmp_data, htmp_args, draw_heatmap, col='condition',square=True,
                cmap='coolwarm', vmin=-6, vmax=6, col_wrap=6)
        fg.suptitle('Instentaneous performance', y=0.98)
        fg.figure.subplots_adjust(left=0.15, bottom=0.2, top=0.9)
        if dispvar.get():
            lfdisp = mplEmbedPlot(fg.figure, master=pltwindow)
            pltwindow.add(lfdisp.frame, text='Inst Perf Heatmap')
        if savevar.get():
            fg.savefig(savefilevar.get()+"\\perfheatmap.png", transparent=True)

def generate_htmplf(data, pltwindow):
    if htmplck.get():
        try:
            data['inst means filtered']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No instentaneous analysis')
            errlab.grid()
            return
        gpcols = ['condition', 'Time (hr)', 'Time (ms)','stimulus']
        means_full_means = get_means(data['inst means filtered'], gpcols)
        htmp_data = means_full_means.sort_values(['Time (hr)', 'Time (ms)'])
        htmp_data = htmp_data[(htmp_data['Time (ms)'] < 2000) & (htmp_data['Time (ms)'] > -700)]
        htmp_args=['Time (ms)', 'Time (hr)', 'lick']
        fg = draw_facet(htmp_data, htmp_args, draw_heatmap, col='stimulus', row='condition',
                col_order=['stimulus', 'blank'],square=True,
                cmap='viridis', vmin=0, vmax=10)
        fg.suptitle('Instentaneous licking', y=0.98)
        fg.figure.subplots_adjust(left=0.1, bottom=0.2, top=0.9)
        if dispvar.get():
            fdisp = mplEmbedPlot(fg.figure, master=pltwindow)
            pltwindow.add(fdisp.frame, text='Inst Lick Freq Heatmap')
        if savevar.get():
            fg.savefig(savefilevar.get()+"\\lckfreqheatmap.png", transparent=True)

def generate_nthpart(data, pltwindow):
    if last20val.get():
        try:
            last20 = data['part means']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No part analysis')
            errlab.grid()
            return
        gpcols = ['condition', 'Time (hr)', 'Time (ms)','stimulus']
        last20["Day"] = last20["Day"].apply(day_to_label)
        g = plot_pair_strip(last20, y="lick",x="stimulus", order=["stimulus", "blank"],
                        hue="stimulus", palette=["green", "red"], ylabel='Lick frequency (Hz)',
                        hue_order=["stimulus", "blank"],
                        row="condition", col="Day",jitter=False, s=7, aspect=0.5, legend=False, connect=True)
        g.figure.subplots_adjust(left=0.2, bottom=0.05)
        for (r, c), ax in g.axes_dict.items():
            ax.set_title(f'{r} {c}', fontsize=10)
        if dispvar.get():
            fdisp = mplEmbedPlot(g.figure, master=pltwindow)
            pltwindow.add(fdisp.frame, text='Nth Part')
        if savevar.get():
            g.savefig(savefilevar.get()+"\\nthpart.png", transparent=True)

def generate_trialshr(data, pltwindow):
      if trialsbyhour.get():
        try:
            trials = data['fixed trial counts']
        except KeyError:
            errdisp = tk.Toplevel()
            errlab = tk.Label(errdisp, text='No fixed analysis')
            errlab.grid()
            return
        gpcols = ['condition', 'Time (hr)', 'Time (ms)','stimulus']
        tots = sum_trials(trials, ['condition', 'animal', 'Time (hr)'],[], "trial no")
        xlim = [tots['Time (hr)'].min() - 2, tots['Time (hr)'].max() + 2]
        ymax = math.ceil(tots['trial no'].max()/100.0)*100
        g = plot_trial_hr(tots, xlim=xlim,ylim=[0, ymax], title=None, ms=7, lw=1,err_kws={"lw":0.5}, ylabel='Trial no')
        g.figure.subplots_adjust(left=0.15, bottom=0.15)
        if dispvar.get():
            fdisp = mplEmbedPlot(g.figure, master=pltwindow)
            pltwindow.add(fdisp.frame, text='Trials by hour')
        if savevar.get():
            g.savefig(savefilevar.get()+"\\trialsbyhour.png", transparent=True)

def generate_plots():
    f = loadfilevar.get()
    if f == '':
        loadfileerrmsg.set("Please select data to plot.")
        loadfileerr.grid(row=1,column=1)
        return
    if savevar.get() and savefilevar.get() == '':
        savefileerrmsg.set("Please select a folder to save plots.")
        savefileerr.grid(row=2, column=0)
        return
    data = dict()
    for type_name in os.listdir(f):
        splt_name = type_name.split('.')[0].split('_')
        splt_name = [re.search('[a-z]+',splt_name[0])[0]] + splt_name[2:]
        data[' '.join(splt_name)] = pd.read_csv(f + '\\' + type_name, dtype=dtype_map, header=0)

    plt_wind = tk.Toplevel()
    plt_window = ttk.Notebook(plt_wind)
    plt_window.pack()
    generate_ant_perf(data,plt_window)
    generate_ant_lf(data,plt_window)
    generate_fullpfhr(data,plt_window)
    genereate_fullckhr(data,plt_window)
    generate_htmplf(data,plt_window)
    generate_htmppf(data,plt_window)
    generate_nthpart(data,plt_window)
    generate_trialshr(data,plt_window)

analyplot = BooleanVar(value=False)
def analyze_and_plot():
    analyplot.set(True)
    loadfilevar.set(add_file_label_var.get())
    analysis()
    return

plot_lf = ResizeEqualLabelFrame(root,cols=1, rows=5, text="2. Plotting")

load_lf = ResizeEqualLabelFrame(plot_lf,cols=2, rows=6, text='Data to plot')
loadfilevar = StringVar()
loadfileerrmsg = StringVar()
loadfileerr = WrappingLabel(load_lf, textvariable=loadfileerrmsg, foreground='red')
loadbutt = EnterButton(load_lf, text='Select Folder', command=(lambda :open_directory(loadfilevar, loadfileerr)))
loadfile = WrappingLabel(load_lf,textvariable=loadfilevar)

ant_lf = ResizeEqualLabelFrame(plot_lf, rows=2, cols=1,text='Fixed Window Plots')
ant_lf_var = BooleanVar(value=True)
ant_lf_cb = EnterCheckbutton(ant_lf, text='Fixed window lick frequency', variable=ant_lf_var,onvalue=True, offvalue=False)
ant_perf_var = BooleanVar(value=True)
ant_perf = EnterCheckbutton(ant_lf, text='Fixed window performance', variable=ant_perf_var,onvalue=True, offvalue=False)

analplot = BooleanVar(value=False)

full_lf = ResizeEqualLabelFrame(plot_lf,rows=4, cols=1, text='Instentaneous Plots')

fulllckhr = BooleanVar(value=True)
fulllckhrcb = EnterCheckbutton(full_lf, text='Instentaneous lick frequency by hour', variable=fulllckhr,onvalue=True, offvalue=False)
fullpfhr = BooleanVar(value=True)
fullpfhrcb = EnterCheckbutton(full_lf, text='Instentaneous performance by hour', variable=fullpfhr,onvalue=True, offvalue=False)
htmplck = BooleanVar(value=True)
htmplckcb = EnterCheckbutton(full_lf, text='Instentaneous lick frequency heatmap', variable=htmplck,onvalue=True, offvalue=False)
htmppf = BooleanVar(value=True)
htmppfcb = EnterCheckbutton(full_lf, text='Instentaneous performance heatmap', variable=htmppf,onvalue=True, offvalue=False)

trial_lf = ResizeEqualLabelFrame(plot_lf,rows=3, cols=1, text='Trial Plots')
last20val = BooleanVar(value=True)
last20cb = EnterCheckbutton(trial_lf, text='Nth fraction of day stimulus vs blank', variable=last20val,onvalue=True, offvalue=False)
trialsbyhour = BooleanVar(value=True)
trialsbyhourcb = EnterCheckbutton(trial_lf, text='Trials by hour', variable=trialsbyhour, onvalue=True, offvalue=False)

other_lf = ResizeEqualLabelFrame(plot_lf,cols=1,rows=1, text='Other Plots')

save_lf = ResizeEqualLabelFrame(plot_lf,cols=2, rows=4, text='Save Options')
savevar = BooleanVar(value=True)
savecb = EnterCheckbutton(save_lf, text='Save plots', variable=savevar,onvalue=True, offvalue=False, command=(lambda:disablesavebutt()))

savefilevar = StringVar()
savefileerrmsg = StringVar()
savefileerr = WrappingLabel(save_lf, textvariable=savefileerrmsg, foreground='red')
savebutt = EnterButton(save_lf, text='Select Folder', command=(lambda :open_directory(savefilevar, savefileerr)))
savelab = WrappingLabel(save_lf, textvariable=savefilevar, width=25)
def disablesavebutt():
    if savevar.get():
        savebutt.state(['!disabled'])
    else:
        savebutt.state(['disabled'])

dispvar = BooleanVar(value=True)
dispcb = EnterCheckbutton(save_lf, text='Display plots', variable=dispvar,onvalue=True, offvalue=False)

run_plots_butt = EnterButton(plot_lf, text='Generate plots', command=generate_plots)

run_all_butt = EnterButton(root, text='Analyze and plot', command=analyze_and_plot)

#scb = ttk.Scrollbar(root, orient='vertical', command=root.yview)
#scb.grid

defaultpad=2
defoutpad=5
# grid plot buttons

plot_lf.grid(row=0,column=1, padx=defoutpad, pady=defoutpad,sticky='nsew')

load_lf.grid(row=0,column=0, padx=defoutpad, pady=defoutpad,sticky='ew')
loadbutt.grid(row=0,column=0, padx=defaultpad, pady=defaultpad,sticky='nsew')
loadfile.grid(row=0,column=1, padx=defaultpad, pady=defaultpad,sticky='nsew')

ant_lf.grid(row=1,column=0, padx=defoutpad, pady=defoutpad,sticky='nsew')
ant_perf.grid(row=1,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')
ant_lf_cb.grid(row=0,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')

full_lf.grid(row=2,column=0, padx=defoutpad, pady=defoutpad,sticky='nsew')
fullpfhrcb.grid(row=1,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')
fulllckhrcb.grid(row=0,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')
htmppfcb.grid(row=3,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')
htmplckcb.grid(row=2,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')

trial_lf.grid(row=3,column=0, padx=defoutpad, pady=defoutpad,sticky='nsew')
last20cb.grid(row=0,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')
trialsbyhourcb.grid(row=1,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')

#other_lf.grid()

save_lf.grid(row=6,column=0, padx=defoutpad, pady=defoutpad,sticky='nsew')
savebutt.grid(row=1,column=0, padx=defaultpad, pady=defaultpad,sticky='nsew')
savelab.grid(row=1,column=1, padx=defaultpad, pady=defaultpad,sticky='nsew')
savecb.grid(row=0,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')
dispcb.grid(row=3,column=0, padx=defaultpad, pady=defaultpad,sticky='nsw')

run_plots_butt.grid(row=7,column=0, padx=defoutpad, pady=defoutpad,sticky='nsew')

# grid buttons for analysis column
analysis_lf.grid(row=0, column=0, padx=defoutpad, pady=defoutpad,sticky='ewns')
files_lf.grid(row=0, column=0, padx=defaultpad, pady=defaultpad,sticky='ewns')

files_butt.grid(row=0, column=0, sticky='nsew',padx=defaultpad, pady=defaultpad)
files_lab.grid(row=0, column=1, sticky='ew', padx=defaultpad, pady=defaultpad)

metadata_butt.grid(row=2, column=0, sticky='nsew',padx=defaultpad, pady=defaultpad)
metadata_lab.grid(row=2, column=1, sticky='ew',padx=defaultpad, pady=defaultpad)

keep_en.grid(row=7, column=1, sticky='w',padx=defaultpad, pady=defaultpad)
keep_lab.grid(row=7, column=0,sticky='ewns',padx=defaultpad, pady=defaultpad)

values_en.grid(row=8, column=1, sticky='w',padx=defaultpad, pady=defaultpad)
values_lab.grid(row=8, column=0, sticky='ewns',padx=defaultpad, pady=defaultpad)

time_en.grid(row=5, column=1, padx=defoutpad, pady=defoutpad,sticky='w')
time_lab.grid(row=5, column=0, padx=defoutpad, pady=defoutpad,sticky='ewns')

analy_lf.grid(row=3, column=0, padx=defoutpad, pady=defoutpad,sticky='nsew')

inst_cbutt.grid(row=1, column=0,sticky='nsew',padx=defaultpad, pady=defaultpad, columnspan=2)

fixed_cbutt.grid(row=2, column=0, sticky='nsew',padx=defaultpad,pady=defaultpad, columnspan=2)
startlabel.grid(row=3, column=0,padx=[50,defaultpad], pady=defaultpad,sticky='ewns')
fixed_anal_start_en.grid(row=3, column=1,padx=defaultpad, pady=defaultpad,sticky='w')
endlabel.grid(row=4, column=0,padx=[50,defaultpad], pady=defaultpad,sticky='nsew')
fixed_anal_end_en.grid(row=4, column=1,padx=defaultpad, pady=defaultpad,sticky='nsw')

nthfrac_cbutt.grid(row=6, column=0, columnspan=3,padx=defaultpad, pady=defaultpad,sticky='nsew')
en_num.grid(row=7, column=1,padx=defaultpad, pady=defaultpad,sticky='w')
nth_frac_lab.grid(row=7, column=0,padx=[50,defaultpad], pady=defaultpad,sticky='ewns')

trial_thresh_lf.grid(row=4, column=0,padx=defoutpad, pady=defoutpad,sticky='w')
min_en.grid(row=0, column=1,padx=defaultpad, pady=defaultpad,sticky='w')
min_lab.grid(row=0, column=0,padx=defaultpad, pady=defaultpad,sticky='ewns')
min_bl_en.grid(row=1, column=1,padx=defaultpad, pady=defaultpad,sticky='w')
min_bl_lab.grid(row=1, column=0,padx=defaultpad, pady=defaultpad,sticky='ewns')
min_st_en.grid(row=2, column=1,padx=defaultpad, pady=defaultpad,sticky='w')
min_st_lab.grid(row=2, column=0,padx=defaultpad, pady=defaultpad,sticky='ewns')

add_lf.grid(row=5, column=0,padx=defoutpad, pady=defoutpad,sticky='w')
add_cbutt.grid(row=2, column=0, columnspan=2,padx=defaultpad, pady=defaultpad,sticky='ew')
add_file_butt.grid(row=0, column=0,padx=defaultpad, pady=defaultpad,sticky='nsew')
add_file_label.grid(row=0, column=1,padx=defaultpad, pady=defaultpad,sticky='nsew')

type_lf.grid(row=6, column=0,padx=defaultpad, pady=defaultpad, sticky='w')
animal_rad.grid(row=0, column=0, padx=[50,defaultpad], pady=defaultpad,sticky='nsew')
condition_rad.grid(row=1, column=0,padx=[50,defaultpad], pady=defaultpad,sticky='nsew')
name_lab.grid(row=2, column=0, padx=[50,defaultpad], pady=defaultpad,sticky='nse')
name_en.grid(row=2, column=1, padx=defaultpad, pady=defaultpad,sticky='nsew')

run_anal_butt.grid(row=20, column=0,padx=defoutpad, pady=defoutpad)

run_all_butt.grid(row=1, column=0, padx=defoutpad, pady=defoutpad)
for i in range(0, 6):
    if i < 2:
        root.grid_columnconfigure(i, weight=1)
        files_lf.grid_columnconfigure(i, weight=1)
    if i < 3:
        root.grid_rowconfigure(i, weight=1)
    analy_lf.grid_rowconfigure(i, weight=1)
    analy_lf.grid_columnconfigure(i, weight=1)
    files_lf.grid_rowconfigure(i, weight=1)
if __name__ == '__main__':
    # Start the main loop
    root.mainloop()
