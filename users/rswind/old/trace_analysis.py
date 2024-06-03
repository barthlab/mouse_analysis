#!python3
# last edit 09/05/2023 Rachel Swindell

import warnings
#to suppress seaborn palette warnings
warnings.filterwarnings("ignore", category=UserWarning)
#to suppress seaborn error estimation NaN warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import os
import argparse
import sys

import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter.simpledialog import Dialog
from tkinter import messagebox

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import loader
import analysis_functions

pd.options.mode.chained_assignment = None  # default='warn'

# TODO: output number of blank and water trials in an epoch in addition to the total number of trials
# TODO: implement last 20% analysis
# TODO: implement inter-trial interval analysis and epoch analysis

# TODO: add GUI that allows script to be run again after 1st run if run with GUI (not if run from command line)
# TODO: add option to output by-animal data?

# TODO: acc time in hours or days?

# adapted from QueryDialog in simpledialog.py from tkinter v 8.6
class MultiPromptDialog(Dialog):
    def __init__(self, title, prompts, initalvalues, types, parent = None):
        # TODO: error if length of prompt, initalvalues and types not the same

        self.prompts = prompts
        self.initialvalues = initalvalues
        self.types = types
        Dialog.__init__(self, parent, title)

    def body(self, master):
        self.entries = []
        for i in range(len(self.prompts)): 
            w = Label(master, text=self.prompts[i], justify=LEFT)
            w.grid(row=i, padx=5, sticky=W)

            self.entries.append(Entry(master, name=f"entry{i}"))
            self.entries[i].grid(row=i, column=1,padx=5, sticky=E)

            if self.initialvalues[i] is not None:
                self.entries[i].insert(0, self.initialvalues[i])
                self.entries[i].select_range(0, END)

        return self.entries[0]
    
    def apply(self):
        self.result = []
        for i in range(len(self.types)):
            try:
                self.result.append(self.types[i](self.entries[i].get()))
            except:
                print("Invalid type. Please input a value with type {self.types[i]}")
        return self.result


def select_directory(prompt):
    directory = filedialog.askdirectory(title=prompt)
    if directory:
        return directory
    else:
        sys.exit("No directory selected")

def get_parameters(title, prompts, initalvalues, types):
    d = MultiPromptDialog(title, prompts, initalvalues, types)
    return d.result
 
# adapted from code written by Alex Kiroff
def parse_args():
    '''Read, interpret, and return console commands.'''
    default_acclimation = "ACC days"
    default_bin_time = 4  # hours
    default_trial_bin = 300 # ms
    default_last_x_percent = 20 # percent
    default_trials = 0
    default_water_trials = 1
    default_blank_trials = 1
    default_condition = "SAT1"

    parser = argparse.ArgumentParser(description="Advanced mouse behavior analysis")
    parser.add_argument("condition_path", type=str, nargs="?", help="Path to directory with data files")
    parser.add_argument("output_path", type=str, nargs="?", help="Path to directory to output analysis to")
    parser.add_argument("--noUI", action="store_true", help="Run from command line, if present")
    parser.add_argument("-r", "--rolling", action="store_true", help="Use rolling window instead of licking frequency bins, if present")
    parser.add_argument("-n", "--condition_name", type=str, default=default_condition, help=f"Name of condition type (default: {default_condition})")
    parser.add_argument("-t", "--trial_bin", type=float, default=default_trial_bin, help=f"Licking frequency bin size in ms (default: {default_trial_bin})")
    parser.add_argument("-m", "--metadata", type=str, help="Path to file with animal metadata")
    parser.add_argument("-a", "--acclimation_time", type=str, default=default_acclimation, help=f"Name of column in metadata file with acclimation times (default: {default_acclimation})")
    parser.add_argument("-b", "--bin_time", type=float, default=default_bin_time, help=f"Bin time in hours (default: {default_bin_time})")
    parser.add_argument("-p", "--last_x_percent", type=float, default=default_last_x_percent, help=f"Percentage of data to consider by day (0 < last_x_percent < 100) (default: {default_last_x_percent})")
    parser.add_argument("-i", "--min_trials", type=float, default=default_bin_time, help=f"Minimum number of trials in a bin (default: {default_trials})")
    parser.add_argument("-w", "--min_water_trials", type=float, default=default_bin_time, help=f"Minimum number of water trials in a bin (default: {default_water_trials})")
    parser.add_argument("-k", "--min_blank_trials", type=float, default=default_bin_time, help=f"Minimum number of blank trials in a bin (default: {default_blank_trials})")
    parsed = parser.parse_args()

    # Check if all csv_files are valid files
    '''for csv_file in parsed.csv_files:
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
    '''
    return parsed

# adapted from code written by Alex Kiroff
def get_user_input():
    '''Read, interpret, and return parameters provided by user.'''
    csv_directory = select_directory("Select directory with CSV files")
    print(f'Directory with behavior files: {csv_directory}')

    analysis_directory = select_directory("Select directory for analysis files")
    print(f'Directory to output analysis to: {analysis_directory}')

    metadata_file = filedialog.askopenfilename()
    if metadata_file == '':
        print(f'No metatdata file provided')
        sys.exit(1)
    print(f'Metadata file: {metadata_file}')

    # condition_name, acc_col_name
    title = "Training Metadata"
    prompts = ["Condition Name", "Acclimation Column Name (in metadata file)"]
    initialvalues = ["SAT", "ACC days"]
    types = [str, str]
    metadata_params = get_parameters(title, prompts, initialvalues, types)
    if metadata_params == None:
        print(f'No metadata information provided')
        sys.exit(1)
    print(f'Condition: {metadata_params[0]}, Acclimation column name: {metadata_params[1]}')

    # bin_time_hours, last_x_percent, freq_window, freq_bin
    title = "Bin Size Parameters"
    prompts = ["Bin size (min)", "Frequency Rolling Average window size (ms)", "Frequency Bin Size (ms)", "Last percent of day"]
    initialvalues = [240, 300, 100, 20]
    types = [int, int, int, int]
    bin_params = []
    while len(bin_params) != len(initialvalues):
        bin_params = get_parameters("Input integers: " + title, prompts, initialvalues, types)
        if bin_params == None:
            print(f'No bin size information provided.')
            sys.exit(1)
    print(f'bin size: {bin_params[0]} min, rolling window: {bin_params[1]} ms, freq bin: {bin_params[2]} ms, last percent: {bin_params[3]}')

    # min_trials, min_water_trials, min_blank_trials
    title = "Trial number thresholds per bin"
    prompts = ["Minimum total trials", "Minimum water trials", "Minimum blank trials"]
    initialvalues = [10, 1, 1]
    types = [int, int, int]
    min_trials_nos = []
    while len(min_trials_nos) != len(initialvalues):
        min_trials_nos = get_parameters("Input integers: " + title, prompts, initialvalues, types)
        if min_trials_nos == None:
            print(f'No minimum trial information provided.')
            sys.exit(1)
    print(f'min number of trials: {min_trials_nos[0]}, min water trials: {min_trials_nos[1]}, min blank trials: {min_trials_nos[2]}')

    return csv_directory, analysis_directory, metadata_file, metadata_params, bin_params, min_trials_nos

def generate_performance_plot(perf_data):
    '''Plot average performance (stimulus - blank) and return the figure.
    
    Creates an individual plot for each condition. For a condition, plots average performance (stimulus licking frequecny - blank licking frequency)
    vs the trial time for each time bin in the training time.'''
    # plot all timebins average performance trace on the same plot
    g = sns.relplot(data=perf_data,kind="line", x="Time (ms)", y="lick",col="condition", 
                    hue="Time (hr)", palette="coolwarm", errorbar="se",err_style="bars", legend="full")

    # add lines at air puff and water delivery
    for ax in g.axes.flat:
        ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
        ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
        ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)
        ax.set_ylim([-10, 10])
        ax.set_ylabel("Performance")
    return g.fig

def generate_lickfreq_plot(data):
    '''Plot average lick frequency for stimulus and blank trials and return the figure.
     
      Creates indiviudal plots for each time bin across the training time. Plots average frequency values vs the trial time for each time bin.'''
    g = sns.relplot(data=data,kind="line",x="Time (ms)", y="lick", col="Time (hr)",col_wrap=7,
                    hue="stimulus", palette=["green", "red"], hue_order=["water", "blank"], errorbar="se",err_style="bars", legend="full")

    # add lines at air puff and water delivery
    for ax in g.axes.flat:
        ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
        ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
        ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)
        ax.set_ylim([0, 12])
        ax.set_ylabel("Lick Freq (hz)")
    return g.fig

if __name__ == "__main__":
    default_acc_time = 2 # Animals are assumed to be acclimated for 2 days if no acclimation time is provided
    
    args = parse_args()
    if args.noUI:
        csv_directory = args.condition_path
        analysis_directory = args.output_path
        metadata_file = args.metadata
        condition_name = args.condition_name
        acc_col_name = args.acclimation_time
        time_bin = args.bin_time * 60 # convert to mins

        if args.rolling:
            freq_window = args.trial_bin
            freq_bin = 100

        else:
            freq_window = 100
            freq_bin = args.trial_bin

        last_percent = args.last_x_percent
        min_trials = args.min_trials
        min_water_trials = args.min_water_trials
        min_blank_trials = args.min_blank_trials

    else:
        csv_directory, analysis_directory, metadata_file, metadata_params, bin_params, min_trials_nos = get_user_input()
        condition_name, acc_col_name = metadata_params[0], metadata_params[1]
        time_bin, freq_window, freq_bin, last_percent = bin_params[0], bin_params[1], bin_params[2], bin_params[3]
        min_trials, min_water_trials, min_blank_trials = min_trials_nos[0], min_trials_nos[1], min_trials_nos[2]

    if not args.noUI: 
        print("Loading files...")

    metadata = pd.read_excel(metadata_file)
    df = loader.make_condition_df(csv_directory, condition_name, metadata, acc_col_name, default_acc_time)

    
    if not args.noUI: 
        print("Calculating Statistics... ", end='')

    data, prev_blank, prev_water = lickfreq_analysis(df, freq_window, freq_bin, time_bin, (not args.noUI))
    mean_statistics, counts, performance = aggregate_analysis(data, min_trials, min_blank_trials, min_water_trials)
    # prev_blank_stats, prev_blank_counts, prev_blank_perf = aggregate_analysis(prev_blank, min_trials, min_blank_trials, min_water_trials)
    # prev_water_stats, prev_water_counts, prev_water_perf = aggregate_analysis(prev_water, min_trials, min_blank_trials, min_water_trials)
    
    if not args.noUI: print("Writing data to csv...\n", end='')

    output_dir = f'{analysis_directory}/{condition_name}/stats/{condition_name}'
    output_prev = f'{analysis_directory}/{condition_name}/prev_stats/{condition_name}'

    cols = ["condition", "sex", "age", "strain", "animal", "stimulus", "Time (hr)", "Time (ms)", "lick"]
    mean_statistics.to_csv(f'{output_dir}_lick_frequency.csv', columns=cols, index=False)
    prev_blank_stats.to_csv(f'{output_prev}_prevblank_lick_frequency.csv', columns=cols, index=False)
    prev_water_stats.to_csv(f'{output_prev}_prevwater_lick_frequency.csv', columns=cols, index=False)

    cols = ["condition", "sex", "age", "strain", "animal", "Time (hr)", "Time (ms)", "trial no"]
    counts.to_csv(f'{output_dir}_trial_counts.csv', columns=cols, index=False)
    prev_blank_counts.to_csv(f'{output_prev}_prevblank_trial_counts.csv', columns=cols, index=False)
    prev_water_counts.to_csv(f'{output_prev}_prevwater_trial_counts.csv', columns=cols, index=False)

    cols = ["condition", "sex", "age", "strain", "animal", "Time (hr)", "Time (ms)", "lick"]
    performance.to_csv(f'{output_dir}_performance.csv', columns=cols, index=False)
    prev_blank_perf.to_csv(f'{output_prev}_prevblank_performance.csv', columns=cols, index=False)
    prev_water_perf.to_csv(f'{output_prev}_prevwater_performance.csv', columns=cols, index=False)

    if not args.noUI: print("Making plots... ", end='', flush=True)

    trialcount_fig = generate_trialcount_plot(counts)
    trialcount_fig.savefig(f'{output_dir}_num_trials.png')

    if not args.noUI: print("... ", end='', flush=True)

    perf_fig = generate_performance_plot(performance)
    perf_fig.savefig(f'{output_dir}_performance.png')

    if not args.noUI: print("... ", end='', flush=True)

    lickfreq_fig = generate_lickfreq_plot(mean_statistics)
    lickfreq_fig.savefig(f'{output_dir}_lickfreq.png')

    if not args.noUI: print("...\nDone.\n", end='')