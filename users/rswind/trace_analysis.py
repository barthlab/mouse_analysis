#!python3
# last edit 7/6/2023 Rachel Swindell

import re
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

# to run:
# python3 trace_analysis.py --metadata metadata-file --trial_bin_size bin-size --bin_size bin-size [--group] path-to-data
# maintain both user-input and console-run options

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
                print(self.initialvalues[i])
                self.entries[i].select_range(0, END)

        return self.entries[0]
    
    def apply(self):
        self.result = []
        for i in range(len(self.types)):
            self.result.append(self.types[i](self.entries[i].get()))
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


# adapted from code in mouse_analysis.py written by Alex Kiroff
# TODO: fix to the parameters my script need
def get_user_input():

    csv_directory = select_directory("Select directory with CSV files")

    analysis_directory = select_directory("Select directory for analysis files")

    metadata_file = filedialog.askopenfile() 

    # condition_name, acc_col_name
    title = "Training Metadata"
    prompts = ["Condition Name", "Acclimation Column Name (in metadata file)"]
    initialvalues = [None, None]
    types = [str, str]
    metadata_params = get_parameters(title, prompts, initialvalues, types)

    # bin_time_hours, last_x_percent, freq_window, freq_bin
    title = "Bin Size Parameters"
    prompts = ["Bin size (hours)", "Frequency Rolling Average window size (ms)", "Frequency Bin Size (ms)", "Last percent of day"]
    initialvalues = [4, 300, 100, 20]
    types = [int, int, int, int]
    bin_params = get_parameters(title, prompts, initialvalues, types)

    # min_trials, min_water_trials, min_blank_trials
    title = "Trial number thresholds per bin"
    prompts = ["Minimum total trials", "Minimum water trials", "Minimum blank trials"]
    initialvalues = [10, 1, 1]
    types = [int, int, int]
    min_trials_nos = get_parameters(title, prompts, initialvalues, types)

    return csv_directory, analysis_directory, metadata_file, metadata_params, bin_params, min_trials_nos

# adapted from code in mouse_analysis.py written by Alex Kiroff
# TODO: fix to the data/file structure my script needs
def save_data_to_csv(final_statistics, acclimation_time, bin_time, file_path):
    # Column labels
    labels = ["Time Bin", "Avg Water Lick Freq (Hz)", "Avg Blank Lick Freq (Hz)", "Performance", "Num Datapoints"]

    # Open the CSV file for writing
    with open(file_path, 'w+', newline='') as csv_file:
        writer = csv.writer(csv_file)

        # Write column labels
        writer.writerow(labels)

        # Write data rows
        row_time = -1 * acclimation_time + bin_time / 2
        for row in final_statistics:
            avg_lick_frq_water, avg_lick_frq_blank, performance, num_datapoints = row
            data_row = [row_time, avg_lick_frq_water, avg_lick_frq_blank, performance, num_datapoints]
            writer.writerow(data_row)
            row_time += bin_time


# analysis: traces over entire trial -> rolling window with bin size flexibiliy 
#           bin over training time -> bin size flexibility
#           drop trials -> number of trials flexibility
def calculate_statistics(data, freq_window, freq_bin, time_bin, min_trials, min_water, min_blank):
    # do not modify loaded data
    aa = data.copy()

    # calculated time to air puff
    index = ["condition", "animal","trial no"]
    data = analysis_functions.puff_delta(aa, index)

    # rolling window average licking frequency
    # TODO: (moderately) slow (26s on 4 conditions, 10+ animals per condition)
    keep = ["timestamp", "age", "sex", "strain", "acc", "delivery delta"]
    values = ["lick", "poke"]
    index = ["condition", "animal", "trial no", "trial type"]
    data = analysis_functions.rolling_frequency_average(data, freq_window, values, keep, index)

    # resampling for trial alignment
    keep = ["timestamp", "lick", "poke", "age", "sex", "strain", "acc"]
    index = ["condition", "animal", "trial no", "trial type"]
    data = analysis_functions.resample_align(data, freq_bin, "delivery delta", keep, index)

    #label start of trial at each sample
    index = ["condition", "animal", "trial no"]
    data = analysis_functions.get_trial_start(data, index, "timestamp")

    #bin by time (4h bins)
    index = ["condition", "animal", "trial no", "trial type", "delivery delta"]
    data = analysis_functions.bin_by_time(data, time_bin, "h", index, "trial time")

    # calulate time bin relative to start of SAT training
    index = ["condition", "animal", "trial no"]
    data = analysis_functions.delta(data, index, "trial time")

    # drop bins with fewer than 10 total trials, or with no trials of one kind
    index = ["condition", "animal", "delta","delivery delta"]
    key = "trial type"
    data = analysis_functions.drop_bins(data, min_trials, min_blank, min_water, index, key)

    #convert time bins and trial time to float representations
    data["Time (hr)"] = data["delta"].astype("timedelta64[h]")
    data["Time (ms)"] = data["delivery delta"].astype("timedelta64[ms]")

    # number of trials for each timebin
    index = ["condition", "animal", "delta", "delivery delta"]
    keep = ["age", "sex", "strain", "acc", "Time (hr)", "Time (ms)"]
    counts = analysis_functions.trial_counts(data, index, keep, "trial no")
    
    # mean licking frequencys for each animals for each time bin 
    keep = ["age", "sex", "strain", "acc","Time (hr)", "Time (ms)"]
    index = ["condition", "animal", "trial type", "delta", "delivery delta"]
    value = ["lick"]
    data_mean = analysis_functions.mean_bin(data_mean, index, value, keep)

    #threshold trials to 200ms before to 2000ms after air puff
    data_mean = analysis_functions.thresh(data_mean, -200, 2000)
    
    #calculate performance for each animal for each time bin
    index = ["condition", "animal", "delta", "delivery delta"]
    keep = ["age", "sex", "strain", "Time (hr)", "Time (ms)"]
    perf = analysis_functions.performance(data, index, keep)

    return (data, data_mean, counts, perf)

def generate_trialcount_plot(counts_data):
    cond = (counts_data["Time (hr)"] < 24) & (counts_data["Time (hr)"] > -24) 
    g = sns.catplot(counts_data[cond], x="Time (hr)", y="trial no", col="condition", kind='bar', color="grey")
    for ax in g.axes.flat:
        ax.set_ylabel("Number of Trials")
        ax.set_ylim([0, 200])

def generate_performance_plots(perf_data):
### Plotting
    # average performance trace across all timebins
    # only plot last day of acclimation and first day of SAT
    cond = (perf_data["Time (hr)"] < 24) & (perf_data["Time (hr)"] > -24)

    # plot all timebins average performance trace on the same plot
    g = sns.relplot(data=perf_data[cond],kind="line", x="Time (ms)", y="lick",col="condition", 
                    hue="Time (hr)", palette="coolwarm", errorbar="se",err_style="bars", legend="full")

    # add lines at air puff and water delivery
    for ax in g.axes.flat:
        ax.axhline(y=0, xmin=0, xmax=1, ls="-", lw=0.75,color="black", zorder=0)
        ax.axvline(x=0, ymin=0, ymax=1, ls="--", color="lightgrey", zorder=0)
        ax.axvline(x=1000, ymin=0, ymax=1, ls="--", color="navy", alpha=0.5, zorder=0)
        ax.set_ylim([-10, 10])
        ax.set_ylabel("Performance")

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
    # average licking frequency by timebin across all animals
    index = ["condition", "delta", "delivery delta", "trial type"]
    keep = ["age", "sex", "strain", "Time (hr)", "Time (ms)"]

def generate_licking_plots(data):
    index = ["condition", "delta", "delivery delta", "trial type"]
    keep = ["age", "sex", "strain", "Time (hr)", "Time (ms)"]

    gp = data.groupby(index)
    lick_avg = gp["lick"].mean()
    keep = gp[keep].first()
    lick_avg = pd.concat([lick_avg, keep], axis=1).reset_index()

    # time to postive performance and magniutude of maximum licking frequency by trial type and timebin
    cond = (lick_avg["Time (hr)"] < 20) & (lick_avg["Time (hr)"] > -24) & (lick_avg["Time (ms)"] > 0) & (lick_avg["Time (ms)"] < 1000)
    pt_puff = lick_avg[cond].sort_values("lick").groupby(["condition", "delta", "trial type"]).first().reset_index()

    sns.catplot(data=pt_puff, kind="bar", x="condition", y="Time (ms)", col="Time (hr)", col_wrap=5, hue="trial type")
    sns.catplot(data=pt_puff, kind="bar", x="condition", y="lick", col="Time (hr)", col_wrap=5, hue="trial type")

if __name__ == "__main__":
    default_acc_time = 2

    csv_directory, analysis_directory, metadata_file, metadata_params, bin_params, min_trials_nos = get_user_input()

        
    metadata = pd.read_excel(metadata_file)
    df = loader.make_condition_df(csv_directory, analysis_directory, metadata, metadata_params[1], default_acc_time)

    freq_window = bin_params[1]
    freq_bin = bin_params[2]
    time_bin = bin_params[0]
    statistics, mean_statistics, counts, performance = calculate_statistics(df, freq_window, freq_bin, time_bin, min_trials_nos[0], 
                                                                            min_trials_nos[1], min_trials_nos[2])
    
    generate_trialcount_plot(counts)
    generate_performance_plots(performance)
    generate_licking_plots(statistics)

    # write data to file
    # for data, counts, performance: write Animal, Time (hr), Time (ms), lick/trial no(for counts)
    # name file with condition and type of data
    # write 3 files (data frames are not really alignable to put columns together meaningfully): 
    #     condition_lickfreq.xlsx, condition_trialcounts.xlsx, condition_performance.xlsx
