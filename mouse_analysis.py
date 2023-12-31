#!python3
version = "v1.1"
# TODO why are we printing 0 for performance and the other values when we don't have that many results? for some partitions we are printing 0 trials, even though we should have trials for those partitions, why is that?

import argparse
import csv
import datetime
import os
import sys
import shutil

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import pandas as pd
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox



def validate_pos_float(input_str):
    try:
        float_val = float(input_str)
        if float_val >= 0:
            return float_val
        else:
            return None
    except ValueError:
        return None


def validate_last_x_percent(input_str):
    float_val = validate_pos_float(input_str)
    if float_val is None or float_val <= 0 or float_val > 100:
        return None
    return float_val


def select_directory(prompt):
    directory = filedialog.askdirectory(title=prompt)
    if directory:
        return directory
    else:
        sys.exit("No directory selected")


def get_files_in_dir(directory):
    files = []
    for file in os.listdir(directory):
        file_path = os.path.join(directory, file)
        if os.path.isfile(file_path):
            files.append(file_path)
    return files


def get_user_input():
    csv_directory = select_directory("Select directory with CSV files")

    while True:
        acclimation_time = validate_pos_float(input("Enter acclimation time (in hours): "))
        if acclimation_time is not None:
            break
        else:
            print("Invalid input. Please enter a valid floating-point number.")

    while True:
        bin_time = validate_pos_float(input("Enter bin time (in hours): "))
        if bin_time is not None:
            break
        else:
            print("Invalid input. Please enter a valid floating-point number.")

    last_x_percent = 20 

    # get last x% value from user - not used by anyone in lab (delete?)
    # while True:
    #    last_x_percent = validate_last_x_percent(input("Enter the last x percent (0 < x <= 100): "))
    #    if last_x_percent is not None:
    #        break
    #    else:
    #        print("Invalid input. Please enter a valid floating-point number between 0 and 100.")

    analysis_directory = csv_directory + '/' + f"analysis_{acclimation_time}_{bin_time}{version}"
    if os.path.exists(analysis_directory):
        shutil.rmtree(analysis_directory)  
    os.makedirs(analysis_directory)

    return get_files_in_dir(csv_directory), analysis_directory, acclimation_time, bin_time, last_x_percent



def make_df_from_csv_files(csv_files):
    df_list = []
    milli = 1e3

    # Read and concatenate the CSV files into a single DataFrame
    for csv_file in csv_files:
        # Parse the filename to extract the timestamp
        filename = os.path.basename(csv_file)

        # Extract the MM_DD_YY_T_HH_MM_SS part
        timestamp_str = filename.split(".txt")[1]
        try:
            timestamp = datetime.datetime.strptime(timestamp_str, "%m_%d_%y_T_%H_%M_%S").replace(tzinfo=datetime.timezone.utc)
        except:
            try:
                timestamp = datetime.datetime.strptime(timestamp_str, "%m_%d_%y~T~%H_%M_%S").replace(tzinfo=datetime.timezone.utc)
            except:
                sys.exit("File timestamp not recognized.")


        # Convert to milliseconds because that is the unit we are using for everything in the backend
        epoch_time = timestamp.timestamp() * milli

        # Read the CSV file and update the time in the first column
        df = pd.read_csv(csv_file, header=None)
        df[0] *= milli
        df[0] += epoch_time
        df_list.append(df)

    df_combined = pd.concat(df_list)
    df_combined.reset_index(drop=True, inplace=True)

    return df_combined



def calculate_start_indices(df):
    # Find the indices where a 7 appears in the 4th column, this indicates it is the
    # end of a trial
    end_indices = df[df[3] == 7].index

    # The indices after an end indice are candidate_start_indices
    candidate_start_indices = end_indices + 1

    # Calculate the index of the candidate_start_indices that are not also end indices
    candidate_start_indices_not_in_end_indices = ~candidate_start_indices.isin(end_indices)

    # Get the start indices from the candidate_start_indeces
    start_indices = candidate_start_indices[candidate_start_indices_not_in_end_indices]

    return start_indices



def bin_into_trials(df, start_indices):
    trials = []
    prev_start_index = 0

    for start_index in start_indices:
        # Extract the rows between the previous start index and the current start index
        trial = df.iloc[prev_start_index:start_index].reset_index(drop=True)
        trials.append(trial)
        prev_start_index = start_index

    # Add the remaining rows after the last start index as the last trial
    last_trial_data = df.iloc[prev_start_index:].reset_index(drop=True)
    trials.append(last_trial_data)

    return trials



def calculate_start_time(trial):
    # The trial starts on the timestamp of the first row of data collected in the trial
    return trial[0][0]



def calculate_variable_delay(trial):
    # Get the largest number in the 4th column
    # because it is 0 during most of the trial, we need to ignore those
    # and then the rest of the time the number is the variable delay
    # so we can get any of those
    variable_delay_ms = trial[4].max()

    return variable_delay_ms



def calculate_lick_frq_before_water(trial, start_lick_time, end_lick_time):
    # Get the number of trials that were in between start_lick_time and end_lick_time
    # AND that had a 2 in the 2nd column (indicating a lick)
    filtered_trials = trial[(start_lick_time <= trial[0]) & (trial[0] < end_lick_time) & (trial[2] == 2)]

    return len(filtered_trials)



def calculate_is_water_trial(trial):
    # If there is a 9 in the third column that means that the code was in the
    # no water part of the code, so we know that it was not a water trial
    return 9 not in trial[3].values


def check_is_decoupled(df):
    # If there are 6 columns in the data, the file differentiates between water and stimulus trials,
    # there is a meaningful distinction that must be checked for
    return len(df.columns) == 6

def calculate_is_stimulus_trial(trial):
    # If the sixth column exists and doesn't have a 0 that means that the code was in the
    # stimulus trial part of the code, so we know that it was a stimulus trial
    return (0 not in trial[5].values)

def calculate_statistics(trial, is_decoupled):
    milli = 1e3

    # Time before the water comes out where we want to measure lick frq
    anticipatory_licking_window = 300 # milliseconds
    puff_time = 500 # milliseconds
    fixed_delay = 500 # milliseconds

    start_time = calculate_start_time(trial)

    variable_delay = calculate_variable_delay(trial)

    end_lick_time = start_time + variable_delay + puff_time + fixed_delay
    start_lick_time = end_lick_time - anticipatory_licking_window

    lick_frq_before_water = calculate_lick_frq_before_water(trial, start_lick_time, end_lick_time) / (anticipatory_licking_window / milli)

    is_water_trial = calculate_is_water_trial(trial)

    # If not decoupled, a stimulus trial is the same as a water trial
    is_stimulus_trial = is_water_trial

    if (is_decoupled):
        is_stimulus_trial = calculate_is_stimulus_trial(trial)

    return start_time, lick_frq_before_water, is_water_trial, is_stimulus_trial



def bin_by_time(dataframe, ms_bin_time):
    # Cages are changed around noon, so that is when we start the bins
    return dataframe.groupby(pd.Grouper(key="start_time", freq=f"{ms_bin_time}ms", origin="1970-01-01 12:00:00"))



def calculate_final_statistics(bin, is_decoupled):
    # If not decoupled, is_stimulus_trial == is_water_trials, so lck_freq_water == lck_freq_stim
    # If decoupled, they are different, and the relevant comparision is stimulus vs. no stimulus
    lick_frq_stim = []
    lick_frq_blank = []
    trials_per_bin = len(bin)
    num_stimwater = 0
    num_stimdry = 0
    num_blankwater = 0
    num_blankdry = 0

    for _, row in bin.iterrows():
        start_time = row['start_time']
        lick_frq_before_water = row['lick_frq_before_water']
        is_water_trial = row['is_water_trial']
        is_stimulus_trial = row['is_stimulus_trial']

        if is_stimulus_trial:
            lick_frq_stim.append(lick_frq_before_water)
            if is_water_trial:
                num_stimwater += 1
            else:
                num_stimdry += 1
        else:
            lick_frq_blank.append(lick_frq_before_water)
            if is_water_trial:
                num_blankwater += 1
            else:
                num_blankdry += 1
        

    avg_lick_frq_stim = sum(lick_frq_stim) / len(lick_frq_stim) if lick_frq_stim else 0
    avg_lick_frq_blank = sum(lick_frq_blank) / len(lick_frq_blank) if lick_frq_blank else 0

    performance = avg_lick_frq_stim - avg_lick_frq_blank

    if is_decoupled: 
        result = avg_lick_frq_stim, avg_lick_frq_blank, performance, num_stimwater, num_stimdry, num_blankwater, num_blankdry, trials_per_bin
    else:
        result = avg_lick_frq_stim, avg_lick_frq_blank, performance, num_stimwater, num_blankdry, trials_per_bin
    return result



def save_data_to_csv(final_statistics, acclimation_time, bin_time, is_decoupled, file_path):
    # Column labels
    labels = ["Time Bin", "Avg Water Lick Freq (Hz)", "Avg Blank Lick Freq (Hz)", "Performance", 
              "Stimulus + Water trials", "Stimulus + No Water trials", "Blank + Water trials", "Blank + No Water trials","Trial Number"]
    # Drop irrelevant labels for non-decoupled analysis
    if not is_decoupled:
        del labels[5:7]

    # Open the CSV file for writing
    with open(file_path, 'w+', newline='') as csv_file:
        writer = csv.writer(csv_file)

        # Write column labels
        writer.writerow(labels)

        # Write data rows
        row_time = -1 * acclimation_time + bin_time / 2
        for row in final_statistics:
            data_row = list(row)
            data_row.insert(0, row_time)
            writer.writerow(data_row)
            row_time += bin_time



def make_graph(final_statistics, acclimation_time, bin_time, file_path):
    # Extract data columns
    time_bins = []
    avg_lick_frq_stim = []
    avg_lick_frq_blank = []
    num_trials = []
    performance = []
    row_time = -1 * acclimation_time + bin_time / 2
    for row in final_statistics:
        time_bins.append(row_time)
        avg_lick_frq_stim.append(row[0])
        avg_lick_frq_blank.append(row[1])
        performance.append(row[2])
        num_trials.append(row[3])
        row_time += bin_time

    # Create the figure and axes
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(8, 10))

    # First plot: Licking (hertz)
    axs[0].plot(time_bins, avg_lick_frq_stim, color='green', label="Stimulus Lick Freq (Hz)", marker='o')
    axs[0].plot(time_bins, avg_lick_frq_blank, color='red', label="Blank Lick Freq (Hz)", marker='o')
    axs[0].set_ylabel("Licking (Hz)")
    axs[0].set_ylim([0,12])
    
    # Second plot: Number of Trials (per bin)
    axs[1].plot(time_bins, num_trials, color='black', label="Number of Trials")
    axs[1].set_ylabel("Number of Trials")
    axs[1].set_ylim([0, 700])
    axs[1].fill_between(time_bins, 0, num_trials, facecolor='grey', alpha=0.5)

    # Third plot: Performance
    axs[2].plot(time_bins, performance, color='black', label="Performance", marker='o')
    axs[2].set_ylabel("Performance")
    axs[2].set_xlabel("Training time (hours)")
    axs[2].set_ylim([-10, 10])
    axs[2].axhline(y=0, xmin=0,xmax=1, ls='--', lw=0.75, color='black',zorder=0)
    
    for ax in axs:
        ax.legend()
        ax.tick_params(labelbottom=True)
        ax.xaxis.set_major_locator(MultipleLocator(12))
        ax.xaxis.set_minor_locator(MultipleLocator(6))
        ax.xaxis.set_major_formatter('{x:.0f}')
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure as an image file
    plt.savefig(file_path)
    plt.close()



if __name__ == "__main__":

    csv_files, analysis_dir, acclimation_time_hours, bin_time_hours, last_x_percent = get_user_input()

    bin_time_ms = bin_time_hours * 60 * 60 * 1000

    df = make_df_from_csv_files(csv_files)
    
    is_decoupled = check_is_decoupled(df)

    start_indices = calculate_start_indices(df)

    trials = bin_into_trials(df, start_indices)

    trial_statistics = []
    for trial in trials:
        if not trial.empty:
            trial_statistics.append(calculate_statistics(trial, is_decoupled))
    dataframe = pd.DataFrame(trial_statistics, columns=["start_time", "lick_frq_before_water", "is_water_trial", "is_stimulus_trial"])
    dataframe["start_time"] = pd.to_datetime(dataframe["start_time"], unit="ms")


    # To use for bins
    binned_trial_statistics = bin_by_time(dataframe, bin_time_ms)

    final_bin_statistics = []
    for _, bin_data in binned_trial_statistics:
        statistics = calculate_final_statistics(bin_data, is_decoupled)
        final_bin_statistics.append(statistics)

    file_prefix = os.path.basename(csv_files[0]).split(".txt")[0]
    save_data_to_csv(final_bin_statistics, acclimation_time_hours, bin_time_hours, is_decoupled, f"{analysis_dir}/{file_prefix}_{bin_time_hours}_bins{version}.csv")
    make_graph(final_bin_statistics, acclimation_time_hours, bin_time_hours, f"{analysis_dir}/{file_prefix}_{bin_time_hours}_bins_plot{version}.png")


    # To use for the last x%
    day_in_milliseconds = 24 * 60 * 60 * 1000
    day_in_hours = 24
    day_binned_trial_statistics = bin_by_time(dataframe, day_in_milliseconds)

    final_day_statistics = []
    for _, day_data in day_binned_trial_statistics:
        start_index = int(len(day_data) * (1 - last_x_percent / 100))
        last_x_percent_of_day = day_data[start_index:]
        statistics = calculate_final_statistics(last_x_percent_of_day, is_decoupled)
        final_day_statistics.append(statistics)

    save_data_to_csv(final_day_statistics, acclimation_time_hours, day_in_hours, is_decoupled, f"{analysis_dir}/{file_prefix}_last20percent{version}.csv")
    #make_graph(final_day_statistics, acclimation_time_hours, bin_time_hours, f"{analysis_dir}/day.png")
