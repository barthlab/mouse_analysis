#!python3

import argparse
import csv
import datetime
import os

import matplotlib.pyplot as plt
import pandas as pd



def parse_args():
    default_acclimation_time = 48  # hours
    default_bin_time = 4  # hours
    default_last_x_percent = 20 # percent

    parser = argparse.ArgumentParser(description="Calculate mouse training statistics")
    parser.add_argument("--animal_name", required=True, type=str, help="Name of the animal, to name the output files")
    parser.add_argument("--csv_files", required=True, type=str, nargs="+", help="In order paths to the CSV files")
    parser.add_argument("--acclimation_time", type=float, default=default_acclimation_time, help=f"Acclimation time in hours (default: {default_acclimation_time})")
    parser.add_argument("--bin_time", type=float, default=default_bin_time, help=f"Bin time in hours (default: {default_bin_time})")
    parser.add_argument("--last_x_percent", type=float, default=default_last_x_percent, help=f"Percentage of data to consider (0 < last_x_percent < 100) (default: {default_last_x_percent})")

    parsed = parser.parse_args()

    # Check if all csv_files are valid files
    for csv_file in parsed.csv_files:
        if not os.path.isfile(csv_file):
            parser.error(f"The provided CSV file '{csv_file}' does not exist.")

    # Check if acclimation_time is greater than 0
    if parsed.acclimation_time < 0:
        parser.error("Acclimation time must be greater than 0.")

    # Check if bin_time is greater than 0
    if parsed.bin_time < 0:
        parser.error("Bin time must be greater than 0.")

    # Check if last_x_percent is in the range (0, 100)
    if not 0 < parsed.last_x_percent < 100:
        parser.error("last_x_percent must be in the range (0, 100).")

    return parsed



def make_df_from_csv_files(csv_files):
    df_list = []
    milli = 1e3

    # Read and concatenate the CSV files into a single DataFrame
    for csv_file in csv_files:
        # Parse the filename to extract the timestamp
        filename = os.path.basename(csv_file)

        # Extract the MM_DD_YY_T_HH_MM_SS part
        timestamp_str = filename.split(".txt")[1]
        timestamp = datetime.datetime.strptime(timestamp_str, "%m_%d_%y_T_%H_%M_%S").replace(tzinfo=datetime.timezone.utc)

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
    # If there is a 9 in the third row that means that the code was in the
    # no water part of the code, so we know that it was not a water trial
    return 9 not in trial[3].values



def calculate_statistics(trial):
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

    return start_time, lick_frq_before_water, is_water_trial



def bin_by_time(dataframe, ms_bin_time):
    # Cages are changed around noon, so that is when we start the bins
    return dataframe.groupby(pd.Grouper(key="start_time", freq=f"{ms_bin_time}ms", origin="1970-01-01 12:00:00"))



def calculate_final_statistics(bin):
    lick_frq_water = []
    lick_frq_blank = []
    trials_per_bin = len(bin)

    for _, row in bin.iterrows():
        start_time = row['start_time']
        lick_frq_before_water = row['lick_frq_before_water']
        is_water_trial = row['is_water_trial']

        if is_water_trial:
            lick_frq_water.append(lick_frq_before_water)
        else:
            lick_frq_blank.append(lick_frq_before_water)

    avg_lick_frq_water = sum(lick_frq_water) / len(lick_frq_water) if lick_frq_water else 0
    avg_lick_frq_blank = sum(lick_frq_blank) / len(lick_frq_blank) if lick_frq_blank else 0

    performance = avg_lick_frq_water - avg_lick_frq_blank

    return avg_lick_frq_water, avg_lick_frq_blank, performance, trials_per_bin



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



def make_graph(final_statistics, acclimation_time, bin_time, file_path):
    # Extract data columns
    time_bins = []
    avg_lick_frq_water = []
    avg_lick_frq_blank = []
    num_trials = []
    performance = []
    row_time = -1 * acclimation_time + bin_time / 2
    for row in final_statistics:
        time_bins.append(row_time)
        avg_lick_frq_water.append(row[0])
        avg_lick_frq_blank.append(row[1])
        performance.append(row[2])
        num_trials.append(row[3])
        row_time += bin_time

    # Create the figure and axes
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(8, 10))

    # First plot: Licking (hertz)
    axs[0].plot(time_bins, avg_lick_frq_water, color='green', label="Water Lick Freq (Hz)")
    axs[0].plot(time_bins, avg_lick_frq_blank, color='red', label="Blank Lick Freq (Hz)")
    axs[0].set_ylabel("Licking (Hz)")
    axs[0].legend()
    axs[0].set_xticks(time_bins)
    axs[0].set_xticklabels(time_bins)

    # Second plot: Number of Trials (per bin)
    axs[1].plot(time_bins, num_trials, color='blue', label="Number of Trials")
    axs[1].set_ylabel("Number of Trials")
    axs[1].legend()
    axs[1].set_xticks(time_bins)
    axs[1].set_xticklabels(time_bins)

    # Third plot: Performance
    axs[2].plot(time_bins, performance, color='purple', label="Performance")
    axs[2].set_ylabel("Performance")
    axs[2].set_xlabel("Training time (hours)")
    axs[2].legend()
    axs[2].set_xticks(time_bins)
    axs[2].set_xticklabels(time_bins)

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure as an image file
    plt.savefig(file_path)
    plt.close()



if __name__ == "__main__":

    args = parse_args()

    animal_name = args.animal_name
    csv_files = args.csv_files
    acclimation_time = args.acclimation_time
    bin_time = args.bin_time
    bin_time_ms = bin_time * 60 * 60 * 1000
    last_x_percent = args.last_x_percent

    df = make_df_from_csv_files(csv_files)

    start_indices = calculate_start_indices(df)

    trials = bin_into_trials(df, start_indices)

    trial_statistics = []
    for trial in trials:
        if not trial.empty:
            trial_statistics.append(calculate_statistics(trial))
    dataframe = pd.DataFrame(trial_statistics, columns=["start_time", "lick_frq_before_water", "is_water_trial"])
    dataframe["start_time"] = pd.to_datetime(dataframe["start_time"], unit="ms")


    # To use for bins
    binned_trial_statistics = bin_by_time(dataframe, bin_time_ms)

    final_bin_statistics = []
    for _, bin_data in binned_trial_statistics:
        statistics = calculate_final_statistics(bin_data)
        final_bin_statistics.append(statistics)

    save_data_to_csv(final_bin_statistics, acclimation_time, bin_time, f"{animal_name}_bin.csv")
    make_graph(final_bin_statistics, acclimation_time, bin_time, f"{animal_name}_bin.png")


    # To use for the last x%
    day_in_milliseconds = 24 * 60 * 60 * 1000
    day_binned_trial_statistics = bin_by_time(dataframe, day_in_milliseconds)

    final_day_statistics = []
    for _, day_data in day_binned_trial_statistics:
        start_index = int(len(day_data) * (1 - last_x_percent / 100))
        last_x_percent_of_day = day_data[start_index:]
        statistics = calculate_final_statistics(last_x_percent_of_day)
        final_day_statistics.append(statistics)

    save_data_to_csv(final_day_statistics, acclimation_time, bin_time, f"{animal_name}_day.csv")
    #make_graph(final_day_statistics, acclimation_time, bin_time, f"{animal_name}_day.png")
