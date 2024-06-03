### 1. Folder and metadata format

The code expects behavior data presented in this format:

>     Condition1
>     |--Animal1
>        |--temptest.txt
>        |--temptest.txt
>     |--Animal2
>        |--temptest.txt
>        |--temptest.txt
>     |--Animal3
>        |--temptest.txt
>        |--temptest.txt
>     |--etc.

And requires a metadata excel sheet presemted in this format (columns in brackets are the default optional columns included, additional default columns can be inculded):

| Animal ID | acc | [Sex] | [Age] | [Date] | [Cage] | ... |
| --------- | --- | ----- | ----- | -----  | ------ | --- |
| Animal1   | 2   | F     | 25    | 2024-05-04 | Kourtney | |
| Animal2 | 1 | M | 30 | 2024-03-27 | Dream | |
| Animal3 | 2 | M | 27 | 2024-02-23 | Stormi | |
| etc. | ...  |   |    |            | | |

Columns can be in any order. `Animal ID` and `acc` are required. The names in the animal ID column must match the names of the folders containing behavior data. Additional columns can and animals can be present in the metadata file and will be ignored.


### `lickfrequency_pipepline.py`
    usage: lickfrequency_pipeline.py [-h] [-f FORMATTED_OUTPUT]
                                 [--output_name OUTPUT_NAME] [-o ANALYSIS_OUTPUT]
                                 [-a] [-m METADATA] [-d] [--time_bin TIME_BIN]
                                 [-W ANALYSIS_WINDOW ANALYSIS_WINDOW]
                                 [-I [INSTENTANEOUS_BIN]] [-M [MIN_TRIALS]]
                                 [-b [MIN_BLANK]] [-s [MIN_STIMULUS]]
                                 [-k COLUMNS [COLUMNS ...]] [-n NAME] [-v VALUES]
                                 animal [animal ...]

Formats and analyzes given data.

`-f`: If present, write unanalized formatted data (output of `format_arduino_data.py`) to given file.

`--output_name`: if present, name the output files with this `OUTPUT_NAME`. This can be different than name given by `-n`.

`-o`: Directory to write analysis to.

`-a`: If present, add to previous analysis. All parameters except input directory must be the same between the two analyses for correct function.

`-m`: path to metadata file. Required (see above).

`-d`: If present, treat directory as condition. Otherwise treat as single animal.

`--time_bin`: Duration of time bin in hours (default `4`).

`-W`: two values that give the start (inclusive) and end (exclusive) of the fixed analysis window (default `(700, 1000]`).

`-I`: size of instentaneous bin window (default `100`, minimum sample rate for arduino).

`-M`: Minimum number of trials to keep a bin (default `10`).

`-b`: Minimum number of blank trials to keep a bin (default `1`).

`-s`: Minimum number of stimulus trials to keep a bin (default `1`).

`-k`: Set of additional columns to keep from metadata file. Optional, but must be consistent across analyses if `-a` is set.

`-n`: Name of condition used in analysis. If `output_name` is not present, also used to name file.

`-v`: Can be `lick` to analyze lick frequency, `poke` to analyze poke frequency, or `lick, poke` to analyze both. Default `lick`.

`animal`: input directory for data. Can be a single animal or a condition directory containing multiple subdirectories for multiple animals.


### `lickfrequency_analysis.py`

    usage: lickfrequency_analysis.py [-h] [-t {full,bin,nth_part}] [-g ARGS [ARGS ...]]
                                 [-a] [-o OUT [OUT ...]]
                                 [-m DROP_MINIMUMS DROP_MINIMUMS DROP_MINIMUMS]
                                 [-k COLUMNS [COLUMNS ...]] [-v VALUES [VALUES ...]]
                                 [-i INDEX [INDEX ...]]
                                 files [files ...]

Analyzes given data, must have already been formatted.

`-t`: Analysis type to run. `full` analyzes all data, `bin` selects a given bin, and `nth_part` selects a given set of trials by day.

`-g`: Argument options for analysis. 

- `full` - freq_window, freq_bin, and time_bin. 
    - freq_window: length of rolling window (unused)
    - freq_bin: length of bin in ms to bin samples across trial
    - time_bin: length of bin in hours to bin trials across time
- `bin` - time_bin and analysis_window. 
    - time_bin: length of bin in hours to bin trials across time
    - analysis_window: start and end of fixed analysis window
- `nth_part` - num_parts and nth_part.
    - num_parts: number of parts to split trials into
    - nth_part: which part to take

`-a`: If present, add to previous analysis. All parameters except input directory must be the same between the two analyses for correct function.

`-o`: Directory to write analysis to.

`-m`: Number of total, stimulus, and blank trials to keep bin (default `10`, `1`, `1`).

`-k`: Set of additional columns to keep from metadata file. Optional, but must be consistent across analyses if `-a` is set.

`-v`: Can be `lick` to analyze lick frequency, `poke` to analyze poke frequency, or `lick, poke` to analyze both. Default `lick`.

`-i`: index to group data by during analysis.

`files`: formatted files to analyze.

### `format_arduino_data.py`

    usage: format_arduino_data.py [-h] [-c] [-o OUT] [-a] -m METADATA [METADATA ...]
                              [-n NAMES [NAMES ...]] [-k COLUMNS [COLUMNS ...]]
                              files [files ...]

Formats given data.

`-c`: If present, treat the provided directories as conditions

`-o`: File to write formatted data to.

`-a`: If present, add to previous analysis. All parameters except input directory must be the same between the two analyses for correct function.

`-m`: Path to metadata file(s). Required. If multiple directories are provided in `files`, must match the number of directories provided.

`-n`: Names of conditions. If multiple directories are provided in `files`, must match the number of directories provided.

`-k`: Set of additional columns to keep from metadata file. Optional, but must be consistent across analyses if `-a` is set.

`files`: Input animal or condition directories

### `plot_data.py`

Functions to help with plotting data, no cli function.