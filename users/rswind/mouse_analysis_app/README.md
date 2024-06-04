# HOW TO USE `mouse_analysis_app.py`

# 1. INSTALLATION

## Option 1 - miniconda
1. Install Miniconda from this link: https://docs.anaconda.com/free/miniconda/

2. open a Miniconda prompt

3. run: `conda  create -n mouse_analysis python=3.11.3 numpy=1.26.4 pandas=2.2.1 matplotlib=3.8.4 tk=8.6.12 seaborn=0.13.0 openpyxl`


## Option 2 - python only

1. Intall python from this link: https://www.python.org/downloads/

2. open the command line

3. run `pip install numpy==1.26.4 pandas==2.2.1 matplotlib==3.8.4 seaborn==0.13.0 openpyxl`

# 2. RUNNING 

## Option 1 - Miniconda
1. From the Miniconda prompt, run: `conda activate mouse_analysis`

2. From the Miniconda prompt, run: `python3 [path_to_file]`. You should be able to copy and paste this from the file explorer into the command prompt.

    1. e.g. `python3 C:\Users\swind\Documents\lab_work\Barth\mouse_analysis_app.py`

## Option 2 - python only

1. Double-click on python executable file

or

2. Open and run in any IDE

# 3. USAGE

## 1. Folder and metadata format

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

Columns can be in any order. `Animal ID` and `acc` are required. The names in the animal ID column must match the names of the folders containing behavior data. Additional columns can and animals can be present in the metadata file and will be ignored. Animal names **must** not have underscores in them.

## 2. App layout
	
### 1. Analysis - File and parameter selection

**Time bin (hours):** Length of time bin to group trials

**Metadata columns to keep:** Addiditonal columns to look for in provided metadata file. These values are not modified from their format in the column (e.g. a Date column is left a a string representation). If adding to previous analysis, these should exactly match the columns the columns taken from previous analysis (see [Add to Previous output](#add)).

**Values to analyze:** Can be `lick` to analyze lick frequency, `poke` to analyze poke frequency, or `lick, poke` to analyze both. For most analyses, should be left as `lick`.

### 2. Analysis - Analysis Types

> **Important!** If adding to previous analysis, these selections should match previous analysis selections. Analysis results may be unpredictable if they do not. (see [Add to Previous output](#add))

**Instentaneous analysis:** Analysis of instentaneous lick frequency across the entire trial. At each 100ms interval data is collected throughout a trial, calculates instentaneous lick frequency for that timepoint, then averages lickfrequency across all trials in a timebin for that timepoint.

**Fixed window analysis:** Analysis of lick frequency across given fixed window. Left inclusive, right exclusive. Default `[700, 1000)` matches anticipatory lick frequency analysis.

- **start of fixed window:** Time in milliseconds from start of air puff delivery to start of fixed window. Inclusive.

- **end of fixed window:** Time in milliseconds from start of air puff delivery to end of fixed window. Exclusive.

**Fraction of day analysis:** Summary analysis of fixed window for a fraction of trials on each day. Default `5/5` corresponds to last 20% of trials.

- **Fraction of day:** Fraction representing which set of trials to take for this analysis. 
    - Denomenator represents how many bins to create (e.g. `5` gives bins with 20% of trials per day, and `4` gives bins with 25% of bins per day). 
    - Numerator represents which set of trials to take (e.g. `1/5` represents the first 20% of trials, and `3/4` represents the 3rd quarter of trials).

### 3. Minimum trial thresholds

These represent the minumum number of trials required to keep a bin for filtered analysis. Both filtered and unfiltered data is exported from analysis, but only filtered data is automatically plotted (see [Plotting](#plotting)).

**Minimum number of trials:** minimum number of total trials to keep a bin. Default `10`.

**Minimum number of blank trials:** minimum number of blank trials to keep a bin. Default `1`.

**Minimum number of stimulus trials:** minimum number of stimulus trials to keep a bin. Default `1`.

### 4. Output folder

Folder to output analysis results to.

#### Add to previous analysis {# add}

If selected, this will append the current analysis to a set of previously run analyses. The set of analyses run, the parameters for analysis, and the metadata columns kept between the old and new analysis should exactly match for correct behavior. If they do not, the analysis may seem to run correctly but write data incorrectly, or may fail entirely, or may work for some analyses but not others.

#### 5. File type

This indicates whether the selected folder represents a single animal, or a condition, which contains multiple animals with data to analyze.

**Animal/condition name:** (optional) If provided, this will be the name added to the output files; however, the actual condition or animal name in the analysis will be taken from the folder name.

## 3. Plotting

### 1. Data to plot

Folder to output of previously run analysis.

### 2. Fixed window plots

**Fixed window lick frequency:** If selected, plots stimulus (green) and blank (red) lick frequency by condition.

**Fixed window performance:** If selected, plots performance (lick to stimulus - lick to blank)  by condition.

### 3. Instentaneous plots

**Instentaneous lick frequency by hour:** If selected, plots instenaneous lick frequency across trials, with each timebin as a separate subplot. Makes a separate graph for each condition.

**Instentaneous performance by hour:** If selected, plots instenaneous performance across trials, with each timebin as a separate subplot. Makes a separate graph for each condition.

**Instentaneous lick frequency heatmap:** If selected, plots a heatmap of instentaneous lick frequency, with y-axis as time bin and x-axis as time in trial from air puff start. Blue is lower lick frequency and yellow is higher.

**Instentaneous performance heatmap:** If selected, plots a heatmap of instentaneous performance, with y-axis as time bin and x-axis as time in trial from air puff start. Blue is negative performance, white is 0, and red is positive performance.

### 4. Trial plots

**Nth fraction of day stimulus vs blank:** If selected, plots a stimulus vs blank pair plot by animal. Makes a separate subplot for each day and condition, with conditions as columns and days as rows.

**Trials by hour:** If selected, plots the number of trials per timebin. Makes a separate subplot for each animal.

### 5. Save options

**Save plots:** If selected, will save the plots to the selected folder.

**Display plots:** If selected, will display the plots in the widget.