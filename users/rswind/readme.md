Author: Rachel Swindell  
Last Modified: 09-05-23  

***

### **file list**
 
trace_analysis.py  
all_conditions.ps1  
analysis_functions.py  
loader.py  
v16_analysis.inpb  
Eunsol_plots_notebook.inpb  
Rachel_plots_notebook.inpb  

***


### **trace_analysis.py**

Calculates average lick frequency, performance, and number of trials across entire duration of trial. Data is still binned over the duration of the training.
Designed to be run with an entire set of condition data that includes multiple animals, but can also be run on one animal. 

Can be run from the command line or with simple GUI prompts.

**Command Line Options**

```
usage: trace_analysis.py [-h] [--noUI] [-r] [-n CONDITION_NAME] [-t TRIAL_BIN] [-m METADATA] [-a ACCLIMATION_TIME]  
                         [-b BIN_TIME] [-p LAST_X_PERCENT] [-i MIN_TRIALS] [-w MIN_WATER_TRIALS] [-k MIN_BLANK_TRIALS]  
                         [condition_path] [output_path]  

positional arguments:  
  condition_path        Path to directory with data files  
  output_path           Path to directory to output analysis to  

options:  
-h, --help            show this help message and exit  
  --noUI                Run from command line, if present  
  -r, --rolling         Use rolling window instead of licking frequency bins, if present  
  -n CONDITION_NAME, --condition_name CONDITION_NAME  
                        Name of condition type (default: SAT1)  
  -t TRIAL_BIN, --trial_bin TRIAL_BIN  
                        Licking frequency bin size in ms (default: 300)  
  -m METADATA, --metadata METADATA  
                        Path to file with animal metadata  
  -a ACCLIMATION_TIME, --acclimation_time ACCLIMATION_TIME  
                        Name of column in metadata file with acclimation times (default: ACC days)  
  -b BIN_TIME, --bin_time BIN_TIME  
                        Bin time in hours (default: 4)  
  -p LAST_X_PERCENT, --last_x_percent LAST_X_PERCENT  
                        Percentage of data to consider by day (0 < last_x_percent < 100) (default: 20)  
  -i MIN_TRIALS, --min_trials MIN_TRIALS  
                        Minimum number of trials in a bin (default: 0)  
  -w MIN_WATER_TRIALS, --min_water_trials MIN_WATER_TRIALS  
                        Minimum number of water trials in a bin (default: 1)  
  -k MIN_BLANK_TRIALS, --min_blank_trials MIN_BLANK_TRIALS  
                        Minimum number of blank trials in a bin (default: 1)
```

**Run with GUI**

To run with GUI, don't include any arguments in the command line.

usage: python3 trace_analysis.py

GUI prompt structure:

1. Folder where data is located
	- Expected file structure (code will not work if not structured this way):

  >  ```
  >  condition_folder
  >    |-- animal_folder
  >        |-- temptest.txt
  >    |-- animal_folder
  >        |-- temptest.txt
  >    |-- etc.
  >  ```

1. Folder to output analysis

1. Excel file containing relevant metadata about animals

1. Popup with 2 prompts
	1. Condition Name - name for the condition the animals were trained under (default SAT)
	
	1. Acclimation Name - name of column in metadata file containing acclimation time of animals (default ACC days)
	
1. Popup with 4 prompts
	1. Bin Time - Length of timebins in minutes (default 240 - 4 hours)
	
	1. Frequency Rolling Average window size (ms) (default 100 ms)
	
	1. Frequency Bin size (ms) (default 100 ms)
		- This code can bin licking frequency over the trial into a rolling window or discrete bins. Only one of these two options should be set, and the other should be left at the sample rate (the default 100 ms value)
		
	1. Last percent of day (default 20)
		- The percent of trials to output vaules for per day. Used to correlate behavior with other measures.
		
1. Popup with 4 prompts
	1. Minimum total trials (default 10)
		- the minimum total number of trials in a bin to keep it
	
	1. Minimum blank trials (default 1)
		- the minimum number of blank trials in a bin to keep it
		
	1. Minimum water trials (default 1)
		- the minimum number of water trials in a bin to keep it
	
	- Any bin that has fewer than the specified number of trials will not be included in the output. Set all three values to 0 to keep all bins.
	
All popups will not resolve until valid values are given.
	
### **analysis_functions.py**

Library of helper functions for data analysis.

### **loader.py**

Library of helper functions for loading data into pandas dataframe from `.txt` files and attaching metadata to data.

### **all_conditions.ps1**

This is a powershell script I wrote for my own convenience. It runs `trace_analysis.py` on multiple conditions with the same parameters, and writes the output to a specified folder.
It runs from the command line assuming you are in  the same directory as `trace_analysis.py`. As written, it runs using a rolling window 300 ms average, across 0.5, 1, 2, and 4 hours bins.
All other parameters are default. All paths would need to be modified to be of use on a different system, and parameters can be modified as needed or as of interest.

### **v16_analysis.inpb**

This file is a notebook comparing that the analysis approach I took using Pandas and the old analysis code (SomatosensoryLearning_LickAnalysisv16.m)

### **Eunsol_plots_notebook.inpb and Rachel_plots_notebook.inpb**

These are Jupyter notebooks I used to visualize the analyzed output from `trace_analysis.py`. The file path would need to be modified to be useful on other systems. 
Eunsol_plots_notebook contains plots for the analysis of Eunsol's data, and Rachel_plots_notebook contains pltos for the analysis of my data.

### **video_alignment.inpb**

This is a Jupyter notebook I used to align video recordings of trials and trial lick frequency data. Requires ffmpeg and ffmpeg-python.