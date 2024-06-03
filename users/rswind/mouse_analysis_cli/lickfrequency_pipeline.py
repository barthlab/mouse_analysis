import os
import sys
import argparse
import pandas as pd
import format_arduino_data
import lickfrequency_analysis
import plot_data



if __name__ == '__main__':
    # TODO: file dialog to select files/conditions/animals/metaconditions
    parser = argparse.ArgumentParser()
    parser.add_argument('animal', nargs='+', help='Animal directory(s) to analyze')
    parser.add_argument('-f', '--formatted_output', help='File to write or add formatted data to')
    parser.add_argument('--output_name',help='Name to append to output files.')
    parser.add_argument('-o', '--analysis_output', help='Directory to write analysis output to')
    parser.add_argument('-a', '--add',action='store_true', help='If provided, append additional data to output file. Default behavior is to overwrite data in output file.')
    parser.add_argument('-m', '--metadata', help='Metadata file to use for analysis')
    parser.add_argument('-d', action='store_true', help='If provided, treat given path as a condition directory containing multiple animals')
    parser.add_argument('--time_bin', default=4*60, nargs=1, help='Time bin size for analysis in minutes')
    parser.add_argument('-W', '--analysis_window', nargs=2, default=[700, 1000], help='Analysis window for anticipatory analysis, in ms from air puff initiation', type=int)
    parser.add_argument('-I', '--instentaneous_bin', nargs='?', default=100, help='Bin size for instentaneous analysis, in ms', type=int)
    parser.add_argument('-M', '--min_trials', nargs='?', default=10, help='Miniumum number of trials to keep a bin', type=int)
    parser.add_argument('-b', '--min_blank', nargs='?', default=1, help='Miniumum number of blank trials to keep a bin', type=int)
    parser.add_argument('-s', '--min_stimulus', nargs='?', default=1, help='Miniumum number of stimulus trials to keep a bin', type=int)
    parser.add_argument('-k', '--columns', nargs='+', help='Set of columns to keep for analysis')
    parser.add_argument('-n', '--name', help='Name of condition')
    parser.add_argument('-v', '--values', help="Name of values to analyze", default=['lick'])
    args = parser.parse_args()

    dr = args.animal[0].split('\\')
    if len(dr) == 1:
        dr = args.animal[0].split('/')
    name = args.name
    if not args.name:
        name = dr[-2]
    parent = '\\'.join(dr[:-1])

    # format data
    formatted_res = format_arduino_data.run_formatter(args.animal, args.d, args.formatted_output, args.add, [args.metadata], [args.name], args.columns)
    if formatted_res.empty:
        print('Stopping analysis - no animals loaded.')
        sys.exit(1)
    # run analysis - instentaneous lick frequency, anticipatory lick frequency, last 20% of trials for anticipatory lick frequency
    mins = {
        "min_trials": args.min_trials,  # minimum number of trials in a bin to keep
        "min_blank": args.min_blank,   # minimum number of water trials
        "min_stimulus": args.min_stimulus # minimum number of water trials 
        }

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
    keep = args.columns + ['acc', 'stimulus', 'water', 'type']    

    full_suffix = ['data', 'means', 'means_filtered', 
                   'performance', 'performance_filtered', 
                   ]
    ant_suffix = full_suffix + ['trial_counts', 'trials_summed_day'] 

    last20_suffix = ['data', 'means', 'performance', 'trial_counts']
    
    if args.analysis_output:
        if args.add:
            analysis_output = os.listdir(args.analysis_output)
            analysis_output.sort()
            full_out = [args.analysis_output + '\\' +  analysis_output[i] for i in range(7,12)]
            ant_out = [args.analysis_output + '\\' +  analysis_output[i] for i in range(0,7)]
            last20_out = [args.analysis_output + '\\' +  analysis_output[i] for i in range(12,16)]
        else:
            output_name = name
            if args.output_name:
                output_name = args.output_name
            full_out = [args.analysis_output + '\\full_'  + output_name + '_' + suff + '.txt' for suff in full_suffix]
            ant_out = [args.analysis_output + '\\anticipatory_'  + output_name + '_' + suff + '.txt' for suff in ant_suffix]
            last20_out = [args.analysis_output + '\\last20_'  + output_name + '_' + suff  + '.txt' for suff in last20_suffix]

    # ignore columns provided in 'keep' that are not present in data
    for c in keep:
        if c not in formatted_res.columns:
            keep.remove(c)
        
    print('Full lick frequency analysis for %s' % name)
    full_analysis = lickfrequency_analysis.run_analysis(formatted_res, 'full', 
                        full_keys, mins, keep, args.values, index, 
                        [None, args.instentaneous_bin, args.time_bin], full_out, args.add)

    print('Anticipatory lick frequency analysis for %s' % name)
    ant_analysis = lickfrequency_analysis.run_analysis(formatted_res, 'bin', 
                        keys, mins, keep, args.values, index, [args.analysis_window, args.time_bin], ant_out, args.add)
    print(f'Last 20% of trials analysis for {name}')
    last20 = lickfrequency_analysis.run_analysis(ant_analysis[0], 'nth_part', 
                        keys, mins, keep, args.values, index, [5, 5], last20_out, args.add)
    
    # # for reference: all super-conds I have
    # conds = [    
    #     'C:\\Users\\swind\\Documents\\lab_work\\Barth\\behavior_analysis\\behavior_data\\raw_data\\Eunsol\\all_animals\\dcz',
    #     'C:\\Users\\swind\\Documents\\lab_work\\Barth\\behavior_analysis\\behavior_data\\raw_data\\Eunsol\\all_animals',
    #     # r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\Rachel\good',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\Eunsol\new vs old acc\new acc',
    #    r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\retrain',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\lesion\Surgery done by Kate',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\joe\Joes C57 training data',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\joe\L4 scnn-cre',
    #  ]
    # # and corresponding metadata files
    # metas = [
    #     'C:\\Users\\swind\\Documents\\lab_work\\Barth\\behavior_analysis\\behavior_data\\raw_data\\Eunsol\\eunsol_metadata.xlsx',
    #     'C:\\Users\\swind\\Documents\\lab_work\\Barth\\behavior_analysis\\behavior_data\\raw_data\\Eunsol\\eunsol_metadata.xlsx',
    #     # r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\Rachel\animal summary.xlsx',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\Eunsol\new vs old metadata.xlsx',
    #    r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\retrain_metadata.xlsx',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\lesion\lesion metadata.xlsx',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\joe_metadata.xlsx',
    #     r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\L4_metadata.xlsx'
    # ]
    # cond1 = r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\raw_data\joe\joe 6ohda'
    # meta1 = r'C:\Users\swind\Documents\lab_work\Barth\behavior_analysis\behavior_data\joe_metadata.xlsx'

