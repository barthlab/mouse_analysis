% Somatosensory learning: Lick Analysis Code Version 16 outputs
% individual Trials from last 20% as well as plotting data.
% Just save the multiple split files and the code will automatically merge them and proceed with analysis.
% This code outputs all trials, and analyzes the anticipatory lick frequency for water versus no-water trials. 
% The output analysis includes lick frequency binned at every four hours and the last 20% each day.
% In addition, the code also creates a last 20% excel file containing the details of individual trials of the last 20%
% window, and combines the instantaneous lick frequency, and calculates lick probability of either 0, 1, 2, or 3 licks in
% that 0.3s anticipatory window. 
% Initially created by: Sarah Bernhard
% Modified by: Alexander Hsu, Finalized Date: 01/31/2019
% Modified by: Ajit Ray, Finalized: 10/22/2019
%% 

clear all;
clc;

%% Locate the folder containing the file, import .txt files
%%% Please note that this modified code only works if it meets the following conditions:
%%% 1. You have save the .txt file from chrome using right click and save
%%% link as.
%%% 2. There is only one .txt file in that particular input folder
fprintf('Please select animal folder to import the .txt files:');
filepathIn = uigetdir;
%file name format: (user added txt)temptest.txt12_29_
txtfilebasename = dir(sprintf('%s%s%s%s', filepathIn, '/', '*','txt'));
%% Following conditional statement allows you to merge split files automatically
if length(txtfilebasename) > 1
    temptest= [];
    for i = 1:length(txtfilebasename)
        temptest_ind{i} = table2array(readtable(sprintf('%s%s', filepathIn, '/',txtfilebasename(i).name)));
        textfilenamestr_ind{i} = cellstr(txtfilebasename(i).name);
        Key   = '_';
        Index{i} = strfind(textfilenamestr_ind{i}, Key);
        %get start time from file name
        timeStarted_ind{i} = datetime(str2num(strcat('20',textfilenamestr_ind{i}{1}(Index{i}{1}(3)-2:(Index{i}{1}(3))-1))), ...
            str2num(textfilenamestr_ind{i}{1}(Index{i}{1}(1)-2:Index{i}{1}(1)-1)), ...
            str2num(textfilenamestr_ind{i}{1}(Index{i}{1}(2)-2:Index{i}{1}(2)-1)) ...
            ,str2num(textfilenamestr_ind{i}{1}(Index{i}{1}(5)-2:Index{i}{1}(5)-1)), ...
            str2num(textfilenamestr_ind{i}{1}(Index{i}{1}(6)-2:Index{i}{1}(6)-1)), ...
            str2num(textfilenamestr_ind{i}{1}(Index{i}{1}(6)+1:Index{i}{1}(6)+2)));
        if i == 1
            temptest=[temptest;temptest_ind{i}];
        elseif i >= 2
            timediff{i-1}=seconds(timeStarted_ind{i}-timeStarted_ind{1});
            temptest_ind{i}(:,1) = temptest_ind{i}(:,1)+ timediff{i-1};
            temptest=[temptest;temptest_ind{i}];
        end
    end
    timeStarted = timeStarted_ind{1};
    textfilenamestr = textfilenamestr_ind{1}{1};
    T=array2table(temptest);
    writetable(T,[textfilenamestr(1:end-27),'-merged',textfilenamestr(end-26:end)]);
else
    temptest = table2array(readtable(sprintf('%s%s', filepathIn, '/',txtfilebasename(1).name)));
    textfilenamestr = txtfilebasename(1).name;
    Key   = '_';
    Index = strfind(textfilenamestr, Key);
    %get start time from file name
    timeStarted = datetime(str2num(strcat('20',textfilenamestr(Index(3)-2:Index(3)-1))), ...
        str2num(textfilenamestr(Index(1)-2:Index(1)-1)),str2num(textfilenamestr(Index(2)-2:Index(2)-1)) ...
        ,str2num(textfilenamestr(Index(5)-2:Index(5)-1)),str2num(textfilenamestr(Index(6)-2:Index(6)-1)), ...
        str2num(textfilenamestr(Index(6)+1:Index(6)+2)));
end

%% Select the output animal folder to store lick analysis excel files
%%% Select the folder where you want all five excel files and one figure file
%%% to be stored. Five files include : (1)FILE WITH ORIGINAL NAME - Excel file contents include 
%%% day-wise and type-wise raw data split across sheets; (2)FILE WITH 
%%% ORIGINAL NAME_AnalysisFinal- Excel file contents include colored rows showing
%%% day-wise and time bin-wise split licking and trial data in one sheet; (3) FILE WITH
%%% ORIGINAL NAME_DataFinal - Excel file contents include day-wise and
%%% time-bin wise split licking and trial data along with respective raw data
%%% spread across sheets; (4) FILE WITH ORIGINAL NAME_ForPlotting - Excel
%%% file containing arranged data for directly pasting into original graph
%%% template files for plotting; (5) FILE WIHT ORIGINAL NAME_Last20Data -
%%% Excel file containing sheets carrying information about raw and
%%% processed lick data for the last 20% of trials split day-wise; (6) FILE
%%% WITH ORIGINAL NAME_Plots - Matlab figure file showing the plots of lick
%%% frequence, trial and performance as three separate plots in one figure.

fprintf('\n');
fprintf('Please select folder to store analysis files:');
fprintf('\n');
filepathOutResults = uigetdir;

%%% To identify the number of days of acclimation for plotting
accprompt = 'How many days did you acclimate the animal? Please input digit and press the ENTER key. : ';
accdays = input(accprompt);
fprintf('\n');

%%% To check whether lick port is functional and/or animal has engaged in
%%% any licks.
%there are samples where lick column is non zero
lickporttest = any(temptest(:,3));
if lickporttest == 0
    fprintf('Cannot analyze the file. No licks recorded. Either animal did not lick or lick-port is non-functional.');
    fprintf('\n');
else
    %records what time data started at
    hourRecorded = hour(timeStarted);
    minuteRecorded = minute(timeStarted);
    secRecorded = second(timeStarted);

    %time started as a fraction of an hour
    hourStarted = hourRecorded + minuteRecorded/60 + secRecorded/3600;

    %puts hour into table
    temptest(:,6) = hourStarted(:,1);
    %converts millis into days since start - first column records time since file initiation in ms; 
    % this converts it into a fraction of a day
    daysSinceStart = temptest(:,1)/86400;

    %puts days since start into table
    temptest(:,7) = daysSinceStart(:,1);
    
    %converst days since start into hours since start - same as above but with hours instead of days
    hoursSinceStart = temptest(:,7)*24;
    %puts hours since start into table
    temptest(:,8) = hoursSinceStart(:,1);
    randDelay = temptest(:,5);
    %records the time(in hours) - time the sample was recorded as a fraction of an hour (e.g. 12.568 = 12:36:55)
    temptest(:,9) = hourStarted + hoursSinceStart;
    temptest(:,10) = randDelay(:,1);
    filledMillis = temptest(:,10);
    poke = temptest(:,2);
    %Setting up to fix the poke
    temptest(:,11) = poke(:,1);
    initialPoke = temptest(:,11);
    %fixes the poke, if has a 1 and 1 after then make second 1 into 0 - column 2 records if there was an ir break at that sample
    %so an extended ir break = multiple "pokes" recorded even though there was only one poke
    for i = 1:length(poke)-1
        if ((poke(i+1,1) == 1) && (poke(i,1) == 1))
            initialPoke(i+1,1) = 0;
        end
    end
    %fixing the delay so there are no 0s in column, if 0 before then use number
    %in previous column
    for i = 1:length(filledMillis)
        if filledMillis(i,1) == 0
            filledMillis(i,1) = filledMillis((i-1),1);
        end
    end
    %column 11 of temptest holds poke that records initiation and not length of each poke recorded
    temptest(:,11) = initialPoke(:,1);
    %column 10 of temptest holds delay that has been filled to each sample in a trial
    temptest(:,10) = filledMillis(:,1);
    phase = temptest(:,4);
    
    %% Fill temptest to store useful information
    temptest(:,13) = 0;
    goPoke = temptest(:,13);
    %deciding whether it is a go cycle, if has a 7 and then a 3
    for i = 1:length(phase)-1
        if ((phase(i+1,1) == 3) && (phase(i,1) == 7))
            goPoke(i+1,1) = 1;
        end
    end
    %taking into account the first cycle
    if phase(1,1) == 3
        goPoke(1,1) = 1;
    end
    temptest(:,14) = 0;
    noGoPoke = temptest(:,14);
    %deciding whether it is a no go cycle, if has a 7 and then a 9
    for i = 1:length(phase)-1
        if((phase(i+1,1) == 9) && (phase(i,1) == 7))
            noGoPoke(i+1,1) = 1;
        end
    end
    %taking into account the first cycle
    if phase(1,1) == 9
        noGoPoke(1,1) = 1;
    end
    %column 13 has a 1 at the start of every water trial; 14 has a 1 at the start of every blank
    temptest(:,13) = goPoke(:,1);
    temptest(:,14) = noGoPoke(:,1);
    millis = temptest(:,1);
    %calculating the adjusted millis by adding rand delay and millis
    %adding the random delay (in s) to the time of every sample - the relevant calculation is adding to the sample at the start
    %of a trial, since it allows calculation of the start of puff
    temptest(:,15) = millis+(filledMillis/1000);
    adjustedMillis = temptest(:,15);
    temptest(:,16) = 0;
    adjustedMillisGoPhase = temptest(:,16);
    %if it is a go cycle then write the adjusted millis plus 1
    %add 1 second to the adjusted start time for only the first sample of water trials = start + rand delay + puff + delay = start of water delivery
    for i = 1:length(goPoke)
        if goPoke(i,1) == 1
            adjustedMillisGoPhase(i,1) = adjustedMillis(i,1)+1;
        end
    end
    temptest(:,16) = adjustedMillisGoPhase(:,1);
    
    %% Locate the time-stamp for last 20% trial for each training day
    %%% ALEX ADDED THIS!!!
    %%% first get both go and no go trials identified by transitions
    temptest(:,17) = 0;
    GoANDnoGoPoke = temptest(:,17);
    for i = 1:length(phase)-1
        if((phase(i+1,1) == 9) && (phase(i,1) == 7)) || ((phase(i+1,1) == 3) && (phase(i,1) == 7))
            GoANDnoGoPoke(i+1,1) = 1;
        end
    end
    %%% Including the first cycle
    if phase(1,1) == 9 || phase(1,1) == 3
        GoANDnoGoPoke(1,1) = 1;
    end
    %col 17 has 1 at the start of every trial
    temptest(:,17) = GoANDnoGoPoke(:,1);
    temptest(:,18) = 0;
    %% Create a column where the timestamp of initiation duplicates for entire block of trial
    for iii = 1:length(GoANDnoGoPoke)
        if temptest(iii,17) == 1
            temptest(iii,18) = temptest(iii,1);
        elseif temptest(iii,17) == 0
            temptest(iii,18) = temptest(iii-1,18);
        end
    end
    
    %% Copy initiation time to all times where the trial block has at least a lick in it <- WRONG

    %extends 1s in col 13 and col 14 from start of trial to every sample in trial 
    %(all samples in water trial have 1, all samples in blank trial have 0 in col 13, vice versa in 14)
    for j = 1:length(temptest(:,18))-1
        if temptest(j,18) == temptest(j+1,18)
            if temptest(j,13) == 1
                temptest(j+1,13) = 1;
            end
            if temptest(j,14) == 1
                temptest(j+1,14) = 1;
            end
            continue
        end
    end
    events_tot = length(unique(nonzeros(temptest(:,18))));
    % stores trial time starts + delay = time of start of puff across every sample in a trial
    temptest(:,21) = temptest(:,18) + temptest(:,10)/1000;
    % lick time that followed that time start + delay - 
    for licks = 1:length(temptest(:,3))
        %if there was a lick 
        if temptest(licks,3) == 2 
           %col 22 holds the time of that trial - (start of trial + random delay) = time from start of puff to lick
            temptest(licks,22) = temptest(licks,1) - (temptest(licks,18) + temptest(licks,10)/1000);
        end
    end
    % split into go versus no-go trials - select only the samples that have a 1 in col 13 (water) or 1 in col 14 (blank)
    temptest_go = temptest;
    temptest_nogo = temptest;
    ind_nogo1 = find(temptest(:,14) == 1);
    temptest_go(ind_nogo1,:) = [];
    ind_go1 = find(temptest(:,13) == 1);
    temptest_nogo(ind_go1,:) = [];

    % events_go = length(unique(nonzeros(temptest_go(:,18))));
    % events_nogo = length(unique(nonzeros(temptest_nogo(:,18))));
    
    %% Find what time the last 20% of all trials for a given day starts
    %%% split into days since start, find respective trial ranges, and get the
    %%% last 20% trials of each day (So far, it rounds to the nearest time of
    %%% the last 20%, but can change to ceil or floor depending on what you
    %%% need). You can also change the last % by changing 0.8 (1-percent).
    timestamp = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    % unique_timestamp = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    % last20prctstart_time = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    
    %col 7 is time of sample as fraction of days since start of training
    %each elt of timestamp countains the start time of all samples that occured on a given day, going from noon to noon
    %repeats start time of trial for each sample of trial
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        for idx = 1:length(temptest(:,9))
            if temptest(idx,9) >= 12+24*(d-1) && temptest(idx,9) < 12+24*(d)
                timestamp{d} = [timestamp{d};temptest(idx,18)];
            end
        end
        %if there are at least 3 trials in a day, start time of last 20% of trials is value of elt that is 80% through trials of day
        if numel(unique(timestamp{d})) >= 3
            unique_timestamp{d} = unique(timestamp{d});
            last20prctstart_time(d) = (unique_timestamp{d}(round(0.8*length(unique_timestamp{d}))+1,1));
        else
            last20prctstart_time(d) = 0;
        end
    end
    
    %% START OF LVG : lvg = water trial samples
    %%% transfer all data to mulitple columns so can split up into days
    %%% split up into days based on the hour of that event
    lvgfinal = zeros(1,6*(ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))));
    
    %for each day
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        %for all water samples
        for idx = 1:length(temptest_go(:,9))
            %if the sample was on that day (noon to noon) (col 9 is sample start time as a fraction of an hour)
            if temptest_go(idx,9) >= 12+24*(d-1) && temptest_go(idx,9) < 12+24*(d)
                %write col 1, 7, 9, 21, 18, and 22 into lvgfinal, splitting days into separate cols 
                %i.e. day 1 = cols 1-6; day 2 = cols 7-12, etc.
                lvgfinal(idx+1,1+6*(d-1):6+6*(d-1)) = temptest_go(idx,[1 7 9 21 18 22]);
            end
        end
    end
    %as far as i can tell, the code in this loop does nothing
    %for each day
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        %get the samples for that day
        hh = lvgfinal(2:end,1+6*(d-1):6+6*(d-1));
        %remove all samples that are 0
        hh(all(hh==0,2),:)=[];
        [i,ii] = size(lvgfinal(2:end,1:6));
        %replace all samples with 0
        hh(i,ii) = 0;
        %pad lvg with 0 samples so each day has same number of samples
        lvgfinal(2:end,1+6*(d-1):6+6*(d-1)) = hh(:,1:6);
    end
    
    %% START OF LVF = blank trial samples, same as water samples
    %kept columns = 1 (sample time in ms), 7 (sample time in hours), 9 (sample time as a fraction of an hour, including start time), 
    % 21 (time of puff start across every sample in a trial), 18 (start of trial in ms across every sample of trial) and 22 (time from start of puff to lick)
    %%% transfer all data to mulitple columns so can split up into days
    lvffinal = zeros(1,6*(ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))));
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        for idx = 1:length(temptest_nogo(:,9))
            if temptest_nogo(idx,9) >= 12+24*(d-1) && temptest_nogo(idx,9) < 12+24*(d)
                lvffinal(idx+1,1+6*(d-1):6+6*(d-1)) = temptest_nogo(idx,[1 7 9 21 18 22]);
            end
        end
    end
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        hh = lvffinal(2:end,1+6*(d-1):6+6*(d-1));
        hh(all(hh==0,2),:)=[];
        [i,ii] = size(lvffinal(2:end,1:6));
        hh(i,ii) = 0;
        lvffinal(2:end,1+6*(d-1):6+6*(d-1)) = hh(:,1:6);
    end
    
    %% Converting into an excel file that contained the following headers
    %%% ALEX'S MODIFICATION MOSTLY STARTS HERE!!!
    %this block of code is what spits out the warnigns- when warnings show, code has gotten at least this far
    header = {'Day 1','','','','','','Day 2','','','','','' ...
        ,'Day 3','','','','','','Day 4','','','','','' ...
        ,'Day 5','','','','','','Day 6','','','','','' ...
        ,'Day 7','','','','','','Day 8','','','','','' ...
        ,'Day 9','','','','','','Day 10','','','','',''  ...
        ,'Day 11','','','','','', 'Day 12','','','','','' ...
        , 'Day 13','','','','','', 'Day 14','','','','','' ...
        , 'Day 15','','','','','', 'Day 16','','','','','' ...
        , 'Day 17','','','','','', 'Day 18','','','','','' ...
        , 'Day 19','','','','','', 'Day 20','','','','','' };
    subheader = {'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0' ...
        'Lick Time','Days','Time','Recent T0','Next T0','L0','Lick Time','Days','Time','Recent T0','Next T0','L0'};
    lvgwheader = [header(1:length(lvgfinal(1,:)));subheader(1:length(lvgfinal(1,:))); num2cell(lvgfinal)];
    lvfwheader = [header(1:length(lvgfinal(1,:)));subheader(1:length(lvffinal(1,:))); num2cell(lvffinal)];
    xlsfilename = [textfilenamestr(1:end-23),'.xlsx'];
    xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename),lvgwheader,'water');
    xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename),lvfwheader,'no water');
    
    %% Splitting the cells by time
    %%% this also protects against trials splitting across two time bins
    lvgby4hour_cell = cell(6,(length(lvgfinal(1,:)))/6);
    lvgby4hour_cell(:,:) = {zeros(1,6)};
    lvfby4hour_cell = cell(6,(length(lvffinal(1,:)))/6);
    lvfby4hour_cell(:,:) = {zeros(1,6)};
    %splits the day split data for water trials by manually checking if a sample is in a given time bin; uses the 3rd col (sample time as fraction of an hour
    for d = 1:(length(lvgfinal(1,:)))/6
        for i = 2:length(lvgfinal)
            if lvgfinal(i,3+6*(d-1)) >= 12+24*(d-1) && lvgfinal(i,3+6*(d-1)) < 16+24*(d-1)
                lvgby4hour_cell{1,d} = [lvgby4hour_cell{1,d};lvgfinal(i,1+6*(d-1):6+6*(d-1))];
                if d > 1
                    if ~isempty(find(lvgby4hour_cell{1,d}(2,5) == lvgby4hour_cell{6,d-1}(end,5), 1))
                        idx_across = find(lvgby4hour_cell{2,d}(:,5) == lvgby4hour_cell{1,d}(2,5));
                        lvgby4hour_cell{6,d-1}(end+idx_across-1,:) = lvgby4hour_cell{1,d}(idx_across,:);
                        lvgby4hour_cell{1,d}(2:max(idx_across),:) = [];
                    end
                end
            elseif lvgfinal(i,3+6*(d-1)) >= 16+24*(d-1) && lvgfinal(i,3+6*(d-1)) < 20+24*(d-1)
                lvgby4hour_cell{2,d} = [lvgby4hour_cell{2,d};lvgfinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvgby4hour_cell{2,d}(2,5) == lvgby4hour_cell{1,d}(end,5), 1))
                    idx_across = find(lvgby4hour_cell{2,d}(:,5) == lvgby4hour_cell{2,d}(2,5));
                    lvgby4hour_cell{1,d}(end+idx_across-1,:) = lvgby4hour_cell{2,d}(idx_across,:);
                    lvgby4hour_cell{2,d}(2:max(idx_across),:) = [];
                end
            elseif lvgfinal(i,3+6*(d-1)) >= 20+24*(d-1) && lvgfinal(i,3+6*(d-1)) < 24+24*(d-1)
                lvgby4hour_cell{3,d} = [lvgby4hour_cell{3,d};lvgfinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvgby4hour_cell{3,d}(2,5) == lvgby4hour_cell{2,d}(end,5), 1))
                    idx_across = find(lvgby4hour_cell{3,d}(:,5) == lvgby4hour_cell{3,d}(2,5));
                    lvgby4hour_cell{2,d}(end+idx_across-1,:) = lvgby4hour_cell{3,d}(idx_across,:);
                    lvgby4hour_cell{3,d}(2:max(idx_across),:) = [];
                end
            elseif lvgfinal(i,3+6*(d-1)) >= 24+24*(d-1) && lvgfinal(i,3+6*(d-1)) < 28+24*(d-1)
                lvgby4hour_cell{4,d} = [lvgby4hour_cell{4,d};lvgfinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvgby4hour_cell{4,d}(2,5) == lvgby4hour_cell{3,d}(end,5), 1))
                    idx_across = find(lvgby4hour_cell{4,d}(:,5) == lvgby4hour_cell{4,d}(2,5));
                    lvgby4hour_cell{3,d}(end+idx_across-1,:) = lvgby4hour_cell{4,d}(idx_across,:);
                    lvgby4hour_cell{4,d}(2:max(idx_across),:) = [];
                end
            elseif lvgfinal(i,3+6*(d-1)) >= 28+24*(d-1) && lvgfinal(i,3+6*(d-1)) < 32+24*(d-1)
                lvgby4hour_cell{5,d} = [lvgby4hour_cell{5,d};lvgfinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvgby4hour_cell{5,d}(2,5) == lvgby4hour_cell{4,d}(end,5), 1))
                    idx_across = find(lvgby4hour_cell{5,d}(:,5) == lvgby4hour_cell{5,d}(2,5));
                    lvgby4hour_cell{4,d}(end+idx_across-1,:) = lvgby4hour_cell{5,d}(idx_across,:);
                    lvgby4hour_cell{5,d}(2:max(idx_across),:) = [];
                end
            elseif lvgfinal(i,3+6*(d-1)) >= 32+24*(d-1) && lvgfinal(i,3+6*(d-1)) < 36+24*(d-1)
                lvgby4hour_cell{6,d} = [lvgby4hour_cell{6,d};lvgfinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvgby4hour_cell{6,d}(2,5) == lvgby4hour_cell{5,d}(end,5), 1))
                    idx_across = find(lvgby4hour_cell{6,d}(:,5) == lvgby4hour_cell{6,d}(2,5));
                    lvgby4hour_cell{5,d}(end+idx_across-1,:) = lvgby4hour_cell{6,d}(idx_across,:);
                    lvgby4hour_cell{6,d}(2:max(idx_across),:) = [];
                end
            end
        end
%same as above but for blank trials
        for i = 2:length(lvffinal)
            if lvffinal(i,3+6*(d-1)) >= 12+24*(d-1) && lvffinal(i,3+6*(d-1)) < 16+24*(d-1)
                lvfby4hour_cell{1,d} = [lvfby4hour_cell{1,d};lvffinal(i,1+6*(d-1):6+6*(d-1))];
                if d > 1
                    if ~isempty(find(lvfby4hour_cell{1,d}(2,5) == lvfby4hour_cell{6,d-1}(end,5), 1))
                        idx_across = find(lvfby4hour_cell{2,d}(:,5) == lvfby4hour_cell{1,d}(2,5));
                        lvfby4hour_cell{6,d-1}(end+idx_across-1,:) = lvfby4hour_cell{1,d}(idx_across,:);
                        lvfby4hour_cell{1,d}(2:max(idx_across),:) = [];
                    end
                end
            elseif lvffinal(i,3+6*(d-1)) >= 16+24*(d-1) && lvffinal(i,3+6*(d-1)) < 20+24*(d-1)
                lvfby4hour_cell{2,d} = [lvfby4hour_cell{2,d};lvffinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvfby4hour_cell{2,d}(2,5) == lvfby4hour_cell{1,d}(end,5), 1))
                    idx_across = find(lvfby4hour_cell{2,d}(:,5) == lvfby4hour_cell{2,d}(2,5));
                    lvfby4hour_cell{1,d}(end+idx_across-1,:) = lvfby4hour_cell{2,d}(idx_across,:);
                    lvfby4hour_cell{2,d}(2:max(idx_across),:) = [];
                end
            elseif lvffinal(i,3+6*(d-1)) >= 20+24*(d-1) && lvffinal(i,3+6*(d-1)) < 24+24*(d-1)
                lvfby4hour_cell{3,d} = [lvfby4hour_cell{3,d};lvffinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvfby4hour_cell{3,d}(2,5) == lvfby4hour_cell{2,d}(end,5), 1))
                    idx_across = find(lvfby4hour_cell{3,d}(:,5) == lvfby4hour_cell{3,d}(2,5));
                    lvfby4hour_cell{2,d}(end+idx_across-1,:) = lvfby4hour_cell{3,d}(idx_across,:);
                    lvfby4hour_cell{3,d}(2:max(idx_across),:) = [];
                end
            elseif lvffinal(i,3+6*(d-1)) >= 24+24*(d-1) && lvffinal(i,3+6*(d-1)) < 28+24*(d-1)
                lvfby4hour_cell{4,d} = [lvfby4hour_cell{4,d};lvffinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvfby4hour_cell{4,d}(2,5) == lvfby4hour_cell{3,d}(end,5), 1))
                    idx_across = find(lvfby4hour_cell{4,d}(:,5) == lvfby4hour_cell{4,d}(2,5));
                    lvfby4hour_cell{3,d}(end+idx_across-1,:) = lvfby4hour_cell{4,d}(idx_across,:);
                    lvfby4hour_cell{4,d}(2:max(idx_across),:) = [];
                end
            elseif lvffinal(i,3+6*(d-1)) >= 28+24*(d-1) && lvffinal(i,3+6*(d-1)) < 32+24*(d-1)
                lvfby4hour_cell{5,d} = [lvfby4hour_cell{5,d};lvffinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvfby4hour_cell{5,d}(2,5) == lvfby4hour_cell{4,d}(end,5), 1))
                    idx_across = find(lvfby4hour_cell{5,d}(:,5) == lvfby4hour_cell{5,d}(2,5));
                    lvfby4hour_cell{4,d}(end+idx_across-1,:) = lvfby4hour_cell{5,d}(idx_across,:);
                    lvfby4hour_cell{5,d}(2:max(idx_across),:) = [];
                end
            elseif lvffinal(i,3+6*(d-1)) >= 32+24*(d-1) && lvffinal(i,3+6*(d-1)) < 36+24*(d-1)
                lvfby4hour_cell{6,d} = [lvfby4hour_cell{6,d};lvffinal(i,1+6*(d-1):6+6*(d-1))];
                if ~isempty(find(lvfby4hour_cell{6,d}(2,5) == lvfby4hour_cell{5,d}(end,5), 1))
                    idx_across = find(lvfby4hour_cell{6,d}(:,5) == lvfby4hour_cell{6,d}(2,5));
                    lvfby4hour_cell{5,d}(end+idx_across-1,:) = lvfby4hour_cell{6,d}(idx_across,:);
                    lvfby4hour_cell{6,d}(2:max(idx_across),:) = [];
                end
            end
        end
    end
    
    %% find max length to pad with zeros
    %pad each 4h bin cell with 0s so they can be concatenated columnwise into one matrix
    lvgby4hour_cell_pad = cell(6,(length(lvgfinal(1,:)))/6);
    lvgby4hour_cell_pad(:,:) = {zeros(1,6)};
    for d = 1:(length(lvgfinal(1,:)))/6
        maxlen(d) = max([length(lvgby4hour_cell{1,d}(2:end,1)),length(lvgby4hour_cell{2,d}(2:end,1)), ...
            length(lvgby4hour_cell{3,d}(2:end,1)),length(lvgby4hour_cell{4,d}(2:end,1)), ...
            length(lvgby4hour_cell{5,d}(2:end,1)),length(lvgby4hour_cell{6,d}(2:end,1))]);
        for bi = 1:6
            lvgby4hour_cell_pad{bi,d} = [lvgby4hour_cell{bi,d}(2:end,:); zeros(maxlen(d)-length(lvgby4hour_cell{bi,d}(2:end,1)),6)];
        end
    end
    lvfby4hour_cell_pad = cell(6,(length(lvffinal(1,:)))/6);
    lvfby4hour_cell_pad(:,:) = {zeros(1,6)};
    for d = 1:(length(lvffinal(1,:)))/6
        maxlenf(d) = max([length(lvfby4hour_cell{1,d}(2:end,1)),length(lvfby4hour_cell{2,d}(2:end,1)), ...
            length(lvfby4hour_cell{3,d}(2:end,1)),length(lvfby4hour_cell{4,d}(2:end,1)), ...
            length(lvfby4hour_cell{5,d}(2:end,1)),length(lvfby4hour_cell{6,d}(2:end,1))]);
        for bi = 1:6
            lvfby4hour_cell_pad{bi,d} = [lvfby4hour_cell{bi,d}(2:end,:); zeros(maxlenf(d)-length(lvfby4hour_cell{bi,d}(2:end,1)),6)];
        end
    end
    
    %% Concatenate the four hour bins
    lvg4hrbin_days = cell(1,(length(lvgfinal(1,:)))/6);
    lvf4hrbin_days = cell(1,(length(lvffinal(1,:)))/6);
    %stick 4h bins together columnwise, the same way days were separated earlier, with columns of nans to separate bins
    for d = 1:(length(lvffinal(1,:)))/6
        lvg4hrbin_days{d} = cat(2, lvgby4hour_cell_pad{1,d},nan(maxlen(d),2),lvgby4hour_cell_pad{2,d},...
            nan(maxlen(d),2),lvgby4hour_cell_pad{3,d}, nan(maxlen(d),2),lvgby4hour_cell_pad{4,d},...
            nan(maxlen(d),2),lvgby4hour_cell_pad{5,d},nan(maxlen(d),2),lvgby4hour_cell_pad{6,d},...
            nan(maxlen(d),2));
        lvf4hrbin_days{d} = cat(2, lvfby4hour_cell_pad{1,d},nan(maxlenf(d),2),lvfby4hour_cell_pad{2,d},...
            nan(maxlenf(d),2),lvfby4hour_cell_pad{3,d}, nan(maxlenf(d),2),lvfby4hour_cell_pad{4,d},...
            nan(maxlenf(d),2),lvfby4hour_cell_pad{5,d},nan(maxlenf(d),2),lvfby4hour_cell_pad{6,d},...
            nan(maxlenf(d),2));
    end

    %% calculate licking  frequency (?)
    lvg_lckpk = cell(1,(length(lvgfinal(1,:)))/6);
    lvf_lckpk = cell(1,(length(lvffinal(1,:)))/6);
    %for each day
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        %if there are trials on that day
        if ~isempty(lvg4hrbin_days{d})
            %for each 4hr bin
            %k = 5 = start time  of trial at all samples in bin
            %k + 1 = col 22 of temptest = time after start of puff to lick
            for k = 5:8:length(lvg4hrbin_days{d}(1,:))
                %number of trials in that bin (water)
                lvg_lckpk{d}(1,k-3) = length(nonzeros(unique(lvg4hrbin_days{d}(:,k))));
                %number of samples with a lick between 0.7s and 1s after start of puff (anticipatory licking window) (water)
                lvg_lckpk{d}(1,k-4) = length(find(lvg4hrbin_days{d}(:,k+1) >= 0.7 & lvg4hrbin_days{d}(:,k+1)<=1));
                %number of anticipatory licks/number of trials = anticipatory licking frequency ???????
                lvg_lckpk{d}(1,k-2) = lvg_lckpk{d}(:,k-4)/lvg_lckpk{d}(:,k-3);
                %same as above for blank trials
                lvf_lckpk{d}(1,k-3) = length(nonzeros(unique(lvf4hrbin_days{d}(:,k))));
                lvf_lckpk{d}(1,k-4) = length(find(lvf4hrbin_days{d}(:,k+1) >= 0.7 & lvf4hrbin_days{d}(:,k+1)<=1));
                lvf_lckpk{d}(1,k-2) = lvf_lckpk{d}(:,k-4)/lvf_lckpk{d}(:,k-3);
            end
        end
        endf
    
    %% separately isolate last 20% for go and no go
    lvg_lckpk_lst20 = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    %for water trials, for each day
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        %if there are samples in that day
        if ~isempty(lvg4hrbin_days{d})
            %for each 4h bin
            for kk = 1:8:length(lvg4hrbin_days{d}(1,:))
                %for each sample
                for idx_2 = 1:length(lvg4hrbin_days{d}(:,1))
                    %if there are enough trials for the last 20% to be set and the sample start time is in the last 20%
                    if last20prctstart_time(d) > 0
                        if ~isempty(find(any(lvg4hrbin_days{d}(:)>=last20prctstart_time(d)), 1))
                            if lvg4hrbin_days{d}(idx_2,kk) >= last20prctstart_time(d)
                                %append the relevant sample as a cell
                                lvg_lckpk_lst20{d} = [lvg_lckpk_lst20{d};num2cell(lvg4hrbin_days{d}(idx_2,kk:kk+5))];
                            end
                        else
                            %otherwise append NaNs
                            lvg_lckpk_lst20{d} = num2cell([nan(1,6)]);
                        end
                    else
                        lvg_lckpk_lst20{d} = num2cell([nan(1,6)]);
                    end
                end
            end
            %taking only the samples in the anticipatory licking window (between .7s and 1s after puff start)
            lvg_lckpk_lst20{d} = cat(1,[num2cell(length(find(cell2mat(lvg_lckpk_lst20{d}(:,6)) >= 0.7 & cell2mat(lvg_lckpk_lst20{d}(:,6))<=1))),...
                num2cell(length(nonzeros(~isnan(unique(cell2mat(lvg_lckpk_lst20{d}(:,5))))))),...
                num2cell(length(find(cell2mat(lvg_lckpk_lst20{d}(:,6)) >= 0.7 & cell2mat(lvg_lckpk_lst20{d}(:,6))<=1))/ ...
                length(nonzeros(~isnan(unique(cell2mat(lvg_lckpk_lst20{d}(:,5))))))), ...
                num2cell(nan(1,3))],lvg_lckpk_lst20{d});
        end
    end
    %same as above for blank trials
    lvf_lckpk_lst20 = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        if ~isempty(lvf4hrbin_days{d})
            for kk = 1:8:length(lvf4hrbin_days{d}(1,:))
                for idx_2 = 1:length(lvf4hrbin_days{d}(:,1))
                    if last20prctstart_time(d) > 0
                        if ~isempty(find(any(lvf4hrbin_days{d}(:)>=last20prctstart_time(d)), 1))
                            if lvf4hrbin_days{d}(idx_2,kk) >= last20prctstart_time(d)
                                lvf_lckpk_lst20{d} = [lvf_lckpk_lst20{d};num2cell(lvf4hrbin_days{d}(idx_2,kk:kk+5))];
                            end
                        else
                            lvf_lckpk_lst20{d} = num2cell([nan(1,6)]);
                        end
                    else
                        lvf_lckpk_lst20{d} = num2cell([nan(1,6)]);
                    end
                end
            end
            lvf_lckpk_lst20{d} = cat(1,[num2cell(length(find(cell2mat(lvf_lckpk_lst20{d}(:,6)) >= 0.7 & cell2mat(lvf_lckpk_lst20{d}(:,6))<=1))),...
                num2cell(length(nonzeros(~isnan(unique(cell2mat(lvf_lckpk_lst20{d}(:,5))))))),...
                num2cell(length(find(cell2mat(lvf_lckpk_lst20{d}(:,6)) >= 0.7 & cell2mat(lvf_lckpk_lst20{d}(:,6))<=1))/ ...
                length(nonzeros(~isnan(unique(cell2mat(lvf_lckpk_lst20{d}(:,5))))))), ...
                num2cell(nan(1,3))],lvf_lckpk_lst20{d});
        end
    end
  %adding NaNs to the end of the licks/poke data?
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        if ~isempty(lvg4hrbin_days{d})
            lvg_lckpk{d} = [lvg_lckpk{d}, nan(1,length(lvg4hrbin_days{d}(1,:))-length(lvg_lckpk{d}(1,:)))];
        end
        if ~isempty(lvf4hrbin_days{d})
            lvf_lckpk{d} = [lvf_lckpk{d}, nan(1,length(lvf4hrbin_days{d}(1,:))-length(lvf_lckpk{d}(1,:)))];
        end
    end
    
    %% Write the events, poke number, and lick/poke into excel file splitting by days and go vs nogo
    header_bin = {'Total Events','Poke Number','Licks/Poke','','','', ...
        '','','Total Events','Poke Number','Licks/Poke','','','', ...
        '','','Total Events','Poke Number','Licks/Poke','','','', ...
        '','','Total Events','Poke Number','Licks/Poke','','','', ...
        '','','Total Events','Poke Number','Licks/Poke','','','', ...
        '','','Total Events','Poke Number','Licks/Poke','','','', ...
        '','','Total Events','Poke Number','Licks/Poke','','',''};
    subheader_bin = {'12:00 - 15:59','','','','','', ...
        '','','16:00 - 19:59','','','','','', ...
        '','','20:00 - 23:59','','','','','', ...
        '','','00:00 - 03:59','','','','','', ...
        '','','04:00 - 07:59','','','','','', ...
        '','','08:00 - 11:59','','','','','', ...
        '','','LAST 20% GO + NO GO','','','','',''};
    lvg4hrbinwheader = {};
    lvf4hrbinwheader = {};
    %set up the water trials to write to excel
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        if ~isempty(lvg4hrbin_days{d})
            if length(lvg4hrbin_days{d}(:,1))-length(lvg_lckpk_lst20{d}(2:end,1)) >= 0
                lvg4hrbinwheader{d} = [header_bin; ...
                    num2cell([lvg_lckpk{d},cell2mat(lvg_lckpk_lst20{d}(1,:))]);...
                    subheader_bin;num2cell([lvg4hrbin_days{d},...
                    [cell2mat(lvg_lckpk_lst20{d}(2:end,:)); ...
                    nan(abs(length(lvg4hrbin_days{d}(:,1))-length(lvg_lckpk_lst20{d}(2:end,1))),6)]])];
            elseif length(lvg4hrbin_days{d}(:,1))-length(lvg_lckpk_lst20{d}(2:end,1)) < 0
                lvg4hrbinwheader{d} = [header_bin; ...
                    num2cell([lvg_lckpk{d},cell2mat(lvg_lckpk_lst20{d}(1,:))]);...
                    subheader_bin;num2cell([[lvg4hrbin_days{d};...
                    nan(abs(length(lvg4hrbin_days{d}(:,1))-length(lvg_lckpk_lst20{d}(2:end,1))),48)],...
                    [cell2mat(lvg_lckpk_lst20{d}(2:end,:))]])];
            end
        else
            lvg4hrbinwheader{d} = [header_bin;num2cell(nan(1,length(header_bin)));subheader_bin];
        end
        %set up the blank trials to write to excel
        if ~isempty(lvf4hrbin_days{d})
            if length(lvf4hrbin_days{d}(:,1))-length(lvf_lckpk_lst20{d}(2:end,1)) >= 0
                lvf4hrbinwheader{d} = [header_bin; ...
                    num2cell([lvf_lckpk{d},cell2mat(lvf_lckpk_lst20{d}(1,:))]);...
                    subheader_bin;num2cell([lvf4hrbin_days{d},...
                    [cell2mat(lvf_lckpk_lst20{d}(2:end,:)); ...
                    nan(abs(length(lvf4hrbin_days{d}(:,1))-length(lvf_lckpk_lst20{d}(2:end,1))),6)]])];
            elseif length(lvf4hrbin_days{d}(:,1))-length(lvf_lckpk_lst20{d}(2:end,1)) < 0
                lvf4hrbinwheader{d} = [header_bin; ...
                    num2cell([lvf_lckpk{d},cell2mat(lvf_lckpk_lst20{d}(1,:))]);...
                    subheader_bin;num2cell([[lvf4hrbin_days{d};...
                    nan(abs(length(lvf4hrbin_days{d}(:,1))-length(lvf_lckpk_lst20{d}(2:end,1))),48)],...
                    [cell2mat(lvf_lckpk_lst20{d}(2:end,:))]])];
            end
        else
            lvf4hrbinwheader{d} = [header_bin;num2cell(nan(1,length(header_bin)));subheader_bin];
        end
    end
    %write data to excel
    %%% d1W, d2W, d3W ... followed by d1N, d2N, d3N ...
    xlsfilename1 = [textfilenamestr(1:end-23),'_DataFinal.xlsx'];
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename1),lvg4hrbinwheader{d}, ...
            sprintf('%s%s%s','d',num2str(d),'W'));
    end
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename1),lvf4hrbinwheader{d}, ...
            sprintf('%s%s%s','d',num2str(d),'N'));
    end
    
    %% for isolating the events, number of pokes, and lick/poke for a separate excel
    lickpoke_time = {};
    lickpoke_water = {};
    lickpoke_nowater = {};
    %for each day
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        lickpoke_time = [cellstr([subheader_bin{1};subheader_bin{9};subheader_bin{17};...
            subheader_bin{25};subheader_bin{33};subheader_bin{41}]);num2cell(nan(1,1));cellstr(subheader_bin{49})];
        if ~isempty(lvg4hrbin_days{d})
            %splitting the num trials/num licks/licks per trial data by 4h time bin for water trials
            lickpoke_water{d} = [num2cell(lvg_lckpk{d}(1:3));num2cell(lvg_lckpk{d}(9:11));...
                num2cell(lvg_lckpk{d}(17:19));num2cell(lvg_lckpk{d}(25:27));...
                num2cell(lvg_lckpk{d}(33:35));num2cell(lvg_lckpk{d}(41:43));num2cell(nan(1,3));lvg_lckpk_lst20{d}(1,1:3)];
        else
            lickpoke_water{d} = nan;
        end
        %same for blank trials
        if ~isempty(lvf4hrbin_days{d})
            lickpoke_nowater{d} = [num2cell(lvf_lckpk{d}(1:3));num2cell(lvf_lckpk{d}(9:11));...
                num2cell(lvf_lckpk{d}(17:19));num2cell(lvf_lckpk{d}(25:27));...
                num2cell(lvf_lckpk{d}(33:35));num2cell(lvf_lckpk{d}(41:43));num2cell(nan(1,3));lvf_lckpk_lst20{d}(1,1:3)];
        else
            lickpoke_nowater{d} = nan;
        end
    end
    %% Individual trial structure and instantaneous lick frequencies for last 20% of each day
    lvg_last20_trial_idx = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    idx_go_trialstart = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvg_indtrials = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvg_lst20_antlck_freq = {};
    lvf_last20_trial_idx = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    idx_nogo_trialstart = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvf_indtrials = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvf_lst20_antlck_freq = {};
    %for each day
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        %% go (water)
        %if there are no trials in last 20%
        if isempty(cell2mat(lvg_lckpk_lst20{d})) == 1 || isempty(cell2mat(lvf_lckpk_lst20{d})) == 1
            lvg_lst20_antlck_freq{d} = NaN;
            lvf_lst20_antlck_freq{d} = NaN;
        
        elseif cell2mat(lvg_lckpk_lst20{d}(1,2)) > 0 && cell2mat(lvf_lckpk_lst20{d}(1,2)) > 0
            lvg_last20_trial_idx{d}(2,1) = 1;
            for j = 1:length(lvg_lckpk_lst20{d}(:,5))-1
                if cell2mat(lvg_lckpk_lst20{d}(j,5)) < cell2mat(lvg_lckpk_lst20{d}(j+1,5))
                    lvg_last20_trial_idx{d}(j+1,1) = 1;
                end
            end
            idx_go_trialstart{d} = find(lvg_last20_trial_idx{d} == 1);
            idx_go_trialstart{d} = [idx_go_trialstart{d};length(lvg_lckpk_lst20{d})];
            for ii = 1:length(idx_go_trialstart{d})-1
                %breaking individual trials into cells
                lvg_indtrials{ii,d} = cell2mat(lvg_lckpk_lst20{d}(idx_go_trialstart{d}(ii):idx_go_trialstart{d}(ii+1)-1,1:6));
            end
            lvg_lst20_antlck_cts = [];
            for i = 1:ii
                %% calculation of licking frequency for last 20% of trials
                %finding samples in anticipatory window from each indiviudual trial
                lvg_lst20_antlck_cts(i) = length(find(lvg_indtrials{i,d}(:,6) >= 0.7 & lvg_indtrials{i,d}(:,6)<=1));
                lvg_lst20_freq(i,d) = lvg_lst20_antlck_cts(i)/0.3;
            end
            lvg_lst20_antlck_freq{d} = lvg_lst20_antlck_cts/0.3;
            %% no go
            lvf_last20_trial_idx{d}(2,1) = 1;
            for j = 1:length(lvf_lckpk_lst20{d}(:,5))-1
                if cell2mat(lvf_lckpk_lst20{d}(j,5)) < cell2mat(lvf_lckpk_lst20{d}(j+1,5))
                    lvf_last20_trial_idx{d}(j+1,1) = 1;
                end
            end
            idx_nogo_trialstart{d} = find(lvf_last20_trial_idx{d} == 1);
            idx_nogo_trialstart{d} = [idx_nogo_trialstart{d};length(lvf_lckpk_lst20{d})];
            for ii = 1:length(idx_nogo_trialstart{d})-1
                lvf_indtrials{ii,d} = cell2mat(lvf_lckpk_lst20{d}(idx_nogo_trialstart{d}(ii):idx_nogo_trialstart{d}(ii+1)-1,1:6));
            end
            lvf_lst20_antlck_cts = [];
            for i = 1:ii
                lvf_lst20_antlck_cts(i) = length(find(lvf_indtrials{i,d}(:,6) >= 0.7 & lvf_indtrials{i,d}(:,6)<=1));
                lvf_lst20_freq(i,d) = lvf_lst20_antlck_cts(i)/0.3;
            end
            lvf_lst20_antlck_freq{d} = lvf_lst20_antlck_cts/0.3;
        end
    end
    
    lvglvfall = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvglvf_lckfreq_all = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvg_indtrials_marked = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    lvf_indtrials_marked = cell(1,ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7)))));
    header_days = cell(1,1);
    header_w = cell(1,1);
    header_n = cell(1,1);
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        if any(~cellfun(@isempty,lvf_indtrials(:,d))) > 0
            header_final = {'T0 for trial','Individual Trial Anticipatory Lick Frequency in Hz ((Licks/Poke)/0.3)',...
                'Time Stamp','Days','Hours','Next T0','T0','L0','Trial Type(Go = 2, NoGo = 1)',...
                'T0 for trial','Individual Trial Anticipatory Lick Frequency in Hz ((Licks/Poke)/0.3)',...
                'Time Stamp','Days','Hours','Next T0','T0','L0','Trial Type(Go = 2, NoGo = 1)'};
            header_final2 = {'Water','No Water'};
            for i = 1:length(nonzeros(~cellfun(@isempty,lvg_indtrials(:,d))))
                lvg_indtrials_marked{i,d} = cat(2,[[lvg_indtrials{i,d}(1,1)-1e-6,lvg_lst20_antlck_freq{d}(i)];...
                    NaN(length(lvg_indtrials{i,d}(:,1))-1,2)],...
                    (cat(2,lvg_indtrials{i,d},2*ones(length(lvg_indtrials{i,d}(:,1)),1))));
            end
            for i = 1:length(nonzeros(~cellfun(@isempty,lvf_indtrials(:,d))))
                lvf_indtrials_marked{i,d} = cat(2,[[lvf_indtrials{i,d}(1,1)-1e-6,lvf_lst20_antlck_freq{d}(i)];...
                    NaN(length(lvf_indtrials{i,d}(:,1))-1,2)],...
                    (cat(2,lvf_indtrials{i,d},ones(length(lvf_indtrials{i,d}(:,1)),1))));
            end
            
            lvglvfall{d} = cat(1,lvg_indtrials_marked{:,d},lvf_indtrials_marked{:,d});
            lvglvfall{d} = sortrows(lvglvfall{d},3);
            lvglvfall{d}(:,10:18) = lvglvfall{d};
            ind_nogo_lst20 = find(lvglvfall{d}(:,9) == 1);
            lvglvfall{d}(ind_nogo_lst20,1:9) = NaN;
            ind_go_lst20 = find(lvglvfall{d}(:,9) == 2);
            lvglvfall{d}(ind_go_lst20,10:18) = NaN;
            lvglvfall_output{d} = [header_final;num2cell(lvglvfall{d})];
            
            if length(lvg_lst20_antlck_freq{d}) >= length(lvf_lst20_antlck_freq{d})
                lvglvf_lckfreq_all{d} = [lvg_lst20_freq(1:length(lvg_lst20_antlck_freq{d}),d),...
                    [lvf_lst20_freq(1:length(lvf_lst20_antlck_freq{d}),d);...
                    NaN(length(lvg_lst20_antlck_freq{d})-length(lvf_lst20_antlck_freq{d}),1)]];
            else
                lvglvf_lckfreq_all{d} = [[lvg_lst20_freq(1:length(lvg_lst20_antlck_freq{d}),d);...
                    NaN(length(lvf_lst20_antlck_freq{d})-length(lvg_lst20_antlck_freq{d}),1)],...
                    lvf_lst20_freq(1:length(lvf_lst20_antlck_freq{d}),d)];
            end
            lvglvf_lckfreq_alloutput{d} = [header_final2;num2cell(lvglvf_lckfreq_all{d})];
            freq0w(d) = length(find(lvg_lst20_freq(1:length(lvg_lst20_antlck_freq{d}),d) == 0))/length(lvg_lst20_antlck_freq{d});
            freq3w(d) = length(find(lvg_lst20_freq(1:length(lvg_lst20_antlck_freq{d}),d) == 1/0.3))/length(lvg_lst20_antlck_freq{d});
            freq6w(d) = length(find(lvg_lst20_freq(1:length(lvg_lst20_antlck_freq{d}),d) == 2/0.3))/length(lvg_lst20_antlck_freq{d});
            freq10w(d) = length(find(lvg_lst20_freq(1:length(lvg_lst20_antlck_freq{d}),d) >= 10))/length(lvg_lst20_antlck_freq{d});
            freq0n(d) = length(find(lvf_lst20_freq(1:length(lvf_lst20_antlck_freq{d}),d) == 0))/length(lvf_lst20_antlck_freq{d});
            freq3n(d) = length(find(lvf_lst20_freq(1:length(lvf_lst20_antlck_freq{d}),d) == 1/0.3))/length(lvf_lst20_antlck_freq{d});
            freq6n(d) = length(find(lvf_lst20_freq(1:length(lvf_lst20_antlck_freq{d}),d) == 2/0.3))/length(lvf_lst20_antlck_freq{d});
            freq10n(d) = length(find(lvf_lst20_freq(1:length(lvf_lst20_antlck_freq{d}),d) >= 10))/length(lvf_lst20_antlck_freq{d});
        else
            header_final2 = {'Water','No Water'};
            lvglvf_lckfreq_alloutput{d} = [header_final2;num2cell(NaN(1,2))];
            freq0w(d) = NaN;
            freq3w(d) = NaN;
            freq6w(d) = NaN;
            freq10w(d) = NaN;
            freq0n(d) = NaN;
            freq3n(d) = NaN;
            freq6n(d) = NaN;
            freq10n(d) = NaN;
        end
        header_lickprob = {'No Lick Probability';'1 Lick Probability';'2 Licks Probability';'3 Licks Probability'};
        header_days = [header_days,{sprintf('%s%s','Day ',num2str(d))}];
        header_w = [header_w,{'Water'}];
        header_n = [header_n,{'No Water'}];
    end
    header_3 = [header_w,header_n];
    header_days2 = [header_days,header_days];
    allfreqw = cat(1,freq0w,freq3w,freq6w,freq10w);
    allfreqn = cat(1,freq0n,freq3n,freq6n,freq10n);
    allfreq_output = [[header_3;header_days2];[[header_lickprob,num2cell(allfreqw)],[header_lickprob,num2cell(allfreqn)]]];
    xlsfilename3 = [textfilenamestr(1:end-23),'_Last20Data.xlsx'];
    
    for d = 1:length(lvglvfall_output)
        if ~isempty(lvglvfall_output{d}) || ~isempty(lvglvfall_output{d})
            xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename3),lvglvfall_output{d}, ...
                sprintf('%s%s','Last 20% Trial-by-Trial day_',num2str(d)));
            xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename3),lvglvf_lckfreq_alloutput{d}, ...
                sprintf('%s%s','Last 20% Ant. Lick Freq. day_',num2str(d)));
        else
            xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename3),{'CANNOT COMPUTE! ONE OR MORE TRIAL TYPE MISSING.'}, ...
                sprintf('%s%s','Last 20% Trial-by-Trial day_',num2str(d)));
            xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename3),{'CANNOT COMPUTE! ONE OR MORE TRIAL TYPE MISSING.'}, ...
                sprintf('%s%s','Last 20% Ant. Lick Freq. day_',num2str(d)));
        end
    end
    
    xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename3),allfreq_output, 'L 20% Ant. Lick Ct Contingency');
    
    
    %% for organizing the output into desired format
    % MIGHT NEED TO FIX
    lvglvf_finaloutput1 = {};
    for d = 1:ceil(min(max(temptest_go(:,7)),max(temptest_nogo(:,7))))
        header_final = {'','Water','','','','No Water','','','','Total Trials',...
            '','Lick Frequency in Hz ((Licks/Poke)/0.3)','','','Performance'};
        subheader_final = {sprintf('%s%s','Day',num2str(d)),'Total Events','Poke Number','Licks/Poke',...
            '','Total Events','Poke Number','Licks/Poke','', ...
            'Water + No Water','','Water','No Water','',...
            'Water - No Water'};
        lvglvf_finaloutput = {};
        if ~isempty(lvg4hrbin_days{d}) && ~isempty(lvf4hrbin_days{d})
            lvglvf_finaloutput = [lvglvf_finaloutput;lickpoke_time,lickpoke_water{d},...
                num2cell(nan(8,1)),lickpoke_nowater{d},num2cell(nan(8,1))];
            lvglvf_finaloutput = [lvglvf_finaloutput,num2cell(cell2mat(lickpoke_water{d}(:,2))+ ...
                cell2mat(lickpoke_nowater{d}(:,2))),num2cell(nan(8,1))];
            lvglvf_finaloutput = [lvglvf_finaloutput,num2cell(cell2mat(lickpoke_water{d}(:,3))/0.3),...
                num2cell(cell2mat(lickpoke_nowater{d}(:,3))/0.3),num2cell(nan(8,1))];
            lvglvf_finaloutput = [lvglvf_finaloutput,num2cell(cell2mat(lickpoke_water{d}(:,3))/0.3 - ...
                cell2mat(lickpoke_nowater{d}(:,3))/0.3)];
            lvglvf_finaloutput = [lvglvf_finaloutput;num2cell(nan(2,15))];
        else
            lvglvf_finaloutput = [num2cell(nan(10,15))];
        end
        lvglvf_finaloutput = [header_final;subheader_final;lvglvf_finaloutput];
        lvglvf_finaloutput1 = [lvglvf_finaloutput1;lvglvf_finaloutput];
    end
    % for it = 1:d
    %     for idx = [3:8 10]
    %         if isnan(cell2mat(lvglvf_finaloutput1(idx+((it-1)*12),15)))
    %             if ~isnan(cell2mat(lvglvf_finaloutput1(idx+((it-1)*12),12)))
    %             lvglvf_finaloutput1(idx+((it-1)*12),15) = lvglvf_finaloutput1(idx+((it-1)*12),12);
    %             elseif ~isnan(cell2mat(lvglvf_finaloutput1(idx+((it-1)*12),13)))
    %             lvglvf_finaloutput1(idx+((it-1)*12),15) = num2cell(-cell2mat(lvglvf_finaloutput1(idx+((it-1)*12),13)));
    %             else
    %             lvglvf_finaloutput1(idx+((it-1)*12),15) = num2cell(0);
    %             end
    %         end
    %     end
    % end
    
    xlsfilename2 = [textfilenamestr(1:end-23),'_AnalysisFinal.xlsx'];
    xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename2),lvglvf_finaloutput1);
    
    %% For plotting lick frequency, performance and trial data as  graph, and creating necessary excel file for origin plotting
    plot_lvglvf_finaloutput = lvglvf_finaloutput1;
    finalmatforplot = [];
    headerforplot = {'Time Bin (hr)', 'Water Lick Freq (Hz)', 'NoWater Lick Freq (Hz)', 'Performance', 'Trial No.'};
    for a = 0:(d-1)
        golickfreqplot = cell2mat(plot_lvglvf_finaloutput(((a*12) + 3):((a*12)+8), 12));
        nogolickfreqplot = cell2mat(plot_lvglvf_finaloutput(((a*12) + 3):((a*12)+8), 13));
        performlickfreqplot = cell2mat(plot_lvglvf_finaloutput(((a*12) + 3):((a*12)+8), 15));
        trialnumplot = cell2mat(plot_lvglvf_finaloutput(((a*12) + 3):((a*12)+8), 10));
        timhrbinplot = transpose([(24*(a+1-accdays)-22) : 4 : (24*(a+1-accdays)-2)]);
        matforplot = [timhrbinplot, golickfreqplot, nogolickfreqplot, performlickfreqplot, trialnumplot];
        finalmatforplot = [finalmatforplot; matforplot];
    end
    plotoutputwithheader = [headerforplot; num2cell(finalmatforplot)];
    xlsfilename4 = [textfilenamestr(1:end-23),'_ForPlotting.xlsx'];
    xlswrite(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename4),plotoutputwithheader);
    
    majorticksonx = [(-24*accdays):12:(24*(d-accdays))];
    minorticksonx = [(-24*accdays):4:(24*(d-accdays))];
        
    plot1 = subplot(3,1,1);
    plot(finalmatforplot(:, 1), finalmatforplot(:, 2), 'g-o','LineWidth',1.5,'MarkerSize',4,...
        'MarkerEdgeColor','g','MarkerFaceColor','g');
    hold on;
    plot(finalmatforplot(:, 1), finalmatforplot(:, 3), 'r-o','LineWidth',1.5,'MarkerSize',4,...
        'MarkerEdgeColor','r','MarkerFaceColor','r');
    xlabel('Training time (hr)');
    ylabel('Licking (Hz)');
    xlim([((-24*accdays)-4) (24*(d-accdays)+4)]);
    ylim([0 10]); %
    xticks(majorticksonx);
    set(gca,'box','off')
    plot1.TickDir = 'out';
    plot1.XMinorTick = 'on';
    plot1.XAxis.MinorTickValues = minorticksonx;
    hold off;
    
    plot2 = subplot(3,1,2);
    h = area(finalmatforplot(:, 1), finalmatforplot(:, 5));
    h(1).FaceColor = [0.8 0.8 0.8];
    xlabel('Training time (hr)');
    ylabel('Trial number');
    xlim([((-24*accdays)-4) (24*(d-accdays)+4)]);
    ylim([0 500]);
    xticks(majorticksonx);
    yticks([0 125 250 375 500]);
    set(gca,'box','off')
    plot2.TickDir = 'out';
    plot2.XMinorTick = 'on';
    plot2.XAxis.MinorTickValues = minorticksonx;
    
    plot3 = subplot(3,1,3);
    zeroperformancex = ([((-24*accdays)-4) (24*(d-accdays)+4)]).';
    zeroperformancey = zeros(size(zeroperformancex));
    plot(zeroperformancex, zeroperformancey, 'b:','LineWidth',0.5);
    hold on;
    plot(finalmatforplot(:, 1), finalmatforplot(:, 4), 'k-o','LineWidth',1.5,'MarkerSize',4,...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    xlabel('Training time (hr)');
    ylabel('Performance');
    xlim([((-24*accdays)-4) (24*(d-accdays)+4)]);
    ylim([-10 10]);
    xticks(majorticksonx);
    yticks([-10 -5 0 5 10]);
    set(gca,'box','off')
    plot3.TickDir = 'out';
    plot3.XMinorTick = 'on';
    plot3.XAxis.MinorTickValues = minorticksonx;
    hold off;
    
    figfilename4 = [textfilenamestr(1:end-23),'_Plots.fig'];
    fullplotfigname4 = strcat(filepathOutResults, '\', figfilename4);
    savefig(fullplotfigname4);
    %% Color code the Excel files
    Excel = actxserver('excel.application');
    WB = Excel.Workbooks.Open(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename1));
    for it = 2:d*2+1
        WB.Worksheets.Item(it).Range('A3:F3').Interior.ColorIndex = 36;
        WB.Worksheets.Item(it).Range('I3:N3').Interior.ColorIndex = 15;
        WB.Worksheets.Item(it).Range('Q3:V3').Interior.ColorIndex = 40;
        WB.Worksheets.Item(it).Range('Y3:AD3').Interior.ColorIndex = 34;
        WB.Worksheets.Item(it).Range('AG3:AL3').Interior.ColorIndex = 38;
        WB.Worksheets.Item(it).Range('AO3:AT3').Interior.ColorIndex = 35;
        WB.Worksheets.Item(it).Range('AW3:BB3').Interior.ColorIndex = 24;
    end
    WB.Save();WB.Close();
    
    WB = Excel.Workbooks.Open(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename2));
    
    for it = 1:d
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(3+((it-1)*12))),':',sprintf('%s%s','O',num2str(3+((it-1)*12))))).Interior.ColorIndex = 36;
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(4+((it-1)*12))),':',sprintf('%s%s','O',num2str(4+((it-1)*12))))).Interior.ColorIndex = 15;
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(5+((it-1)*12))),':',sprintf('%s%s','O',num2str(5+((it-1)*12))))).Interior.ColorIndex = 40;
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(6+((it-1)*12))),':',sprintf('%s%s','O',num2str(6+((it-1)*12))))).Interior.ColorIndex = 34;
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(7+((it-1)*12))),':',sprintf('%s%s','O',num2str(7+((it-1)*12))))).Interior.ColorIndex = 38;
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(8+((it-1)*12))),':',sprintf('%s%s','O',num2str(8+((it-1)*12))))).Interior.ColorIndex = 35;
        WB.Worksheets.Item(1).Range(cat(2,sprintf('%s%s','A',num2str(10+((it-1)*12))),':',sprintf('%s%s','O',num2str(10+((it-1)*12))))).Interior.ColorIndex = 24;
    end
    WB.Save();WB.Close();
    
    WB = Excel.Workbooks.Open(sprintf('%s%s%s',filepathOutResults,'/',xlsfilename3));
    for it = 2:2:2*length(lvglvfall_output)+1
        WB.Worksheets.Item(it).Range('A:I').Interior.ColorIndex = 35;
        WB.Worksheets.Item(it).Range('J:R').Interior.ColorIndex = 38;
    end
    for it = 3:2:2*length(lvglvfall_output)+1                   %sum(~cellfun(@isempty,lvglvfall_output),2)+1
        WB.Worksheets.Item(it).Range('A1:A500').Interior.ColorIndex = 35;
        WB.Worksheets.Item(it).Range('B1:B500').Interior.ColorIndex = 38;
    end
    % for it = it+1
    %     WB.Worksheets.Item(it).Range(strcat('A:',string(char('A' + (length(d)+ 1))))).Interior.ColorIndex = 35;
    %     WB.Worksheets.Item(it).Range(strcat(string(char('A' + (length(d)+ 2))),':',string(char('A' + (2*length(d)+ 2))))).Interior.ColorIndex = 38;
    % end
    WB.Save();WB.Close();
    
    Excel.Quit();
end
