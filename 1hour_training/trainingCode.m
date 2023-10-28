%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Block SAT Training Code

% Created By: Julia Loghinov
% Based On: "SomatosensoryLearning_LickAnalysisV16.m"
% Last Editted: 1/23/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%
% Loading and Initializing Data
%%%%%%%%%%%%%%%%%%%%
% load training data and rename each file as the following:
% FILL IN THIS BLOCK EVERY TIME WITH SPECIFIC DATA

mouseName = "OEI3";
cageName = "Chicago";
data = {Untitled}
days = ["80%"];
colors = ["#ffc75e", "#f5942c", "#cce2ff", "#70aeff", "#0060de","#003273"];

%%%%%%%%%%%%%%%%%%%%
% Now the code will take the trail numbers and plot trials across time
%%%%%%%%%%%%%%%%%%%%
totTrialNum = [];

f=figure;
f.Position = [200, 800, 2000, 275];
for t = 1:length(data)
    temptest = data{t};
    
    time = table2array(temptest(:,1));
    time = time-min(time);
    water = table2array(temptest(:,2));
    lick = table2array(temptest(:,3));
    block = table2array(temptest(:,4));
    delay = table2array(temptest(:,5));

    % List of ALL trail start indices
    TrialStarts = [];
    for i = 2:length(block)
        if (block(i)==3 && block(i-1)~=3) || (block(i)==9 && block(i-1)~=9)
            % i is the index of the trial start in the data, but we also need
            % the actual time
            startTime = time(i);
            TrialStarts = [TrialStarts; i startTime];
        end
    end

    % Based on when the trials were started, we can see how many trials were
    % completed in set time bins

    timeBin = 300; %seconds
    numTrials = [];
    numBins = ceil(max(time)/timeBin);
    for step = 1:numBins
        trialCount = sum(TrialStarts(:,2)>=((step-1)*timeBin) & TrialStarts(:,2)<step*timeBin);
        numTrials = [numTrials; trialCount];
    end
    
    totTrialNum(t) = length(TrialStarts);

    plot(numTrials, 'Color', colors(t),'LineWidth',3);
    hold on;

    disp("Total number of trials:");
    disp(length(TrialStarts));
end

ttn = totTrialNum;
%legend("ACC1 ("+ttn(1)+")","ACC2 ("+ttn(2)+")","SAT1 ("+ttn(3)+")");
% legend("ACC1 ("+ttn(1)+")","ACC2 ("+ttn(2)+")","SAT1 ("+ttn(3)+")");
xlabel("Time (minutes)");
ylabel("Trial Count");
set(gca,'XTick',[0:1:length(numTrials)]+.5);
set(gca,'XTickLabel',0:timeBin/60:max(time)/60);
set(gca,'YGrid','on');
% title("Number of Trials ("+cageName+" - "+mouseName+")");
title("# of Initiated Trials");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
% Extracting Licking Frequency from the data
%%%%%%%%%%%%%%%%%%%%

GoLickFreq = {};
NoGoLickFreq = {};

for t = 1:length(data)
    temptest = data{t};
    
    time = table2array(temptest(:,1));
    time = time-min(time);
    water = table2array(temptest(:,2));
    lick = table2array(temptest(:,3));
    block = table2array(temptest(:,4));
    delay = table2array(temptest(:,5));
    
    % List of GO and NOGO trail start indices & times
    GoStarts = [];
    NoGoStarts = [];
    for i = 2:length(block)
        startTime = time(i);
        if block(i)==3 && block(i-1)~=3
            GoStarts = [GoStarts; i startTime];
        end
        if block(i)==9 && block(i-1)~=9
            NoGoStarts = [NoGoStarts; i startTime];
        end
    end
    
    % Now we will get the licking frequency within time blocks from the Go
    % trials:
    TrialGoData = {};
    
    timeBin = 300; %seconds
    numBins = ceil(max(time)/timeBin);
    ind = 1;
    for step = 1:numBins
        stepLickFreq = [];
        if ind<=length(GoStarts)
                currTime = GoStarts(ind,2);
        end
        while ind<=length(GoStarts) & currTime>=((step-1)*timeBin) & currTime<step*timeBin;
            % In this loop we go to the GO index in the raw data where the
            % trial starts, and from there get the licking frequency, and
            % then put it in an array that will be averaged later and go
            % into the licking frequency plot.
            dataIndex = GoStarts(ind,1);
            numDelay = ceil(delay(dataIndex)/100);
            newInd = dataIndex + numDelay + 8;
            currLickFreq = 1000 * (sum(lick(newInd:newInd+2))/2) / 300; %Hz
            stepLickFreq = [stepLickFreq currLickFreq];
            
            ind = ind + 1;
            if ind<length(GoStarts)
                currTime = GoStarts(ind,2);
            end
        end
        TrialGoData{step} = stepLickFreq;
    end
    
    GoLickFreq{t} = TrialGoData;
    
    % Now we will get the licking frequency within time blocks from the 
    % NOGo trials:
    TrialNoGoData = {};
    
    timeBin = 300; %seconds
    numBins = ceil(max(time)/timeBin);
    ind = 1;
    for step = 1:numBins
        stepLickFreq = [];
        if ind<=height(NoGoStarts)
                currTime = NoGoStarts(ind,2);
        end
        while ind<=height(NoGoStarts) & currTime>=((step-1)*timeBin) & currTime<step*timeBin;
            % In this loop we go to the NoGO index in the raw data where 
            % the trial starts, and from there get the licking frequency,
            % and then put it in an array that will be averaged later and 
            % go into the licking frequency plot.
            disp(".");
            dataIndex = NoGoStarts(ind,1);
            numDelay = ceil(delay(dataIndex)/100);
            newInd = dataIndex + numDelay + 8;
            currLickFreq = 1000 * (sum(lick(newInd:newInd+2))/2) / 300; %Hz
            stepLickFreq = [stepLickFreq currLickFreq];
            
            ind = ind + 1;
            if ind<=height(NoGoStarts)
                currTime = NoGoStarts(ind,2);
            end
        end
        TrialNoGoData{step} = stepLickFreq;
    end
    
    NoGoLickFreq{t} = TrialNoGoData;
    
end

%% Plotting licking freq
% Now that we have the licking frequencies from the Go and NoGo trials we
% can plot them

f=figure;
f.Position = [300, 800, 600, 700];
% figure
for t = 1:length(data)
    
    GoAvg = [];
    GoStd = [];
    for i = 1:length(GoLickFreq{t})
        GoAvg = [GoAvg mean(GoLickFreq{t}{i})];
        GoStd = [GoStd std(GoLickFreq{t}{i})];
    end
    
    NoGoAvg = [];
    NoGoStd = [];
    for i = 1:length(NoGoLickFreq{t})
        NoGoAvg = [NoGoAvg mean(NoGoLickFreq{t}{i})];
        NoGoStd = [NoGoStd std(NoGoLickFreq{t}{i})];
    end
    
    subplot(length(data),1,t);
    plot(GoAvg,'g-o','LineWidth',1.5,'MarkerSize',10,...
        'MarkerEdgeColor','g','MarkerFaceColor','g'); hold on;
    plot(NoGoAvg,'r-o','LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor','r');
    
    xlabel("Time (minutes)");
    ylabel("Licking Freq. (Hz)");
    set(gca,'XTick',[0:1:length(numTrials)]+.5);
    set(gca,'XTickLabel',0:timeBin/60:max(time)/60);
    xlim([0,14]);
    ylim([0,11]);
    title(days(t));
%     sgtitle("Licking Freq. ("+cageName+" - "+mouseName+")");
    sgtitle("Anticipatory Licking Frequency (Hz) Across 1 Hour of SAT");
    
end

%% Plotting average licking freq between Go and NoGo trials

f=figure;
f.Position = [300, 800,1000,200];
% figure
for t = 1:length(data)
    
    AllGoFreqs = [];
    for i = 1:length(GoLickFreq{1,t})
        AllGoFreqs = [AllGoFreqs cell2mat(GoLickFreq{1,t}(i))];
    end
    
    AllNoGoFreqs = [];
    for i = 1:length(NoGoLickFreq{1,t})
        AllNoGoFreqs = [AllNoGoFreqs cell2mat(NoGoLickFreq{1,t}(i))];
    end
    
    subplot(1,length(data),t);
    plot([1,2],[mean(AllGoFreqs), mean(AllNoGoFreqs)],'black','LineWidth',2); hold on;
    scatter(1,mean(AllGoFreqs),100,'Filled','g'); hold on;
    scatter(2,mean(AllNoGoFreqs),100,'Filled','r'); hold on;
    
     ylabel("Licking Freq. (Hz)");
     set(gca,'XTick',[1,2]);
     set(gca,'XTickLabel',["Stim", "Blank"]);
     xlim([.5,2.5]);
     ylim([0,10]);
     title(days(t));
%      sgtitle("Stim vs. Blank Anticipatory Licking Frequency ("+mouseName+", all)");
     sgtitle("Average Stim vs. Blank Anticipatory Licking Frequency (Hz)");

end

% Last 20% of trials
f=figure;
f.Position = [300, 800,1000,200];
for t = 1:length(data)
    
    AllGoFreqs = [];
    for i = 1:length(GoLickFreq{1,t})
        AllGoFreqs = [AllGoFreqs cell2mat(GoLickFreq{1,t}(i))];
    end
    percent20lenGo = floor(length(AllGoFreqs)*.2);
    GoFreq20 = AllGoFreqs(end-percent20lenGo : end);
    
    AllNoGoFreqs = [];
    for i = 1:length(NoGoLickFreq{1,t})
        AllNoGoFreqs = [AllNoGoFreqs cell2mat(NoGoLickFreq{1,t}(i))];
    end
    percent20lenNoGo = floor(length(AllNoGoFreqs)*.2);
    NoGoFreq20 = AllNoGoFreqs(end-percent20lenNoGo : end);
    
    subplot(1,length(data),t);
    plot([1,2],[mean(GoFreq20), mean(NoGoFreq20)],'black','LineWidth',2); hold on;
    scatter(1,mean(GoFreq20),100,'Filled','g'); hold on;
    scatter(2,mean(NoGoFreq20),100,'Filled','r'); hold on;
    
    ylabel("Licking Freq. (Hz)");
    set(gca,'XTick',[1,2]);
    set(gca,'XTickLabel',["Go", "NoGo"]);
    xlim([.5,2.5]);
    ylim([0,10]);
    title(days(t));
    sgtitle("Go vs. NoGo Licking Frequency ("+mouseName+", last 20% trials)"); 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% T-Test stuff
% %                   Acc1    Acc2    SAT1   SAT2    SAT3     SAT4
% GoAllAverages = [ 2.6061, 3.2197, 3.9394, 9.0646, 8.8546, 9.4889;... %M1
%                   0.8676, 3.6585, 5.7386, 9.6588, 9.8718, 9.6992;... %F2
%                   2.2917, 2.1587, 1.3262, 6.9928, 9.2385, 9.7535;... %M1
%                   3.2143, 2.6716, 2.1377, 3.9312, 6.3810, 9.0999;... %F4
%                   3.5606, 2.7387, 2.5556, 6.3333, 8.3213, 8.2263;... %M5
%                   1.9850, 1.1282, 1.2145, 2.5000, 2.9806, 5.5271];   %F5
% 
% Go20Averages = [ 2.3188, 4.1667, 5.4667, 8.4615, 6.5409, 8.1967;... %M1
%                  0.6250, 3.2353, 6.9444, 9.8077, 9.6226, 9.7619;... %F2
%                  2.8571, 3.7879, 3.6667, 8.7719, 9.3617, 9.7354;... %F3
%                  3.0556, 2.3810, 2.8070, 3.5088, 4.6512, 8.9583;... %F4
%                  3.3333, 2.7193, 3.4667, 5.1773, 7.9532, 7.2727;... %M5
%                  2.9825, 1.1111, 1.1111, 2.6984, 4.0171, 6.2500];   %F5
%              
% NoGoAllAverages = [0.7937, 2.4242, 4.3860, 7.4242, 6.3978, 5.0397;... %M1
%                    1.9444, 3.7778, 3.4815, 7.9333, 6.5641, 2.3810;... %F2
%                    2.5926, 3.4848, 0.9333, 6.1212, 7.1631, 7.5203;... %F3
%                    2.8571, 2.7778, 3.4667, 5.8333, 7.0175, 8.3077;... %F4
%                    3.1250, 2.2222, 1.8627, 5.3595, 7.1144, 6.6667;... %M5
%                    1.9608, 1.7460, 1.5238, 1.6667, 3.1624, 5.5556];   %F2
%                
% NoGo20Averages = [0.6667, 2.0000, 4.0741, 6.9048, 3.0769, 0.3704;... %M1
%                   1.1111, 4.7619, 3.6667, 8.1818, 2.8571, 1.1905;... %F2
%                   1.3333, 1.3333, 0.0000, 8.6111, 6.0000, 6.4706;... %F4
%                   3.3333, 2.3810, 3.3333, 4.5455, 6.6667, 6.9048;... %F4
%                   0.8333, 2.9167, 1.2500, 4.5455, 6.9048, 5.2083;... %M5
%                   0.0000, 2.9630, 2.9167, 2.3077, 2.9630, 5.6667];   %F2
%        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Getting the data to fill in the above chart, but only 1 mouse at a time
% GoAllAverages = [];
% Go20Averages = [];
% NoGoAllAverages = [];
% NoGo20Averages = [];
% 
% for d = 1:length(days);
%     
%     goPart = GoLickFreq{d};
%     allGoArray = [];
%     for i = 1:length(goPart)
%         allGoArray = [allGoArray, goPart{i}];
%     end
%     allGoAvg = mean(allGoArray);
%     numOf20PercentGo = round(length(allGoArray)*.2);
%     last20AvgGo = mean(allGoArray(end-numOf20PercentGo:end));
%     
%     nogoPart = NoGoLickFreq{d};
%     allNoGoArray = [];
%     for i = 1:length(nogoPart)
%         allNoGoArray = [allNoGoArray, nogoPart{i}];
%     end
%     allNoGoAvg = mean(allNoGoArray);
%     numOf20PercentNoGo = round(length(allNoGoArray)*.2);
%     last20AvgNoGo = mean(allNoGoArray(end-numOf20PercentNoGo:end));
%     
%     GoAllAverages = [GoAllAverages, allGoAvg];
%     Go20Averages = [Go20Averages, last20AvgGo];
%     NoGoAllAverages = [NoGoAllAverages, allNoGoAvg];
%     NoGo20Averages = [NoGo20Averages, last20AvgNoGo];
%     
% end
% GoAllAverages
% Go20Averages
% NoGoAllAverages
% NoGo20Averages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% pValuesAll = [];
% pValues20 = [];
% for day = 1:6
%     pAll = signrank(GoAllAverages(:,day),NoGoAllAverages(:,day));
% %     [pAll, h] = ttest(GoAllAverages(:,day),NoGoAllAverages(:,day));
%     pValuesAll(day) = pAll;
%     
%     p20 = signrank(Go20Averages(:,day),NoGo20Averages(:,day));
% %     [p20, h] = ttest(Go20Averages(:,day),NoGo20Averages(:,day));
%     pValues20(day) = p20;
% end
%               
% % Combined All Trials
% f=figure;
% f.Position = [300, 800,1000,200];
% for t = 1:length(data);
%     
%     subplot(1,length(data),t);
%     for mouse = 1:6
%         if mouse == 5 | mouse == 6
%             plot([1,2],[GoAllAverages(mouse,t), NoGoAllAverages(mouse,t)],'k:','LineWidth',2); hold on;
%         else
%             plot([1,2],[GoAllAverages(mouse,t), NoGoAllAverages(mouse,t)],'black','LineWidth',2); hold on;
%         end
%         scatter(1,GoAllAverages(mouse,t),100,'Filled','g'); hold on;
%         scatter(2,NoGoAllAverages(mouse,t),100,'Filled','r'); hold on;
%     end
%     
%     ylabel("Licking Freq. (Hz)");
%     set(gca,'XTick',[1,2]);
%     set(gca,'XTickLabel',["Go", "NoGo"]);
%     xlim([.5,2.5]);
%     ylim([0,10]);
%     title(days(t)+", p="+pValuesAll(t));
%     sgtitle("Go vs. NoGo Licking Frequency (All trials)");
%     
% end
% 
% % Combined 20% Last Trials
% f=figure;
% f.Position = [300, 800,1000,200];
% for t = 1:length(data);
%     
%     subplot(1,length(data),t);
%     for mouse = 1:6
%         if mouse == 5 | mouse == 6
%             plot([1,2],[Go20Averages(mouse,t), NoGo20Averages(mouse,t)],'k:','LineWidth',2); hold on;
%         else
%             plot([1,2],[Go20Averages(mouse,t), NoGo20Averages(mouse,t)],'black','LineWidth',2); hold on; 
%         end
%         scatter(1,Go20Averages(mouse,t),100,'Filled','g'); hold on;
%         scatter(2,NoGo20Averages(mouse,t),100,'Filled','r'); hold on;
%     end
%     
%     ylabel("Licking Freq. (Hz)");
%     set(gca,'XTick',[1,2]);
%     set(gca,'XTickLabel',["Go", "NoGo"]);
%     xlim([.5,2.5]);
%     ylim([0,10]);
%     title(days(t)+", p="+pValues20(t));
%     sgtitle("Go vs. NoGo Licking Frequency (last 20% trials)");
%     
% end
%            