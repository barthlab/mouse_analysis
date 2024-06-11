% export data to originlab to plot normalized curve
% initial created date: 11/15/2019
% modified date: 12/08/2019
% input requirements: Fall from suite2p / arduino xlsx file / arduino time
%                     point xlsx file
% note: change the saving place


ops.xoff(:,[3059:3061,6120:6122,9181:9183,12242:12244,15303:15305,18364:18366]) = []; % M8
ops.yoff(:,[3059:3061,6120:6122,9181:9183,12242:12244,15303:15305,18364:18366]) = []; % M8
F(:,[3059:3061,6120:6122,9181:9183,12242:12244,15303:15305,18364:18366]) = []; % M8
Fneu(:,[3059:3061,6120:6122,9181:9183,12242:12244,15303:15305,18364:18366]) = []; % M8
F(:,[3059:3061,6120:6122,9181:9183,12242:12244,15303:15305,18364:18366,21425:21427,24486:24488]) = []; % M9
%% find iscells from suite2p; calculate baseline upon local minima; one-day data per matrix; one-cell data per colume

imaging_day = 1; % imaging session

iscell_index = find(iscell(:,1)==1);

base_per_day = {};
for k = 1:length(iscell_index)
    data_per_day1 = reshape(F(iscell_index(k),:),[size(F,2)/imaging_day],imaging_day);        % change the number of 10-minute session
    data_per_day2 = reshape(Fneu(iscell_index(k),:),[size(Fneu,2)/imaging_day],imaging_day);
    data_per_day = data_per_day1 - 0.7*data_per_day2; 
%     data_per_day = reshape(F(iscell_index(k),:),[size(F,2)/30],30);
%     data_per_day = reshape(Fneu(iscell_index(k),:),[size(Fneu,2)/32],32);
    
    frame_num = size(data_per_day,1);

    base_per_day{end+1} = data_per_day;
end

% cell shift

cellshift_x = reshape(ops.xoff,[size(F,2)/imaging_day],imaging_day); % change the number of imaging days
cellshift_y = reshape(ops.yoff,[size(F,2)/imaging_day],imaging_day);
cellshift_vector = sqrt(cellshift_x.^2 + cellshift_y.^2);

% cellshift_vector_z3 = zscore(cellshift_vector);
cellshift_vector2 = {cellshift_vector};
cellshift_vector_z = repelem(cellshift_vector2,length(iscell_index));

%% align arduino time point with matrix
[~,sheet_name1]=xlsfinfo('Arduino.xlsx');
for k=1:numel(sheet_name1)
  [arduino_parameter{k},arduino{k}]=xlsread('Arduino.xlsx',sheet_name1{k});
end

[~,sheet_name2]=xlsfinfo('Arduino time point.xlsx')
for l=1:numel(sheet_name2)
  arduino_time{l}=xlsread('Arduino time point.xlsx',sheet_name2{l});
end


range_puff_all_page = {};
range_blank_all_page = {};
puff_merge = {};
blank_merge = {};
c_all = {};
pixel_puff_all_page = {};
pixel_blank_all_page = {};
puff_pixel_merge = {};
blank_pixel_merge = {};

for z = 1: length(iscell_index)
    for w = 1:k
    

        framerate_row = 1: length(base_per_day{z});
        framerate = framerate_row'*600000/length(base_per_day{z});
        intensity_base = base_per_day{z}(:,w);
        pixel = cellshift_vector_z{z}(:,w);

%     sampling_freq = [];
%     for i = 2:lengthclear(intensity_output(:,1))
%         sampling_freq(i) = intensity_output(i) - intensity_output(i-1);
%     end


        trial_time = arduino_time{w}(:,1);

        b=1;
        for i = 1:length(trial_time)
            for j = 1:length(framerate)-1
                a1 = abs(framerate(b)-trial_time(i));
                a2 = abs(framerate(j+1)-trial_time(i));
                if a1 > a2
                    b = j+1;
                end
            end
            c(i) = b;       
        end

        range = [];
        pixel_range = [];
        for xyz = 1:length(c)
            range_pre = [];
            pixel_pre = [];
            celem = c(xyz);
            for n = (celem-51):1:(celem+50)                                
                range_pre = [range_pre, intensity_base(n)];
                pixel_pre = [pixel_pre, pixel(n)];
            end
            range(xyz, :) = [range_pre];
            pixel_range(xyz, :) = [pixel_pre];
        end
        range_new = range';
        pixel_new = pixel_range';
%         if z == 1
%             if w == 3 || w == 5
%                 disp(range_new(1,1));
%             end
%         end
        

        puff = find(strcmp(arduino{w},'real'));
        blank = find(strcmp(arduino{w},'fake'));

        range_puff = [];
        pixel_puff = [];
        for v = 1:length(puff)
            range_puff = [range_puff, range_new(:,puff(v))];
            pixel_puff = [pixel_puff, pixel_new(:,puff(v))];
        end

        range_blank = [];
        pixel_blank = [];
        for q = 1:length(blank)
            range_blank = [range_blank, range_new(:,blank(q))];
            pixel_blank = [pixel_blank, pixel_new(:,blank(q))];
        end
        range_puff_all_page{end+1} = range_puff;
        pixel_puff_all_page{end+1} = pixel_puff;
        range_blank_all_page{end+1} = range_blank; 
        pixel_blank_all_page{end+1} = pixel_blank;
        
%         filename_puff = sprintf('PuffIntensity_%d',z)
%         xlswrite(fullfile('E:\Mo Zhu - 2019\S1 L23 awake imaging\191106 Emx-cre 2m\Data analysis\intensity output', filename_puff), range_puff, w, 'A1');
%     
%         filename_blank = sprintf('BlankIntensity_%d',z)
%         xlswrite(fullfile('E:\Mo Zhu - 2019\S1 L23 awake imaging\191106 Emx-cre 2m\Data analysis\intensity output', filename_blank), range_blank, w, 'A1');
        c_all{end+1} = c;

        
    end
     
end

all_frame = {};
for q = 1:length(c_all)/length(iscell_index)
    trial_frame = c_all{q};
    all_frame{end+1} = trial_frame;
end

for g = 1:length(range_puff_all_page)/2
            
    puff_page = [range_puff_all_page{g*2-1}';range_puff_all_page{g*2}']; 
    puff_page_pixel = [pixel_puff_all_page{g*2-1}';pixel_puff_all_page{g*2}']; 
            
    puff_merge{end+1} = puff_page';
    puff_pixel_merge{end+1} = puff_page_pixel';
            
end
        
for g = 1:length(range_blank_all_page)/2
            
    blank_page = [range_blank_all_page{g*2-1}';range_blank_all_page{g*2}'];  
    blank_page_pixel = [pixel_blank_all_page{g*2-1}';pixel_blank_all_page{g*2}']; 
            
    blank_merge{end+1} = blank_page';
    blank_pixel_merge{end+1} = blank_page_pixel';
            
end

%calculate the change in pixel shift

d_shift = [];
d_shift_all = {};
for num4 = 1:length(puff_pixel_merge)
    pixel_merge_mat = puff_pixel_merge{num4};
    for num5 = 1:size(pixel_merge_mat,2)
        d_shift_col = abs(diff(pixel_merge_mat(:,num5)));
        d_shift = [d_shift, d_shift_col];    
    end
    d_shift_all{end+1} = d_shift;
    d_shift = [];   
end

d_shift_b = [];
d_shift_all_b = {};
for num6 = 1:length(blank_pixel_merge)
    pixel_merge_mat_b = blank_pixel_merge{num6};
    for num7 = 1:size(pixel_merge_mat_b,2)
        d_shift_col_b = abs(diff(pixel_merge_mat_b(:,num7)));
        d_shift_b = [d_shift_b, d_shift_col_b];    
    end
    d_shift_all_b{end+1} = d_shift_b;
    d_shift_b = [];   
end


%% correlation between Ca signal and pixel shift
% corr_all is correlation of every trial; corr_by_cell is averaged
% correlation of every cell

% average movement during 1s response window
ave_shift_all  = [];
ave_shift_trial2 = [];
ave_shift_trial_all = {};

for num8 = 1:(length(d_shift_all)/length(iscell_index))
    mat_shift = d_shift_all{num8};
    ave_shift_trial = mean(mat_shift(52:60,:));
    ave_shift = mean(mean(mat_shift(52:60,:)));
    ave_shift_all = [ave_shift_all, ave_shift];
    
    ave_shift_trial2 = [ave_shift_trial2,ave_shift_trial];
    ave_shift_trial_all{end+1} = ave_shift_trial2;
    ave_shift_trial2 = [];
end


% correlation

corr_cell = [];
corr_p_cell = [];
corr_all = {};
corr_p_all = {};
corr_ave_all = [];
corr_p_ave_all = [];
for num1 = 1:length(puff_merge)
    puff_mat = puff_merge{num1};
    pixel_mat = d_shift_all{num1};
        for num2 = 1:size(puff_mat,2)
%             [corr2,p_val] = corrcoef(puff_mat(47:64,num2),pixel_mat(47:64,num2));
            [corr2,p_val] = corrcoef(puff_mat(1:101,num2),pixel_mat(1:101,num2));

            corr = corr2(2);
            corr_p = p_val(2);
            corr_cell = [corr_cell, corr];
            corr_p_cell = [corr_p_cell, corr_p];
            corr_p_ave = length(find(corr_p_cell<0.05))/length(corr_p_cell);
            corr_ave = nanmean(corr_cell);
        end
        corr_all{end+1} = corr_cell;
        corr_p_all{end+1} = corr_p_cell;
        corr_cell = [];
        corr_p_cell = [];
        
        corr_ave_all = [corr_ave_all, corr_ave];
        corr_p_ave_all = [corr_p_ave, corr_p_ave_all];
end

corr_by_cell = reshape(corr_ave_all, [length(iscell_index), length(corr_ave_all)/length(iscell_index)]);
corr_p_by_cell = reshape (corr_p_ave_all, [length(iscell_index), length(corr_ave_all)/length(iscell_index)]);

corr_sheet = {corr_by_cell,corr_p_by_cell};

for num3 = 1:length(corr_sheet)
    xlswrite('pixel shift corr P6_2.xlsx',corr_sheet{num3},num3)
end

%% Spontaneous activity for duplicated cells

event_num = [];
event_all = {};
for x = 1:length(base_per_day)  
    for e = 1:2:size(base_per_day{x},2)  %only for duplicated cells!!! change to 1:size(base_per_day{x},2) for training data!!!
        base_mat2 = base_per_day{x};
        base_mat = base_mat2(1:459,:);
        pks_cell = findpeaks(base_mat(:,e));
        event_cell = length(find(pks_cell>mean(base_mat(:,e))+2*std(base_mat(:,e))));
        event_num(end+1) = event_cell;
        event_normed = event_num/event_num(1);
    end
    event_all{end+1} = event_normed;
    event_num = [];
end
event_mat = reshape(cell2mat(event_all),[length(cell2mat(event_all))/length(event_all),length(event_all)])';

%% Spontaneous activity for training dataset

% Spon activity before stimulus train

event_num = [];
event_num3 = [];
event_all = {};
event_all3 = {};
for x = 1:length(base_per_day)  
    for e = 1:size(base_per_day{x},2)
        base_mat2 = base_per_day{x};
        base_mat = base_mat2(1:459,:);
        base_mat3 = base_mat2(end-459:end,:);
        pks_cell = findpeaks(base_mat(:,e));
        pks_cell3 = findpeaks(base_mat3(:,e));
        event_cell = length(find(pks_cell>mean(base_mat(:,e))+2*std(base_mat(:,e))));
        event_cell3 = length(find(pks_cell3>mean(base_mat3(:,e))+2*std(base_mat3(:,e))));
        event_num(end+1) = event_cell;
        event_num3(end+1) = event_cell3;
%         event_normed = event_num/event_num(1);
    end
    event_all{end+1} = event_num;
    event_all3{end+1} = event_num3;
    event_num = [];
    event_num3 = [];
end

event_mat = reshape(cell2mat(event_all),[length(cell2mat(event_all))/length(event_all),length(event_all)])';
event_mat3 = reshape(cell2mat(event_all3),[length(cell2mat(event_all3))/length(event_all3),length(event_all3)])';

event_mat4 = [];
for a = 1:size(event_mat,1)
    event_cell = reshape(event_mat(a,:),2,size(event_mat,2)/2);
    event_day = sum(event_cell);
    event_mat4 = [event_mat4;event_day];
end

event_mat5 = [];
for a = 1:size(event_mat,1)
    event_cell3 = reshape(event_mat3(a,:),2,size(event_mat3,2)/2);
    event_day3 = sum(event_cell3);
    event_mat5 = [event_mat5;event_day3];
end

% event_normed_all = [];
% for a = 1:size(event_mat,1)
%     event_cell = reshape(event_mat(a,:),2,size(event_mat,2)/2);
%     event_day = sum(event_cell);
%     event_normed = event_day/mean(event_day(2:6)); %change baseline days: ACC2-6
%     event_normed_all = [event_normed_all;event_normed];
%     event_normed = [];
% end
% event_normed_all3 = [];
% for a = 1:size(event_mat3,1)
%     event_cell3 = reshape(event_mat3(a,:),2,size(event_mat3,2)/2);
%     event_day3 = sum(event_cell3);
%     event_normed3 = event_day3/mean(event_day3(2:6)); %change baseline days: ACC2-6
%     event_normed_all3 = [event_normed_all3;event_normed3];
%     event_normed3 = [];
% end

% xlswrite('spontaneous activity pre_P6_1.xlsx',event_normed_all);
% xlswrite('spontaneous activity post_P6_1.xlsx',event_normed_all3);

xlswrite('spon activity pre_nonnormed P6_2.xlsx',event_mat4);
xlswrite('spon activity post_nonnormed P6_2.xlsx',event_mat5);

% Spon activity during stimulus train

event_freq_all = [];
for num8 = 1:length(puff_merge)
    puff_rmv = puff_merge{num8};
    puff_rmv(52:56,:)  = [];  % remove signals from 1s response time window
    puff_rmv2 = puff_rmv(:);
    puff_rmv_pks = findpeaks(puff_rmv2);
    puff_rmv_event = length(find(puff_rmv_pks>mean(puff_rmv2)+2*std(puff_rmv2)));
    event_freq = puff_rmv_event/(size(puff_rmv,2)*((102/5.11)-(5/5.11)));
    event_freq_all(end+1) = event_freq;
end

spon_during = reshape(event_freq_all,[imaging_day/2, length(iscell_index)])';
xlswrite('spon activity during_freq P6_3.xlsx',spon_during);


%% Identify responsive cell

trial_mat2 = [];
puff_merge2 = {};
prestim_all= {};
poststim_all = {};
response2 = [];
response_all = {};
prestim_sd_all = {};
pos2 = [];
position_all = {};
baseline_raw_all = {};
for s = 1: length(puff_merge)
        for d = 1:size(puff_merge{s},2)
            baseline_raw = mean(puff_merge{s}(47:51,:));
            trial_mat = puff_merge{s}(:,d);
            dF_F = (trial_mat - mean(trial_mat(47:51)))/mean(trial_mat(47:51));
            trial_mat2 = [trial_mat2,dF_F];
            prestim = trial_mat2(47:51,:);
            prestim_base = mean(prestim);
            prestim_sd = std(prestim);
            poststim = trial_mat2(52:56,:);
            [poststim_base,position] = max(poststim);

            for e = 1:length(prestim_base)
                if poststim_base(e) > prestim_base(e)+2*prestim_sd(e)
                    response = 1;
                    pos = position(e);
                else
                    response = 0;
                    pos = 0;
                end
            end
            response2 = [response2,response];
            pos2 = [pos2,pos];
        end
        puff_merge2{end+1} = trial_mat2;
        prestim_all{end+1} = prestim_base;
        prestim_sd_all{end+1} = prestim_sd;
        poststim_all{end+1} = poststim_base;
        response_all{end+1} = response2;
        position_all{end+1} = pos2;
        baseline_raw_all{end+1} = baseline_raw;
        trial_mat2 = [];
        prestim_base = [];
        poststim_base = [];
        response2 = [];   
        pos2 = [];
        baseline_raw2=[];
end

% test whether the raw baseline is smaller than 0.25
test3 = [];
for i = 1:length(baseline_raw_all)
    test1 = baseline_raw_all{i};
%     test2 = test1(:)';
    test3 = [test3, test1];
end
test3 = test3';

test4 = [];
for i = 1:length(response_all)
    test5 = response_all{i};
%     test6 = test5(:)';
    test4 = [test4, test5];
end
test4 = test4';

ind = find(test3>0 & test3<0.25);
trial_remove = test4(ind);

% find position of the value in cell array based on index
len_all = [];
for i = 1:330
    len = length(baseline_raw_all{i});
    len_all(end+1) = len;
end
sum(len_all)'

% %rearrange peak position
% position_all2 = [];
% for p = 11:16:length(position_all)  %change sequence base on #of training days and which day to extract!!!
%     pos3 = position_all{p};
%     position_all2 = [position_all2,pos3];
% end
% position_all3 = nonzeros(position_all2);
% % position_med = median(position_all3);

% find AVE peak position for each cell
pos_ave_all = [];
for i = 1:length(position_all)
    pos_cell = position_all{i};
    pos_ave = mean(nonzeros(pos_cell));
    pos_ave_all(end+1) = pos_ave;
end
pos_ave_rearrange = reshape(pos_ave_all, [length(pos_ave_all)/length(iscell_index), length(iscell_index)]);
% xlswrite('peak time bin delay P6_1.xlsx',pos_ave_rearrange);


perctrial_all = [];
for num = 1:length(response_all)
    percent_trial = sum(response_all{num})/length(response_all{num});
    perctrial_all = [perctrial_all,percent_trial];
end
perctrial_re = reshape(perctrial_all,[16,length(response_all)/16])'; %change # of imaging day!!!
xlswrite('responsive trials percent P6_2.xlsx',perctrial_re);

resp_trial_all = {};
for num2 = 1:length(puff_merge2)
    resp_mat=[];
    for num3 = 1:length(response_all{num2})
        if response_all{num2}(num3)==1   %change between 1 and 0 to decide responsive / suppressed trials
            resp_trial = puff_merge2{num2}(:,num3);
        elseif response_all{num2}(num3)==0
               resp_trial = NaN(size(puff_merge2{num2}(:,num3)));
        end
        resp_mat = [resp_mat,resp_trial];
    end
    resp_trial_all{end+1} = resp_mat;
end


puff_merge4 = {};
for p = 1:length(resp_trial_all)
    trial_mat5 = resp_trial_all{p};
    trial_mat6 = mean(trial_mat5,2,'omitnan');
    puff_merge4{end+1} = trial_mat6;
end

ind = length(puff_merge4)/length(iscell_index);
for q = 1:ind
    ACC1_all = [puff_merge4{1:ind:length(puff_merge4)}];
    ACC2_all = [puff_merge4{2:ind:length(puff_merge4)}];
    ACC3_all = [puff_merge4{3:ind:length(puff_merge4)}];
    ACC4_all = [puff_merge4{4:ind:length(puff_merge4)}];
    ACC5_all = [puff_merge4{5:ind:length(puff_merge4)}];
    ACC6_all = [puff_merge4{6:ind:length(puff_merge4)}];
    SAT1_all = [puff_merge4{7:ind:length(puff_merge4)}];
    SAT2_all = [puff_merge4{8:ind:length(puff_merge4)}];
    SAT3_all = [puff_merge4{9:ind:length(puff_merge4)}];
    SAT4_all = [puff_merge4{10:ind:length(puff_merge4)}];
    SAT5_all = [puff_merge4{11:ind:length(puff_merge4)}];
    SAT6_all = [puff_merge4{12:ind:length(puff_merge4)}];
    SAT7_all = [puff_merge4{13:ind:length(puff_merge4)}];
    SAT8_all = [puff_merge4{14:ind:length(puff_merge4)}];
    SAT9_all = [puff_merge4{15:ind:length(puff_merge4)}];
    SAT10_all = [puff_merge4{16:ind:length(puff_merge4)}];

end

puff_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all,SAT1_all,SAT2_all,SAT3_all,SAT4_all,SAT5_all,SAT6_all,SAT7_all,SAT8_all,SAT9_all,SAT10_all};
% puff_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all};

puff_sheet3 = {};
for w = 1:length(puff_sheet)
    mat_day = puff_sheet{w};
    mat_day2 = mat_day(47:64,:);
    mat_day3 = mat_day2';
    puff_sheet3{end+1} = mat_day3;
end

for x = 1:length(puff_sheet3)
    xlswrite('responsive trials only P6_2.xlsx',puff_sheet3{x},x)
end

% sum_all = [];
% for a = 1:length(response_all)
%     sum_cell = sum(response_all{a});
%     sum_all = [sum_all,sum_cell];
% end
% 
% sum_all2 = reshape(sum_all,16,a/16); %change # of day
% sigcell_num = [];
% for e = 1:size(sum_all2,1)
%     sigcell = length(find(sum_all2(e,:)));
%     sigcell_num = [sigcell_num,sigcell];
% end


%% Sparse representation - what fraction of cell on a given trial is responsive

trial_bycell = [];
trial_bycell_all = {};
for j = 1: length(response_all)/length(iscell_index)
    for i = 1:length(iscell_index)
        trial = response_all{(i-1)*(length(response_all)/length(iscell_index))+j};
        trial_bycell = [trial_bycell;trial];
    end
    trial_bycell_all{end+1} = trial_bycell;
    trial_bycell = [];
end

percent_all = [];
for k = 1:length(trial_bycell_all)
    percent_cell = sum(trial_bycell_all{k})/size(trial_bycell_all{k},1);
    percent_cell_ave = mean(percent_cell);
    percent_all(end+1) = percent_cell_ave;
end

%% dF/Fo

trial_mat2 = [];
puff_merge2 = {};
for s = 1: length(puff_merge)
        for d = 1:size(puff_merge{s},2)
            trial_mat = puff_merge{s}(:,d);
            dF_F = (trial_mat-mean(trial_mat(47:51)))/mean(trial_mat(47:51));
            trial_mat2 = [trial_mat2,dF_F];
        end
        puff_merge2{end+1} = trial_mat2;
        trial_mat2 = [];
end
    
puff_merge3 = {};
for p = 1:length(puff_merge2)
    trial_mat3 = puff_merge2{p};
    trial_mat4 = mean(trial_mat3,2);
    puff_merge3{end+1} = trial_mat4;
end

ind = length(puff_merge3)/length(iscell_index);
for q = 1:ind
    ACC1_all = [puff_merge3{1:ind:length(puff_merge3)}];
    ACC2_all = [puff_merge3{2:ind:length(puff_merge3)}];
    ACC3_all = [puff_merge3{3:ind:length(puff_merge3)}];
    ACC4_all = [puff_merge3{4:ind:length(puff_merge3)}];
    ACC5_all = [puff_merge3{5:ind:length(puff_merge3)}];
    ACC6_all = [puff_merge3{6:ind:length(puff_merge3)}];
    SAT1_all = [puff_merge3{7:ind:length(puff_merge3)}];
    SAT2_all = [puff_merge3{8:ind:length(puff_merge3)}];
    SAT3_all = [puff_merge3{9:ind:length(puff_merge3)}];
    SAT4_all = [puff_merge3{10:ind:length(puff_merge3)}];
    SAT5_all = [puff_merge3{11:ind:length(puff_merge3)}];
    SAT6_all = [puff_merge3{12:ind:length(puff_merge3)}];
    SAT7_all = [puff_merge3{13:ind:length(puff_merge3)}];
    SAT8_all = [puff_merge3{14:ind:length(puff_merge3)}];
    SAT9_all = [puff_merge3{15:ind:length(puff_merge3)}];
    SAT10_all = [puff_merge3{16:ind:length(puff_merge3)}];

end

puff_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all,SAT1_all,SAT2_all,SAT3_all,SAT4_all,SAT5_all,SAT6_all,SAT7_all,SAT8_all,SAT9_all,SAT10_all};
% puff_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all};

puff_sheet2 = {};
for w = 1:length(puff_sheet)
    mat_day = puff_sheet{w};
    mat_day2 = mat_day(47:64,:);
    mat_day3 = mat_day2';
    puff_sheet2{end+1} = mat_day3;
end

% for testing only ----------------------
% max_test_all = [];
% for test=1:16
%     mat_test = puff_sheet2{test};
% %     ave_test = mean(mat_test);
% %     max_test = max(ave_test(6:10));
%     max_test = max(mat_test(:,6:10),[],2);
%     max_test_all = [max_test_all,max_test];
% end
% ave_max_test = mean(max_test_all);
%-----------------------------------------
for x = 1:length(puff_sheet2)
    xlswrite('df.f_with Fneu correction P6_2.xlsx',puff_sheet2{x},x)
end


%% dF/Fo for blank trial

trial_mat2 = [];
blank_merge2 = {};
for s = 1: length(blank_merge)
        for d = 1:size(blank_merge{s},2)
            trial_mat = blank_merge{s}(:,d);
            dF_F = (trial_mat-mean(trial_mat(47:51)))/mean(trial_mat(47:51));
            trial_mat2 = [trial_mat2,dF_F];
        end
        blank_merge2{end+1} = trial_mat2;
        trial_mat2 = [];
end
    
blank_merge3 = {};
for p = 1:length(blank_merge2)
    trial_mat3 = blank_merge2{p};
    trial_mat4 = mean(trial_mat3,2);
    blank_merge3{end+1} = trial_mat4;
end

ind = length(blank_merge3)/length(iscell_index);
for q = 1:ind
    ACC1_all = [blank_merge3{1:ind:length(blank_merge3)}];
    ACC2_all = [blank_merge3{2:ind:length(blank_merge3)}];
    ACC3_all = [blank_merge3{3:ind:length(blank_merge3)}];
    ACC4_all = [blank_merge3{4:ind:length(blank_merge3)}];
    ACC5_all = [blank_merge3{5:ind:length(blank_merge3)}];
    ACC6_all = [blank_merge3{6:ind:length(blank_merge3)}];
    SAT1_all = [blank_merge3{7:ind:length(blank_merge3)}];
    SAT2_all = [blank_merge3{8:ind:length(blank_merge3)}];
    SAT3_all = [blank_merge3{9:ind:length(blank_merge3)}];
    SAT4_all = [blank_merge3{10:ind:length(blank_merge3)}];
    SAT5_all = [blank_merge3{11:ind:length(blank_merge3)}];
    SAT6_all = [blank_merge3{12:ind:length(blank_merge3)}];
    SAT7_all = [blank_merge3{13:ind:length(blank_merge3)}];
    SAT8_all = [blank_merge3{14:ind:length(blank_merge3)}];
    SAT9_all = [blank_merge3{15:ind:length(blank_merge3)}];
    SAT10_all = [blank_merge3{16:ind:length(blank_merge3)}];

end

blank_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all,SAT1_all,SAT2_all,SAT3_all,SAT4_all,SAT5_all,SAT6_all,SAT7_all,SAT8_all,SAT9_all,SAT10_all};
% blank_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all};

blank_sheet2 = {};
for w = 1:length(blank_sheet)
    mat_day = blank_sheet{w};
    mat_day2 = mat_day(47:64,:);
    mat_day3 = mat_day2';
    blank_sheet2{end+1} = mat_day3;
end

for x = 1:length(blank_sheet2)
    xlswrite('df.f_Fneu correction_blank_P6_2.xlsx',blank_sheet2{x},x)
end

%% dF/F for responsive blank trials

trial_mat2 = [];
blank_merge2 = {};
prestim_all= {};
poststim_all = {};
response2 = [];
response_all = {};
prestim_sd_all = {};
for s = 1: length(blank_merge)
        for d = 1:size(blank_merge{s},2)
            trial_mat = blank_merge{s}(:,d);
            dF_F = (trial_mat - mean(trial_mat(47:51)))/mean(trial_mat(47:51));
            trial_mat2 = [trial_mat2,dF_F];
            prestim = trial_mat2(47:51,:);
            prestim_base = mean(prestim);
            prestim_sd = std(prestim);
            poststim = trial_mat2(52:56,:);
            poststim_base = max(poststim);
            for e = 1:length(prestim_base)
                if poststim_base(e) > prestim_base(e)+2*prestim_sd(e)
                    response = 1;
                else
                    response = 0;
                end
            end
            response2 = [response2,response];
        end
        blank_merge2{end+1} = trial_mat2;
        prestim_all{end+1} = prestim_base;
        prestim_sd_all{end+1} = prestim_sd;
        poststim_all{end+1} = poststim_base;
        response_all{end+1} = response2;
        trial_mat2 = [];
        prestim_base = [];
        poststim_base = [];
        response2 = [];   
end

perctrial_all = [];
for num = 1:length(response_all)
    percent_trial = sum(response_all{num})/length(response_all{num});
    perctrial_all = [perctrial_all,percent_trial];
end
perctrial_re = reshape(perctrial_all,[16,length(response_all)/16])'; %change # of imaging day!!!
xlswrite('responsive trials_blank_percent_P6_2.xlsx',perctrial_re);

resp_trial_all = {};
for num2 = 1:length(blank_merge2)
    resp_mat=[];
    for num3 = 1:length(response_all{num2})
        if response_all{num2}(num3)==1   %change between 1 and 0 to decide responsive / suppressed trials
            resp_trial = blank_merge2{num2}(:,num3);
        elseif response_all{num2}(num3)==0
               resp_trial = NaN(size(blank_merge2{num2}(:,num3)));
        end
        resp_mat = [resp_mat,resp_trial];
    end
    resp_trial_all{end+1} = resp_mat;
end

blank_merge4 = {};
for p = 1:length(resp_trial_all)
    trial_mat5 = resp_trial_all{p};
    trial_mat6 = mean(trial_mat5,2,'omitnan');
    blank_merge4{end+1} = trial_mat6;
end

ind = length(blank_merge4)/length(iscell_index);
for q = 1:ind
    ACC1_all = [blank_merge4{1:ind:length(blank_merge4)}];
    ACC2_all = [blank_merge4{2:ind:length(blank_merge4)}];
    ACC3_all = [blank_merge4{3:ind:length(blank_merge4)}];
    ACC4_all = [blank_merge4{4:ind:length(blank_merge4)}];
    ACC5_all = [blank_merge4{5:ind:length(blank_merge4)}];
    ACC6_all = [blank_merge4{6:ind:length(blank_merge4)}];
    SAT1_all = [blank_merge4{7:ind:length(blank_merge4)}];
    SAT2_all = [blank_merge4{8:ind:length(blank_merge4)}];
    SAT3_all = [blank_merge4{9:ind:length(blank_merge4)}];
    SAT4_all = [blank_merge4{10:ind:length(blank_merge4)}];
    SAT5_all = [blank_merge4{11:ind:length(blank_merge4)}];
    SAT6_all = [blank_merge4{12:ind:length(blank_merge4)}];
    SAT7_all = [blank_merge4{13:ind:length(blank_merge4)}];
    SAT8_all = [blank_merge4{14:ind:length(blank_merge4)}];
    SAT9_all = [blank_merge4{15:ind:length(blank_merge4)}];
    SAT10_all = [blank_merge4{16:ind:length(blank_merge4)}];

end

blank_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all,SAT1_all,SAT2_all,SAT3_all,SAT4_all,SAT5_all,SAT6_all,SAT7_all,SAT8_all,SAT9_all,SAT10_all};
% blank_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all};

blank_sheet3 = {};
for w = 1:length(blank_sheet)
    mat_day = blank_sheet{w};
    mat_day2 = mat_day(47:64,:);
    mat_day3 = mat_day2';
    blank_sheet3{end+1} = mat_day3;
end

for x = 1:length(blank_sheet3)
    xlswrite('responsive trials_blank_P6_2.xlsx',blank_sheet3{x},x)
end


%% plot individual trials

frame_time = [600000/frame_num/2+600000/frame_num*-51:600000/frame_num:-600000/frame_num/2+600000/frame_num*51];
frame_time2 = frame_time(47:64);

%one trace per figure
for y = 1:length(puff_merge2)
    for z = 1:size(puff_merge2{y},2)
        cellbyday = puff_merge2{y};
        plot(frame_time2,cellbyday(47:64,z), 'Color','k','LineWidth',3)
        xlim([-900,2500]);
        set(gcf, 'Position',  [100, 100, 1200, 600])
        yline(0, 'r--', 'LineWidth', 3);
        saveas(gca,['Cell ',sprintf('%d', y),'_Trial ',sprintf('%d', z)],'jpg')
    end
end


% in same figure with different color
for y = 1:length(puff_merge2)
    c_ind = 0;
    for z = 1:size(puff_merge2{y},2)
        cellbyday = puff_merge2{y};
        hold on
        color = [(220-c_ind),(220-c_ind),(220-c_ind)]./255;
        plot(frame_time2,cellbyday(47:64,z), 'Color',color,'LineWidth',3)
        xlim([-900,2500]);
        c_ind = c_ind+8;
        set(gcf, 'Position',  [100, 100, 1200, 600])
        hold off 
    end
    saveas(gca,sprintf('%d.jpg', y))
    clf
end

for y = 1:length(puff_merge2)
    tiledlayout(1,size(puff_merge2{y},2));
    for z = 1:size(puff_merge2{y},2)
        ax = nexttile
        plot(frame_time2,cellbyday(47:64,z), 'Color','k')
        xlim([-900,2500]);
        linkaxes(ax,'xy')
        set(gcf, 'Position',  [100, 100, 4000, 100])
    end
end

%% area: 10s_before vs. 1s_post puff vs. 9s_after

frame_time = [600000/frame_num/2+600000/frame_num*-51:600000/frame_num:-600000/frame_num/2+600000/frame_num*51];
time_prop = (0-frame_time(1))/(frame_time(length(frame_time)/2+6)-frame_time(length(frame_time)/2+4));

area_under_puff_10s_all = [];
area_under_puff_1s_all = [];
area_under_puff_9s_all = [];
prop_ind_puff_all = [];
prop_ind_puff2 = [];
prop_ind_puff_alltrial = [];
meanMI_day_cell = [];
SEMI_day_cell = [];
ACC_mean_all = [];
ACC_all = [];
SAT1_all = [];
SAT2_all = [];
SAT3_all = [];
SAT4_all = [];
SAT5_all = [];
SAT6_all = [];
SAT7_all = [];
SAT8_all = [];
SAT9_all = [];
SAT10_all = [];
ACC1_all = [];
ACC2_all = [];
ACC3_all = [];
ACC4_all = [];
ACC5_all = [];
ACC6_all = [];
ACC5_6_all = [];

for order = -4:13
for s = 1: length(puff_merge)

        for d = 1:size(puff_merge{s},2)
            AreaMatrix_puff = puff_merge{s};
%             AreaFrame_puff_10s = AreaMatrix_puff(1:size(AreaMatrix_puff,1)/2,d);
            AreaFrame_puff_10s = AreaMatrix_puff(size(AreaMatrix_puff,1)/2-4:size(AreaMatrix_puff,1)/2,d);
            AreaFrame_puff_9s = AreaMatrix_puff(size(AreaMatrix_puff,1)/2+6:end,d);
%             AreaFrame_puff_1s = AreaMatrix_puff(size(AreaMatrix_puff,1)/2+5:size(AreaMatrix_puff,1)/2+6,d);
            AreaFrame_puff_1s = AreaMatrix_puff(size(AreaMatrix_puff,1)/2+order,d);

            positiveArea_puff_1s = sum(AreaFrame_puff_1s(AreaFrame_puff_1s>0));
            negativeArea_puff_1s = sum(AreaFrame_puff_1s(AreaFrame_puff_1s<0));
            positiveArea_puff_9s = sum(AreaFrame_puff_9s(AreaFrame_puff_9s>0));
            negativeArea_puff_9s = sum(AreaFrame_puff_9s(AreaFrame_puff_9s<0));
            positiveArea_puff_10s = sum(AreaFrame_puff_10s(AreaFrame_puff_10s>0));
            negativeArea_puff_10s = sum(AreaFrame_puff_10s(AreaFrame_puff_10s<0));

            area_under_puff_1s = positiveArea_puff_1s + negativeArea_puff_1s;
            area_under_puff_9s = (positiveArea_puff_9s + negativeArea_puff_9s)/9;
            area_under_puff_10s = (positiveArea_puff_10s + negativeArea_puff_10s)/(5/1);%51/2
            
            prop_ind_puff3 = (area_under_puff_1s-area_under_puff_10s)/area_under_puff_10s;
            prop_ind_puff2(end+1) = prop_ind_puff3;
            prop_ind_puff = mean(prop_ind_puff2);
            prop_ind_puff2 = [];
            
            %--------------------------------------------------------------------------------
           
        end

        area_under_puff_1s = [];
        area_under_puff_9s = [];
        area_under_puff_10s = [];
        
        prop_ind_puff_all(end+1) = prop_ind_puff;
%         meanMI_day_cell(end+1) = mean_MI;
%         SEMI_day_cell(end+1) = SE_MI;

end

% testmean_all = [];
% for testnum1 = 1:length(example_all)
% testmean = mean(example_all{testnum1});
% testmean_all(end+1) = testmean;
% end
% re_testmean_all = reshape(testmean_all,[16,45]);
% up_thres= mean(testmean_all)+0.3*(max(testmean_all)-min(testmean_all))
% down_thres= mean(testmean_all)-0.3*(max(testmean_all)-min(testmean_all));

% ttest_puff_per_cell = reshape(ttest_puff,[16,90]);
% ttest_puff_per_cell2 = reshape(ttest_puff2,[16,90]);
prop_ind_puff_rearrange = reshape(prop_ind_puff_all,[16,length(iscell_index)]);
% meanMI_day_cell_rearrage = reshape(meanMI_day_cell,[16,43]);
% SEMI_day_cell_rearrage = reshape(SEMI_day_cell,[16,43]);

% ACC16 = prop_ind_puff_rearrange(1:6,:);
% ACC = reshape(ACC16,[length(iscell_index)*6,1]);
% ACC_mean = mean(mean(ACC16));
% ACC_mean_all(end+1) = ACC_mean;
% ACC_all = [ACC_all,ACC];

ACC1_1 = prop_ind_puff_rearrange(1,:);
ACC1 = reshape(ACC1_1,[length(iscell_index),1]);
ACC1_all = [ACC1_all,ACC1];
ACC2_1 = prop_ind_puff_rearrange(2,:);
ACC2 = reshape(ACC2_1,[length(iscell_index),1]);
ACC2_all = [ACC2_all,ACC2];
ACC3_1 = prop_ind_puff_rearrange(3,:);
ACC3 = reshape(ACC3_1,[length(iscell_index),1]);
ACC3_all = [ACC3_all,ACC3];
ACC4_1 = prop_ind_puff_rearrange(4,:);
ACC4 = reshape(ACC4_1,[length(iscell_index),1]);
ACC4_all = [ACC4_all,ACC4];
ACC5_1 = prop_ind_puff_rearrange(5,:);
ACC5 = reshape(ACC5_1,[length(iscell_index),1]);
ACC5_all = [ACC5_all,ACC5];
ACC6_1 = prop_ind_puff_rearrange(6,:);
ACC6 = reshape(ACC6_1,[length(iscell_index),1]);
ACC6_all = [ACC6_all,ACC6];

% ACC5_6_1 = prop_ind_puff_rearrange(5:6,:);
% ACC5_6 = reshape(ACC5_6_1,[2*length(iscell_index),1]);
% ACC5_6_all = [ACC5_6_all,ACC5_6];

SAT1_1 = prop_ind_puff_rearrange(7,:);
SAT1 = reshape(SAT1_1,[length(iscell_index),1]);
SAT1_all = [SAT1_all,SAT1];
SAT2_1 = prop_ind_puff_rearrange(8,:);
SAT2 = reshape(SAT2_1,[length(iscell_index),1]);
SAT2_all = [SAT2_all,SAT2];
SAT3_1 = prop_ind_puff_rearrange(9,:);
SAT3 = reshape(SAT3_1,[length(iscell_index),1]);
SAT3_all = [SAT3_all,SAT3];
SAT4_1 = prop_ind_puff_rearrange(10,:);
SAT4 = reshape(SAT4_1,[length(iscell_index),1]);
SAT4_all = [SAT4_all,SAT4];
SAT5_1 = prop_ind_puff_rearrange(11,:);
SAT5 = reshape(SAT5_1,[length(iscell_index),1]);
SAT5_all = [SAT5_all,SAT5];
SAT6_1 = prop_ind_puff_rearrange(12,:);
SAT6 = reshape(SAT6_1,[length(iscell_index),1]);
SAT6_all = [SAT6_all,SAT6];
SAT7_1 = prop_ind_puff_rearrange(13,:);
SAT7 = reshape(SAT7_1,[length(iscell_index),1]);
SAT7_all = [SAT7_all,SAT7];
SAT8_1 = prop_ind_puff_rearrange(14,:);
SAT8 = reshape(SAT8_1,[length(iscell_index),1]);
SAT8_all = [SAT8_all,SAT8];
SAT9_1 = prop_ind_puff_rearrange(15,:);
SAT9 = reshape(SAT9_1,[length(iscell_index),1]);
SAT9_all = [SAT9_all,SAT9];
SAT10_1 = prop_ind_puff_rearrange(16,:);
SAT10 = reshape(SAT10_1,[length(iscell_index),1]);
SAT10_all = [SAT10_all,SAT10];

prop_ind_puff_all = [];
end

puff_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all,SAT1_all,SAT2_all,SAT3_all,SAT4_all,SAT5_all,SAT6_all,SAT7_all,SAT8_all,SAT9_all,SAT10_all};
% puff_sheet = {ACC1_all,ACC2_all,ACC3_all,ACC4_all,ACC5_all,ACC6_all,SAT1_all,SAT2_all,SAT3_all,SAT4_all,SAT5_all,SAT6_all};


for x = 1:length(puff_sheet)
    xlswrite('df.f_Fneu.xlsx',puff_sheet{x},x)
end


%% cell shift

cellshift_x = reshape(ops.xoff,[size(F,2)/32],32); % change the number of imaging days
cellshift_y = reshape(ops.yoff,[size(F,2)/32],32);
cellshift_vector = sqrt(cellshift_x.^2 + cellshift_y.^2);

cellshift_vector_z = zscore(cellshift_vector);

% SD = std(cellshift_vector);
% xlswrite('cell_shift P6_1.xlsx',SD);

%------------------
% data_spon1=base_per_day{1,6}(1:100,11);
% data_spon2=cellshift_vector(1:100,11);
% data1 = [data_spon1,data_spon2];
% data_trial1=base_per_day{1,6}(472:2570,11);
% data_trial2=cellshift_vector(472:2570,11);
% data2 = [data_trial1,data_trial2];