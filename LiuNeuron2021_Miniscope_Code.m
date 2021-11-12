%% Miniscope Breathing Analysis
% Breathing: 100Hz -> downsample to 10Hz;
% calcium: 20Hz -> temporal downsampling by 2, into 10Hz, used C01-C43, 43
% cells

close all;
clear;
clc;
cd 'YOUR FOLDER';
Fs = 100; % original sampling frequency
dt = 1/Fs;

%% Load input, LabChart .mat file
Res_file=dir('*.mat');
res_filename = Res_file.name;
importfile(res_filename);
BlockNo = 2;
%8 channels, 3 blocks. Res Voltage = ch 2, Rate = ch 5, Amp = ch 6, TTL =
%ch 4
res(:,2) = data(datastart(2,BlockNo) : dataend(2,BlockNo)); % ch2 is Res Voltage
res(:,3) = data(datastart(4,BlockNo) : dataend(4,BlockNo)); % ch3 is TTL
res(:,4) = data(datastart(5,BlockNo) : dataend(5,BlockNo)); % ch5 is Res Rate
res(:,6) = data(datastart(6,BlockNo) : dataend(6,BlockNo)); % ch6 is amplitude
 
res(res(:,3)<0.05,3) = 0;
res(res(:,3)>1,3) = 1;
res(isnan(res))=0;
res(:,5) = smoothdata(res(:,4),'gaussian',50);
res(:,7) = smoothdata(res(:,6),'gaussian',50);
res(:,1) = (1:size(res,1))*dt; % Time
 
% Find TTL onset
[res_TTL_row_number,~] = find(res(:,3)==1);
 
% Align the start of data to the TTL
res_cut = res((res_TTL_row_number(1)):(res_TTL_row_number(end)),:);
res_cut(:,1) = res_cut(:,1)-res_cut(1,1);
pooled_res_rate(:,1) = res_cut(:,5);
pooled_res_amp(:,1) = res_cut(:,7);
 
% plot the breathing rate for the current file
figure(1);
plot(res_cut(:,1),res_cut(:,5));
ylabel('Res Rate (bpm)','FontSize',12);
title(res_filename, 'Interpreter', 'none','FontSize',12)
hold on;

%% save breathing data
pooled_res_rate = downsample(pooled_res_rate,10);
pooled_res_amp = downsample(pooled_res_amp,10);
pooled_all = horzcat(pooled_res_rate(1:end,:),pooled_res_amp(1:end,:));

writematrix(pooled_all,[erase(res_filename,'.mat') '_all.xlsx']);

%% dff calculation
temp = readmatrix('YOUR MATRIX.xlsx');
sc_calcium = temp(:,4:end);
zscored_dff = zscore(sc_calcium);
time_vec = (1:size(zscored_dff,1))'*dt*10;
plot_xlim = [0,600];

% plot the calcium traces only
for i = 1:size(sc_calcium,2)
    figure(3);
    set(gcf, 'position', [100, 100, 2000, 800]);
    subplot(5, 9, i);
    plot(time_vec,zscored_dff(:,i));
    title(['Cell #',num2str((i))]);
    xlim(plot_xlim);
    hold on;
end
saveas(gcf,'zscored_dff.tif');

%% get breathing data
breathing = temp(:,1:2);

% plot the breathing traces only
for u = 1:2
    figure(4);
    set(gcf, 'position', [100, 500, 2000, 300]);
    subplot(2,1,u);
    plot(time_vec,breathing(:,u));
    xlim(plot_xlim);
    hold on;
end
saveas(gcf,'breathing rate_amp.tif');

%% if there are few cells, plot together
for i = 1:size(sc_calcium,2)
    figure(3);
    set(gcf, 'position', [100, 100, 2000, 800]);
    subplot(size(sc_calcium,2)+2, 1, i);
    plot(time_vec,zscored_dff(:,i));
    title(['Cell #',num2str((i))]);
    xlim(plot_xlim);
    hold on;
end
saveas(gcf,'zscored_dff.tif');

%% Miniscope Pinching PSTH Analysis
%% Get pinching data
close all;
clear;
clc;
cd 'YOUR FOLDER';
Fs = 10; % original sampling frequency
dt = 1/Fs;
n_cell = 43; % used for plotting in the later session
animal = 'M7';

%% FIRST PINCH 
% 1. load files
pinching1 = readmatrix('YOUR MATRIX.xlsx');
close all;
pinching1_zscored = pinching1;
pinching1_zscored(:,2:end) = zscore(pinching1(:,2:end)); % column 1 is the time

% 2. cut based on timestamps
pinchtime1 = [
154.661
271.586
];

dff_zscore_pinch_1 = [];
raw = pinching1_zscored;
pooled_dff = [];
Fs = 10;

for i = 1:length(pinchtime1)
    [~, index] = min(abs(raw(:,1)-pinchtime1(i)));
    dff_zscore_pinch_1{i} = raw((index-5*Fs):(index+10*Fs),:);
    dff_zscore_pinch_1{i}(:,2:end) = dff_zscore_pinch_1{i}(:,2:end)-mean(dff_zscore_pinch_1{i}(1:40,2:end),1); % VERY IMPORTANT NORMALIZATION!!!!!!!
    dff_zscore_pinch_1_notime{i} = dff_zscore_pinch_1{i}(:,2:end);
end

% 3. plot a heatmap and identify the cells
for q = 1:length(pinchtime1)
    matrix_for_heatmap = transpose(dff_zscore_pinch_1{q}(:,2:end));
    [~,idx] = sort(max(matrix_for_heatmap(:,60:120),[],2),'descend');
    matrix_for_heatmap_sort = matrix_for_heatmap(idx,:);
    figure('Position', [10 10 2000 1000]);
    heatmap=imagesc(matrix_for_heatmap_sort,[-1 2]);
    colormap parula;
    xticks(-51:50:99);
    xticklabels({[-10 -5 0 5]});
    yticks(1:1:n_cell);
    yticklabels(num2str(idx));
    ytickangle(0);
    colorbar;
    saveas(gcf,[animal '_SessI_Pinch#',num2str((q)),'_heatmap.tif']);
end

%% SECOND PINCH 
% 1. load files
pinching2 = readmatrix('YOUR MATRIX.xlsx');
%% 
close all;
pinching2_zscored = pinching2;
pinching2_zscored(:,2:end) = zscore(pinching2(:,2:end));

% 2. cut based on timestamps
pinchtime2 = [
38.91
149.969
];

dff_zscore_pinch_2 = [];
raw_2 = pinching2_zscored;
pooled_dff_2 = [];

for i = 1:length(pinchtime2)
    [~, index] = min(abs(raw_2(:,1)-pinchtime2(i)));
    dff_zscore_pinch_2{i} = raw_2((index-5*Fs):(index+10*Fs),:);
    dff_zscore_pinch_2{i}(:,2:end) = dff_zscore_pinch_2{i}(:,2:end)-mean(dff_zscore_pinch_2{i}(1:40,2:end),1); % VERY IMPORTANT NORMALIZATION!!!!!!!
    dff_zscore_pinch_2_notime{i} = dff_zscore_pinch_2{i}(:,2:end);

end

% 3. plot a heatmap and identify the cells
for q = 1:length(pinchtime2)
    matrix_for_heatmap = transpose(dff_zscore_pinch_2{q}(:,2:end));
    [~,idx] = sort(max(matrix_for_heatmap(:,60:120),[],2),'descend');
    matrix_for_heatmap_sort = matrix_for_heatmap(idx,:);
    figure('Position', [10 10 2000 1000]);
    heatmap=imagesc(matrix_for_heatmap_sort,[-1 4]);
    colormap parula;
    xticks(-51:50:99);
    xticklabels({[-10 -5 0 5]});
    yticks(1:1:n_cell);
    yticklabels(num2str(idx));
    ytickangle(0);
    colorbar;
    saveas(gcf,[animal '_SessII_Pinch#',num2str((q)),'_heatmap.tif']);
end

%% THIRD PINCH 
% 1. load files
pinching3 = readmatrix('YOUR MATRIX.xlsx');
%% 
close all;
pinching3_zscored = pinching3;
pinching3_zscored(:,2:end) = zscore(pinching3(:,2:end));
 
% 2. cut based on timestamps
pinchtime3 = [
101.771
175.381
];
 
dff_zscore_pinch_3 = [];
raw_3 = pinching3_zscored;
pooled_dff_3 = [];
 
for i = 1:length(pinchtime3)
    [~, index] = min(abs(raw_3(:,1)-pinchtime3(i)));
    dff_zscore_pinch_3{i} = raw_3((index-5*Fs):(index+10*Fs),:);
    dff_zscore_pinch_3{i}(:,2:end) = dff_zscore_pinch_3{i}(:,2:end)-mean(dff_zscore_pinch_3{i}(1:40,2:end),1); % VERY IMPORTANT NORMALIZATION!!!!!!!
    dff_zscore_pinch_3_notime{i} = dff_zscore_pinch_3{i}(:,2:end);
 
end
 
% 3. plot a heatmap and identify the cells
for q = 1:length(pinchtime3)
    matrix_for_heatmap = transpose(dff_zscore_pinch_3{q}(:,2:end));
    [~,idx] = sort(max(matrix_for_heatmap(:,60:120),[],2),'descend');
    matrix_for_heatmap_sort = matrix_for_heatmap(idx,:);
    figure('Position', [10 10 2000 1000]);
    heatmap=imagesc(matrix_for_heatmap_sort,[-1 4]);
    colormap parula;
    xticks(-51:50:99);
    xticklabels({[-10 -5 0 5]});
    yticks(1:1:n_cell);
    yticklabels(num2str(idx));
    ytickangle(0);
    colorbar;
    saveas(gcf,[animal '_SessIII_Pinch#',num2str((q)),'_heatmap.tif']);
end

%% combine pinching
single_cell_activity_pinch_concat_matrix = [];
dff_zscore_pinch_concat = horzcat(dff_zscore_pinch_1_notime,dff_zscore_pinch_2_notime,dff_zscore_pinch_3_notime);
for c = 1:size(dff_zscore_pinch_concat{1},2) % for all cells
    for m = 1:size(dff_zscore_pinch_concat,2) % for all episodes
        single_cell_activity_pinch{c}(:,m) = dff_zscore_pinch_concat{m}(:,c);
    end
    single_cell_activity_pinch_concat_matrix = horzcat(single_cell_activity_pinch_concat_matrix,single_cell_activity_pinch{c});
end
save([animal '_Pinch_singlecell_average.mat'],'single_cell_activity_pinch');

%% Plot individual cells
close all;
n_episode_pinch = length(pinchtime1) + length(pinchtime2) + length(pinchtime3);

figure(1);
set(gcf, 'position', [100, 100, 1000, 800]);
for p = 1:n_cell
    subplot(5,9,p)
    hold on
    temp = single_cell_activity_pinch{p};
    temp_1 = nanmean(temp,2);
    temp_2 = nanstd(temp,0,2)./sqrt(sum(~isnan(temp(1,:))));
    h = area(-5:0.1:10,[(temp_1-temp_2),(2*temp_2)],-2.5);
    set(h(1),'Edgecolor','none','FaceColor','none');
    set(h(2),'Edgecolor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
    plot(-5:0.1:10,mean(single_cell_activity_pinch{p},2),'k');
    title(['Cell#',num2str((p-1)),' Ensure Ave']);
    ylim([-2.5 4])
end
saveas(gcf,[animal '_Pinch_averaged trace.tif']);

%% Compare maximum activity before and after stimulus
% comparing 0-10 sec data before and after permutation
before_stim_max = max(single_cell_activity_pinch_concat_matrix(1:30,:),[],1);
after_stim_max = max(single_cell_activity_pinch_concat_matrix(50:150,:),[],1);

% t-test for each perm/umperm pair (6 consecutive numbers)
for p = 1:length(before_stim_max)/n_episode_pinch
    [h_beforeafter_max(p), p_beforeafter_max(p)] = ttest(before_stim_max(((p-1)*n_episode_pinch+1):p*n_episode_pinch),after_stim_max(((p-1)*n_episode_pinch+1):p*n_episode_pinch));
end
significant_cellno_beforeafter_max = sum(h_beforeafter_max);
 
% the cells that have significant changes upon pinch uptake is
cell_number_with_significance_beforeafter_max = find(h_beforeafter_max) -1 ;

% comparing 0-10 sec data before and after permutation
before_stim_min = min(single_cell_activity_pinch_concat_matrix(1:30,:),[],1);
after_stim_min = min(single_cell_activity_pinch_concat_matrix(50:150,:),[],1);

% t-test for each perm/umperm pair (6 consecutive numbers)
for p = 1:length(before_stim_min)/n_episode_pinch
    [h_beforeafter_min(p), p_beforeafter_min(p)] = ttest(before_stim_min(((p-1)*n_episode_pinch+1):p*n_episode_pinch),after_stim_min(((p-1)*n_episode_pinch+1):p*n_episode_pinch));
end
significant_cellno_beforeafter_min = sum(h_beforeafter_min);
 
% the cells that have significant changes upon pinch uptake is
cell_number_with_significance_beforeafter_min = find(h_beforeafter_min) -1 ;
cell_number_with_significance_beforeafter_max_and_min = unique(horzcat(cell_number_with_significance_beforeafter_min,cell_number_with_significance_beforeafter_max));
cell_number_with_significance = {cell_number_with_significance_beforeafter_max_and_min cell_number_with_significance_beforeafter_max cell_number_with_significance_beforeafter_min};
save([animal '_cells that respond to pinching.mat'],'cell_number_with_significance');


%% Plot both breathing and pinching heatmaps
%% get breathing data
input = readmatrix('YOUR MATRIX.csv');
% close all;
input_zscored = input;
input_zscored(:,4:end) = zscore(input_zscored(:,4:end)); % column 1 is the time

breathing = input_zscored(:,3)';
breathing = breathing * 1.2-0.2;
single_cell_breathing = input_zscored(:,4:end);

%% breathing heatmap
single_cell_breathing_for_heatmap = transpose(single_cell_breathing);

time1 = 165;
time2 = time1+15;

breathing_cut = breathing(:,time1*Fs:time2*Fs);
matrix_breathing = single_cell_breathing_for_heatmap(:,time1*Fs:time2*Fs);
matrix_for_heatmap_breathing_cut = vertcat(breathing_cut,matrix_breathing);

figure('Position', [200 200 2000 100]);
plot(breathing_cut);

[~,idx] = sort(mean(matrix_breathing,2),'descend');
matrix_breathing_sort = matrix_breathing(idx,:);
matrix_for_heatmap_breathing_sort_cut = vertcat(breathing_cut,matrix_breathing_sort);

% plot cut breathing heatmap
figure('Position', [10 10 2000 1000]);
heatmap=imagesc(matrix_for_heatmap_breathing_sort_cut,[-1.5 3]);
colormap parula;
yticks(1:1:n_cell);
yticklabels(idx);
colorbar;

%% sort the breathing heatmap
[~,idx0] = sort(mean(matrix_breathing,2),'descend');
idx1 = idx0(1:31); % discard the rest
matrix2 = matrix_breathing(idx1,:);

[~,idx2] = sort(mean(matrix2(:,:),2),'descend');
tmp1 = matrix2(idx2,:);
cut_out_matrix_1 = tmp1(1:2,:);
analyze = tmp1(3:end,:);

[~,idx3] = sort(mean(analyze(:,1:40),2),'descend');
tmp2 = analyze(idx3,:);
cut_out_matrix_2 = tmp2(1:8,:);
analyze_2 = tmp2(10:end,:);

[~,idx4] = sort(mean(analyze_2(:,50:100),2),'descend');
cut_out_matrix_3 = analyze_2(idx4,:);

matrix_plot = vertcat(cut_out_matrix_1,cut_out_matrix_2,cut_out_matrix_3);

%% get row order for the breathing heatmap, save the order in "Loc"
for i = 1:size(matrix_plot,1)
 [~,Loc(i)] = ismember(matrix_plot(i,:), matrix_breathing, 'rows');
end
    
%% re-plot pinching heatmap
matrix_for_heatmap_breathing_responsive_only = matrix_for_heatmap(Loc,:);
figure('Position', [10 10 2000 1000]);
heatmap_pinching_breathing_responsive_only=imagesc(matrix_for_heatmap_breathing_responsive_only,[-1 2]);
colormap parula;
axis off;
colorbar;
saveas(gcf,[animal '_pinch_heatmap.tif']);

%% re-plot breathing heatmap
breathing_cut = breathing_cut -0.5;
matrix_for_heatmap_breathing_responsive_only = vertcat(breathing_cut,matrix_plot);
figure('Position', [10 10 2000 1000]);
heatmap_breathing_responsive_only=imagesc(matrix_for_heatmap_breathing_responsive_only,[-1 2]);
colormap parula;
axis off;
colorbar;
saveas(gcf,[animal '_breathing_heatmap.tif']);
