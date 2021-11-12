%% FP analysis

close all;
clear;
clc;
cd 'YOUR FOLDER';
FP_file=dir('*.csv');
[~,ind] = sort([FP_file.datenum]);
FP_file = FP_file(ind);
Fs = 100; % sampling frequency
dt = 1/Fs;

n_max = length(FP_file);
for n=1
    FP_filename = FP_file(n).name;
    dff_table{n} = readtable(FP_filename);
    dff_array{n} = table2array(dff_table{n});
    dff{n}(:,1)=(1:size(dff_array{n},1))*dt;
    dff{n}(:,2:5)=dff_array{n}(:,1:4);
    p = polyfit(dff{n}(:,3),dff{n}(:,2),1);
    yfit = polyval(p,dff{n}(:,3));
    dff{n}(:,6) = (dff{n}(:,2) - yfit)./yfit*100;

    % Find onset
    [TTL_row_number,~] = find(dff{n}(:,4)>=1);

    % Align the start of data to the TTL
    dff_cut{n} = dff{n}(TTL_row_number(1):(TTL_row_number(1)+700*Fs),:);
    dff_cut{n}(:,1) = dff_cut{n}(:,1)-dff_cut{n}(1,1);%-10;
    dff_cut{n}(:,6) = dff_cut{n}(:,6) - mean(dff_cut{n}(1:900,6));
    pooled_FP(:,n) = dff_cut{n}(:,6);
    
    figure(1);
    subplot(n_max,1,n);
    plot(dff_cut{n}(:,1),dff_cut{n}(:,6));
    ylabel('\DeltaF/F0','FontSize',12);
    title(FP_filename, 'Interpreter', 'none','FontSize',12)
    hold on;
end

%% Breathing analysis
for n=1
    Res_file=dir('*.mat');
    res_filename = Res_file(1).name;
    % load respiration file
    
    importfile(res_filename);
    BlockNo = n;
    %8 channels, 3 blocks. Res Voltage = ch 2, Rate = ch 5, Amp = ch 6, TTL =
    %ch 3
    res{n}(:,2) = data(datastart(2,BlockNo) : dataend(2,BlockNo)); % ch2 is Res Voltage
    res{n}(:,3) = data(datastart(3,BlockNo) : dataend(3,BlockNo)); % ch3 is TTL
    res{n}(:,4) = data(datastart(5,BlockNo) : dataend(5,BlockNo)); % ch5 is Res Rate
    res{n}(:,6) = data(datastart(6,BlockNo) : dataend(6,BlockNo)); % ch6 is amplitude

    res{n}(res{n}(:,3)<0.05,3) = 0;
    res{n}(res{n}(:,3)>1,3) = 1;
    res{n}(isnan(res{n}))=0;
    res{n}(:,5) = smoothdata(res{n}(:,4),'gaussian',50);
    res{n}(:,7) = smoothdata(res{n}(:,6),'gaussian',50);
    
    res{n}(:,1) = (1:size(res{n},1))*dt; % Time

    % Find onset
    [res_TTL_row_number,~] = find(res{n}(:,3)==1);

    % Align the start of data to the TTL
    res_cut{n} = res{n}((res_TTL_row_number(1)):(res_TTL_row_number(1)+700*Fs),:);
    res_cut{n}(:,1) = res_cut{n}(:,1)-res_cut{n}(1,1);%-10;
    pooled_res_rate(:,n) = res_cut{n}(:,5);
    pooled_res_amp(:,n) = res_cut{n}(:,7);
    
    figure(2);
    subplot(n_max,1,n);
    plot(res_cut{n}(:,1),res_cut{n}(:,5));
    ylabel('Res Rate (bpm)','FontSize',12);
    title(FP_filename, 'Interpreter', 'none','FontSize',12)
    hold on;
    
    figure('position', [100, 100, 2000, 800]);
    subplot(n_max,1,n)
    yyaxis left;
    plot(dff_cut{n}(:,1),dff_cut{n}(:,6));
    ylabel('\DeltaF/F0 (%)','FontSize',20);

    yyaxis right;
    plot(res_cut{n}(:,1),res_cut{n}(:,5),'LineWidth',0.7);
    xlabel('Time (s)', 'FontSize',20)
    ylabel('Res Rate (bpm)','FontSize',20);
    title(res_filename,'Interpreter','none','FontSize',20);
    legend('Calcium','Res','FontSize',14)
    axis tight;
    hold on;
end

pooled_FP = downsample(pooled_FP,10);
pooled_res_rate = downsample(pooled_res_rate,10);
pooled_res_amp = downsample(pooled_res_amp,10);

writematrix(pooled_FP,[erase(res_filename,'.mat') '_pooled_FP.xlsx']);
writematrix(pooled_res_rate,[erase(res_filename,'.mat') '_pooled_res_rate.xlsx']);
writematrix(pooled_res_amp,[erase(res_filename,'.mat') '_pooled_res_amp.xlsx']);
time_vec = (1:size(pooled_FP,1))'*dt*10;

%% X-corr Analysis
clear;
clc;
cd 'YOUR FOLDER'
CeA_o = readmatrix('YOUR FILE.xlsx');
PBC_o = readmatrix('YOUR FILE.xlsx');

CeA = zscore(CeA_o);
CeA = smoothdata(CeA,'gaussian',100);
PBC = zscore(PBC_o);
PBC = smoothdata(PBC,'gaussian',100);

% Data structure: 4 samples; 4 columns for each animal, column 1: time; column 2: dff; column 3: rate; column 4: amp. 

Fs = 10;
% FP-rate xcorr
for i = 1:size(CeA,2)/4
    [acor_fp_rate_CeA(:,i),lag_fp_rate(:,i)]=xcorr(CeA(:,(4*(i-1)+2)),CeA(:,(4*(i-1)+3)),'coeff'); % plot the cross correlation between normalized data
    lag_fp_rate_CeA(:,i) = lag_fp_rate(:,i)/Fs/60;
end
for i = 1:size(PBC,2)/4
    [acor_fp_rate_PBC(:,i),lag_fp_rate_PBC(:,i)]=xcorr(PBC(:,(4*(i-1)+2)),PBC(:,(4*(i-1)+3)),'coeff'); % plot the cross correlation between normalized data
    lag_fp_rate_PBC(:,i) = lag_fp_rate_PBC(:,i)/Fs/60;
end
% FP-amp xcorr
for i = 1:size(CeA,2)/4
    [acor_fp_amp_CeA(:,i),lag_fp_amp(:,i)]=xcorr(CeA(:,(4*(i-1)+2)),CeA(:,(4*(i-1)+4)),'coeff'); % plot the cross correlation between normalized data
    lag_fp_amp_CeA(:,i) = lag_fp_amp(:,i)/Fs/60;
end
for i = 1:size(PBC,2)/4
    [acor_fp_amp_PBC(:,i),lag_fp_amp_PBC(:,i)]=xcorr(PBC(:,(4*(i-1)+2)),PBC(:,(4*(i-1)+4)),'coeff'); % plot the cross correlation between normalized data
    lag_fp_amp_PBC(:,i) = lag_fp_amp_PBC(:,i)/Fs/60;
end

%% Mechanical - thermal PSTH
%% PLOTTING
cd 'YOUR FOLDER';

matrix_FP = [];
matrix_Res = [];
tmp = dir('*.xlsx');
Fs = 100; % sampling frequency
dt = 1/Fs;

close all;
for i = 1:4
    tmp_name = tmp(i).name;
    matrix_FP{i} = readmatrix(tmp_name,'Sheet','FP');
    matrix_FP_13s{i} = matrix_FP{i}(701:2000,:);
    matrix_FP_13s{i} = matrix_FP_13s{i}-mean(matrix_FP_13s{i}(1:250,:));
    matrix_Res{i} = readmatrix(tmp_name,'Sheet','Res');
    matrix_Res_13s{i} = matrix_Res{i}(701:2000,:);
    matrix_Res_13s{i} = matrix_Res_13s{i}-mean(matrix_Res_13s{i}(1:250,:));
end
time_vec = (1:size(matrix_FP_13s{1},1))*dt-2.5;

figure(1);
yyaxis left;
plot(time_vec,mean(matrix_FP_13s{1},2),'-b','LineWidth',1);
hold on;
plot(time_vec,mean(matrix_FP_13s{2},2),'--b','LineWidth',1);
xlim([-2.5, 10]);
xlabel(['Time from stim (s)']);
ylabel(['\DeltaF/F0 (%)']);
legend('25 degrees','55 degrees');
yyaxis right;
plot(time_vec,mean(matrix_Res_13s{1},2),'-r',time_vec,mean(matrix_Res_13s{2},2),'--r');
ylim([-10, 50]);
ylabel(['Respiration Rate (bpm)']);
xline(0,'--k');
legend('25 degrees','55 degrees');
saveas(gcf,'thermal_25_55.tif');

%% plot thermal w/ error bar
close all;
figure(1);
h1 = shadedErrorBar(time_vec,mean(matrix_FP_13s{1},2),std(matrix_FP_13s{1},0,2)/sqrt(size(matrix_FP_13s{1},2)/2),'lineprops','-g','transparent',1);
h2 = shadedErrorBar(time_vec,mean(matrix_FP_13s{2},2),std(matrix_FP_13s{2},0,2)/sqrt(size(matrix_FP_13s{2},2)/2),'lineprops','-m','transparent',1);

% yyaxis(h1,'right')
xlim([-2.5, 10]);
xlabel(['Time from stim (s)'],'Fontsize',18);
ylabel(['\DeltaF/F0 (%)'],'Fontsize',18);
xline(0,'--k');
legend('25','55');
saveas(gcf,'thermal_FP.tif');

figure(2);
shadedErrorBar(time_vec,mean(matrix_Res_13s{1},2),std(matrix_Res_13s{1},0,2)/sqrt(size(matrix_Res_13s{1},2/2)),'lineprops','--g','transparent',1);
shadedErrorBar(time_vec,mean(matrix_Res_13s{2},2),std(matrix_Res_13s{2},0,2)/sqrt(size(matrix_Res_13s{2},2/2)),'lineprops','--m','transparent',1);
xlim([-2.5, 10]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Arial','fontsize',12)
xlabel('Time from stim (s)','Fontsize',18);
ylim([-10,50]);
ylabel('\DeltaRespiration Rate (bpm)','Fontsize',18);
xline(0,'--k');
legend('25','55');
saveas(gcf,'thermal_Res.tif');


%% AMPLITUDE
cd 'YOUR FOLDER';
% ORDER: 25 55 0G 300G
matrix_Res_Amp = [];
tmp = dir('*.xlsx');
Fs = 100; % sampling frequency
dt = 1/Fs;

close all;
for i = 1:4
    tmp_name = tmp(i).name;
    
    matrix_Res_Amp{i} = readmatrix(tmp_name);
    matrix_Res_Amp_13s{i} = matrix_Res_Amp{i}(701:2000,:);
    matrix_Res_Amp_13s{i} = matrix_Res_Amp_13s{i}-mean(matrix_Res_Amp_13s{i}(1:250,:));
    matrix_Res_Amp_13s_downsample{i} = downsample(matrix_Res_Amp_13s{i},10);
    matrix_Res_Amp_13s_downsample{i} = matrix_Res_Amp_13s_downsample{i}.*1000;
end
time_vec = (1:size(matrix_Res_Amp_13s_downsample{1},1))*dt-2.5;

%% AMPLITUDE AUC 
thermal = readtable('YOUR FILE.xlsx');
thermal_array = table2array(thermal);


AUC_thermal_bl = trapz(thermal_array(thermal_array(:,1)<0,1),thermal_array(thermal_array(:,1)<0,2:end));
AUC_thermal_stim = trapz(thermal_array((thermal_array(:,1)<5)&(thermal_array(:,1)>2.5),1),thermal_array((thermal_array(:,1)<5)&(thermal_array(:,1)>2.5),2:end));
mech = readtable('YOUR FILE.xlsx');
mech_array = table2array(mech);

AUC_mech_bl = trapz(mech_array(mech_array(:,1)<0,1),mech_array(mech_array(:,1)<0,2:end));
AUC_mech_stim = trapz(mech_array((mech_array(:,1)<5)&(mech_array(:,1)>2.5),1),mech_array((mech_array(:,1)<5)&(mech_array(:,1)>2.5),2:end));
