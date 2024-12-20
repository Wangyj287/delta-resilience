
%% subfigure d e

clc
clear all
close all

data = readtable('result_residual.xlsx');
Subsystem1 = {'MK', 'GZ', 'ZY', 'DLS'};
Subsystem2 = {'SS', 'SD', 'SSJ', 'BSW', 'HM'};
Subsystem3 = {'HP', 'DS', 'SSK', 'SSW', 'NS'};
Subsystem4 = {'HS', 'BJ', 'HJ', 'SJK', 'GC'};
Subsystem5 = {'ZD', 'LS', 'FBC', 'NH', 'RQ', 'MA'};
color_all = [0.8392, 0.1529, 0.1569;...
             0.585098039215686,0.570588235294118,0.774509803921569;...
             0.994117647058824,0.405882352941177,0.094117647058824;...
             0.0156862745098039,0.282352941176471,0.419607843137255;...
             0.305882352941177,0.541176470588235,0.317647058823529];
variables = {'lamdvar', 'lamdac1', 'var', 'ac1'};
variable_names = {'\lambda_{Var}', '\lambda_{AR1}', 'Variance', 'AR1'};
offsets = [0.5, 0.5, 0.5,0.5, 0.5];
for varIdx = 1:length(variables)
    
    all_values = data.(variables{varIdx}); 
    normalized_all_values = 50 + 200 * (all_values - min(all_values)) / (max(all_values) - min(all_values)); 
    data.([variables{varIdx} '_normalized']) = normalized_all_values;
    figure;
    hold on;
    set(gcf, "Position", [100 100 800 200]);
    plotData(data, Subsystem1, variables{varIdx}, [variables{varIdx} '_normalized'], color_all(1,:), 'Subsystem 1', 'Column', offsets(1), 0.5);
    plotData(data, Subsystem2, variables{varIdx}, [variables{varIdx} '_normalized'], color_all(2,:), 'Subsystem 2', 'Column', offsets(2), 0.3);
    plotData(data, Subsystem3, variables{varIdx}, [variables{varIdx} '_normalized'], color_all(3,:), 'Subsystem 3', 'Column', offsets(3), 0.3);
    plotData(data, Subsystem4, variables{varIdx}, [variables{varIdx} '_normalized'], color_all(4,:), 'Subsystem 4', 'Column', offsets(4), 0.3);
    plotData(data, Subsystem5, variables{varIdx}, [variables{varIdx} '_normalized'], color_all(5,:), 'Subsystem 5', 'Column', offsets(5), 0.5);
  
    if strcmp(variable_names{varIdx}, 'AR1') || strcmp(variable_names{varIdx}, 'Variance')
        legend('Subsystem 1', 'Subsystem 2', 'Subsystem 3', 'Subsystem 4', 'Subsystem 5', 'Location', 'southeast', 'FontSize', 12);
        legend('boxoff')
    else
        legend('off'); 
    end

    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman'); 
    set(gca, 'linewidth', 2); 
    xlabel('Upstream to downstream', 'FontSize', 14, 'FontName', 'Times New Roman'); 
    ylabel(variable_names{varIdx}, 'FontSize', 14, 'FontName', 'Times New Roman');

 
    set(gca, 'XTick', []);  
    
    ax = gca;
    xlim_data = ax.XLim;
    ylim_data = ax.YLim;
    
    annotation('arrow', [0.13, 0.91], [0.14, 0.14], 'LineWidth',1);

    grid off; 
    
    hold off;
    filename = sprintf('figure_%d.tif', varIdx); 
    print(gcf, filename, '-dtiff', '-r300'); 
end

function plotData(data, subsystem, variable, normalized_variable, color, label, station_column, offset, alpha_value)
   
    indices = ismember(data.(station_column), subsystem);
    ordered_data = data(indices, :);
    [~, order] = ismember(subsystem, ordered_data.(station_column)); 
    ordered_data = ordered_data(order, :); 
    yValues = ordered_data.(variable);
    marker_sizes = ordered_data.(normalized_variable);
    xValues = 1:length(yValues);
    scatter(xValues + offset, yValues, marker_sizes, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
        'LineWidth', 1.5, 'MarkerFaceAlpha', alpha_value, 'MarkerEdgeAlpha', alpha_value);
end


%% subfigure f
figure
[data_all,txt]=xlsread('Residual_trend分析.xlsx',4);
[data_all2,txt]=xlsread('Residual_trend分析-var.xlsx',4);
grouporder={'AR1';'\lambda_A_R_1'};
t=data_all(3:end,[1:2]);
data_all1=data_all(:,3:end);
uu1=data_all1(3:end,[1:12 14:25 27 ]);
data_all22=data_all2(:,3:end);
uu2=data_all22(3:end,[1:12 14:25 27 ]);
num_data = length(data_all1(:,1));
std_data = 0.05;
years = t(:,1) + (t(:,2)-1)/12;
average_data1 = mean(uu1, 2);
average_data2 = mean(uu2, 2);
window_size = 12;
average_data1 = movmean(average_data1, window_size);
average_data2 = movmean(average_data2, window_size);
yyaxis left 
ax = gca;
ax.YColor = 'k'; 
color_mean = color_all(4,:);
color_between = color_all(5,:);
aa1=plot(years, average_data1, 'k-','linewidth',2);
hold on
uu1_up =average_data1 + std(average_data1);
uu1_down =average_data1 - std(average_data1);
for i = 1:length(t)-1
    x = [years(i),years(i+1),years(i+1),years(i)];
    y = [uu1_up(i),uu1_up(i+1),uu1_down(i+1),uu1_down(i)];
    h1 = fill(x,y,'m');
    set(h1,'Facecolor',color_mean,'FaceAlpha',0.3,'EdgeColor','none');
end
hold on;
coeffs1 = polyfit(years, average_data1, 1); 
trend1 = polyval(coeffs1, years); 
trend_line1 = plot(years, trend1, 'r--', 'LineWidth', 2); 
hold on;
ylabel('Resilience') 
ylim([-0.2,1])
yyaxis right
ax.YColor = 'b';  
num_years = 50;
annual_mean_data1 = zeros(num_years, 1);
for i = 1:num_years
    year_indices = (years >= (1966 + i - 1)) & (years < (1966 + i));
    annual_mean_data1(i) = mean(average_data1(year_indices));
end
log_annual_mean_data1 = log(annual_mean_data1);
log_annual_mean_data1=real(log_annual_mean_data1);
t1=1966:2015
aa2 = plot(t1', log_annual_mean_data1, 'b-', 'LineWidth', 2);
hold on;
std_log_annual_mean_data1 =std(log_annual_mean_data1);
uu2_up =log_annual_mean_data1 + std_log_annual_mean_data1;
uu2_down =log_annual_mean_data1 - std_log_annual_mean_data1;
for i = 1:length(t1)-1
    x = [t1(i),t1(i+1),t1(i+1),t1(i)];
    y = [uu2_up(i),uu2_up(i+1),uu2_down(i+1),uu2_down(i)];
    h2 = fill(x,y,'m');
    set(h2,'Facecolor',color_between,'FaceAlpha',0.3,'EdgeColor','none');
end
hold on;
ylabel('Recovery rate \lambda') 
ylim([ -8,0])
coeffs2 = polyfit(t1', log_annual_mean_data1, 1); 
trend2 = polyval(coeffs2, t1'); 
trend_line2 = plot(t1', trend2, 'r--', 'LineWidth', 2); 
legend([aa1, trend_line1, aa2, trend_line2], {'AR1 Data', 'AR1 Trend', '\lambda Data', '\lambda Trend'}, 'location', 'SE');
legend([aa1, aa2], grouporder,'location','SE', 'NumColumns', 2);
legend('boxoff');
set(gca,'FontSize', 14,'FontName', 'Times New Roman'); 
set(gca,'linewidth',2);
xlabel('Year') 
set(gcf,"Position",[10 10 800 200]);
print(gcf, '河道平均AR10911', '-dtiff', '-r300')

%% subfigure g
clear all
rng('default')
data = readtable('Consolidated_Residual_Levels0913.xlsx', 'Sheet', 2);
Subsystem1 = {'MK', 'GZ', 'ZY', 'DLS'};
Subsystem2 = {'SS', 'SD', 'SSJ', 'BSW', 'HM'};
Subsystem3 = {'HP', 'DS', 'SSK', 'SSW', 'NS'};
Subsystem4 = {'HS', 'BJ', 'HJ', 'SJK', 'GC'};
Subsystem5 = {'ZD', 'LS', 'FBC', 'NH', 'RQ', 'MA'};
condition_names = {'Subsystem1', 'Subsystem2', 'Subsystem3', 'Subsystem4', 'Subsystem5'};

data_subsystem{1} = reshape(data{:, Subsystem1}, [], 1);  
data_subsystem{2} = reshape(data{:, Subsystem2}, [], 1);  
data_subsystem{3} = reshape(data{:, Subsystem3}, [], 1);  
data_subsystem{4} = reshape(data{:, Subsystem4}, [], 1);  
data_subsystem{5} = reshape(data{:, Subsystem5}, [], 1); 
all_data = [data_subsystem{1}; data_subsystem{2}; data_subsystem{3}; data_subsystem{4}; data_subsystem{5}];
group_labels = [repmat({'Subsystem1'}, length(data_subsystem{1}), 1); 
                repmat({'Subsystem2'}, length(data_subsystem{2}), 1); 
                repmat({'Subsystem3'}, length(data_subsystem{3}), 1); 
                repmat({'Subsystem4'}, length(data_subsystem{4}), 1); 
                repmat({'Subsystem5'}, length(data_subsystem{5}), 1)];
c =  [0.8392, 0.1529, 0.1569;...
      0.585098039215686,0.570588235294118,0.774509803921569;...
      0.994117647058824,0.405882352941177,0.094117647058824;...
      0.0156862745098039,0.282352941176471,0.419607843137255;...
      0.305882352941177,0.541176470588235,0.317647058823529];
figure('Name', 'Subsystem Violin Plot','WindowStyle','normal');
h = daviolinplot(all_data, 'groups', group_labels, 'xtlabels', condition_names, 'violin', 'full', ...
     'colors', c, 'boxcolors', 'w', 'scatter', 0, 'jitter', 0, 'violinalpha', 0.5, 'scatteralpha', 0.5);  
ylabel('AR1', 'FontName', 'Times New Roman');
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'linewidth', 2);
set(gcf,"Position",[50 50 800 200]);
print(gcf, 'AR1_Violin_Plot', '-dtiff', '-r300');
