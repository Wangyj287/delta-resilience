

%% subfigure i
file_names = {'condition1-1.tif', 'condition2-1.tif', 'condition3-1.tif', 'condition4-1.tif'};
c = [0.8392, 0.1529, 0.1569;...
     0.5851, 0.5706, 0.7745;...
     0.9941, 0.4059, 0.0941;...
     0.0157, 0.2824, 0.4196;...
     0.3059, 0.5412, 0.3176];
all_data = [];
group_labels = [];
condition_names = {'Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'};
for i = 1:length(file_names)
    [data, R] = readgeoraster(file_names{i});
    valid_data = data(data > 0);
    all_data = [all_data; valid_data]; 
    group_labels = [group_labels; i * ones(length(valid_data), 1)]; 
end

figure('Name', 'Subsystem Violin Plot', 'WindowStyle', 'normal');
h = daviolinplot(all_data, 'groups', group_labels, 'xtlabels', condition_names, ...
    'violin', 'full', 'colors', c, 'boxcolors', 'w', 'scatter', 0, 'jitter', 0, ...
    'violinalpha', 0.5, 'scatteralpha', 0.5);
ylabel('AR1', 'FontName', 'Times New Roman');
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'linewidth', 2);
set(gcf, "Position", [50 50 1200 400]);
print(gcf, 'AR1_Violin_Plot1128', '-dtiff', '-r300');
