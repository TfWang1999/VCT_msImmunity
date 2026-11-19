function plot_robustness(export_files, visible_state)
% Plot Supplementary Fig. S8 from the saved robustness analyses.
arguments
    export_files = false
    visible_state = 'on'
end

script_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(script_dir,'results');
abundance_table = readtable(fullfile(results_dir,'abundance_robustness_summary.csv'),'TextType','string');
auc_table = readtable(fullfile(results_dir,'auc_bootstrap_summary.csv'),'TextType','string');
abundance_values = unique(abundance_table.abundance_constant);
regimen_names = unique(abundance_table.regimen,'stable');
score_names = unique(auc_table.score,'stable');

regimen_labels = {'Xel','Xel + RT','Xel + aPD1','Xel + RT + aPD1'};
regimen_colors = [0.00 0.45 0.70; 0.83 0.37 0.00; 0.00 0.60 0.50; 0.65 0.37 0.66];
fig = figure('Color','w','Visible',visible_state,'Units','centimeters','Position',[2 2 21 9]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
hold(ax,'on');
xline(ax,0.5,':','Color',[0.65 0.65 0.65],'HandleVisibility','off');
for k = 1:numel(regimen_names)
    current = sortrows(abundance_table(abundance_table.regimen==regimen_names(k),:), ...
        'abundance_constant');
    plot(ax,current.abundance_constant,current.auc_ms,'-o', ...
        'Color',regimen_colors(k,:),'MarkerFaceColor',regimen_colors(k,:), ...
        'MarkerSize',5,'LineWidth',1.6);
end
set(ax,'XTick',abundance_values,'XLim',[min(abundance_values)-0.02,max(abundance_values)+0.02], ...
    'YLim',[min(abundance_table.auc_ms)-0.015,1.005],'XGrid','on','YGrid','on');
xlabel(ax,'Cell-abundance approximation constant');
ylabel(ax,'AUC');
title(ax,'A  msImmunity discrimination');
legend(ax,regimen_labels,'Location','southoutside','Orientation','horizontal', ...
    'NumColumns',2,'Box','off','FontSize',9);

ax_auc = nexttile(layout);
hold(ax_auc,'on');
offsets = [-0.22 0 0.22];
for m = 1:numel(score_names)
    rows = auc_table.score == score_names(m);
    auc = auc_table.auc(rows);
    errorbar(ax_auc,auc,(1:4)+offsets(m),auc-auc_table.ci_low(rows), ...
        auc_table.ci_high(rows)-auc,'horizontal','o','Color',regimen_colors(m,:), ...
        'MarkerFaceColor',regimen_colors(m,:),'MarkerSize',5,'LineWidth',1.5,'CapSize',7);
end
set(ax_auc,'YTick',1:4,'YTickLabel',regimen_labels,'YDir','reverse', ...
    'YLim',[0.5 4.5],'XLim',[0.55 1.02],'XGrid','on');
xlabel(ax_auc,'AUC');
title(ax_auc,'B  AUC with 95% CI');
legend(ax_auc,cellstr(score_names),'Location','southoutside','Orientation','horizontal', ...
    'Box','off','FontSize',9);
set([ax,ax_auc],'FontName','Arial','FontSize',10,'TickDir','out','Box','off','GridAlpha',0.12);

if export_files
    output = fullfile(results_dir,'robustness_discrimination');
    exportgraphics(fig,[output '.pdf'],'ContentType','vector');
    exportgraphics(fig,[output '.svg'],'ContentType','vector');
    exportgraphics(fig,[output '.png'],'Resolution',300);
end
if strcmp(visible_state,'off')
    close(fig);
end
end
