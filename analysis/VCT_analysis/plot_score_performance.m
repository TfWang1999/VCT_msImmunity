function figures = plot_score_performance(score, export_files, visible_state)
% Plot Fig. 6 (msImmunity and AUC), Fig. S4 (tsImmunity), or Fig. S5 (molsImmunity).
arguments
    score (1,1) string {mustBeMember(score,["msImmunity","tsImmunity","molsImmunity"])} = "msImmunity"
    export_files (1,1) logical = false
    visible_state = 'on'
end
results_dir = fullfile(fileparts(mfilename('fullpath')),'results');
loaded = load(fullfile(results_dir,'vct_analysis_results.mat'),'report');
report = loaded.report;
score_index = find(strcmp(report.score_names,score),1);

panel_size = [8.3005 5.40];
axes_position = [0.205 0.176 0.741 0.685];
titles = {'Xelox','Xelox + RT','Xelox + Anti-PD-1','Xelox + RT + Anti-PD-1'};
point_colors = [0.4 0.6682 0.8446; 0.91 0.5950 0.4588; ...
    0.9574 0.8164 0.4750; 0.6964 0.5104 0.7336];
curve_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
regimen_suffixes = {'xelox','xelox_rt','xelox_anti_pd1','xelox_rt_anti_pd1'};
if score == "msImmunity"
    figures = gobjects(1,5);
    for k = 1:4
        figures(k) = figure('Visible',visible_state,'Color','white','Units','centimeters', ...
            'Position',[2 2 panel_size]);
        ax = axes(figures(k),'Position',axes_position);
        plot_response(ax,report.regimens(k),score_index,score,titles{k}, ...
            'Tumor regression ratio (%)',point_colors(k,:),curve_colors(k,:));
        if export_files
            export_panel(figures(k),fullfile(results_dir,['msImmunity_response_' regimen_suffixes{k}]),panel_size);
        end
    end
    figures(5) = plot_auc(report,visible_state);
    if export_files
        export_panel(figures(5),fullfile(results_dir,'predictive_auc'),[17.0011 7.5034]);
    end
else
    % The two supplementary figures use the Fig. 6 panels in a 2-by-2 grid.
    group_size = 2*panel_size;
    figures = figure('Visible',visible_state,'Color','white','Units','centimeters', ...
        'Position',[2 2 group_size]);
    offsets = [0 0.5; 0.5 0.5; 0 0; 0.5 0];
    for k = 1:4
        position = axes_position/2;
        position(1:2) = position(1:2)+offsets(k,:);
        ax = axes(figures,'Position',position);
        plot_response(ax,report.regimens(k),score_index,score,titles{k}, ...
            'Percentage change of tumor (%)',point_colors(k,:),curve_colors(k,:));
    end
    if export_files
        basenames = {'tsm_IS_final','Module_IS_final'};
        export_panel(figures,fullfile(results_dir,basenames{score_index-1}),group_size);
    end
end
if strcmp(visible_state,'off'); close(figures); end
end

function plot_response(ax, regimen, score_index, score, heading, y_label, point_color, curve_color)
x = regimen.scores(:,score_index);
y = 100*regimen.volume_change(:);
pr_threshold = 100*(0.7^3-1);
pd_threshold = 100*(1.2^3-1);
x_lim = [0 1.05*max(x)];
y_lim = [-100 200];
if max(y)>y_lim(2)
    y_lim(2) = 25*ceil(1.05*max(y)/25);
end
hold(ax,'on');
ends = [x_lim fliplr(x_lim)];
patch(ax,ends,[pr_threshold pr_threshold y_lim(2) y_lim(2)], ...
    [1.0 0.96 0.93],'EdgeColor','none','FaceAlpha',0.4);
patch(ax,ends,[y_lim(1) y_lim(1) pr_threshold pr_threshold], ...
    [0.94 0.97 1.0],'EdgeColor','none','FaceAlpha',0.4);
yline(ax,pr_threshold,'-','Color',[0.2 0.2 0.2],'LineWidth',0.8,'Alpha',0.8);
yline(ax,pd_threshold,'--','Color',[0.2 0.2 0.2],'LineWidth',0.65,'Alpha',0.6);
points = scatter(ax,x,y,5,point_color,'filled','MarkerFaceAlpha',0.65,'MarkerEdgeColor','none');

% Match the 20 quantile bins and the separate upper endpoint in Fig. 6.
edges = quantile(x,linspace(0,1,21));
[~,bins] = histc(x,edges); %#ok<HISTC>
[~,~,occupied_bins] = unique(bins);
x_mean = accumarray(occupied_bins,x,[],@mean);
y_mean = accumarray(occupied_bins,y,[],@mean);
curve = plot(ax,x_mean,y_mean,'Color',curve_color,'LineWidth',1.35);

labels = {'PD','SD','PR'};
label_heights = [0.70 0.345 0.057];
for j = 1:3
    text(ax,0.965,label_heights(j),labels{j},'Units','normalized','FontName','Arial', ...
        'FontSize',8.5,'FontWeight','bold','HorizontalAlignment','right');
end
legend(ax,[curve points],{'Binned mean','Virtual patients'}, ...
    'Location','northeast','Box','off','FontName','Arial','FontSize',8);
title(ax,heading,'FontName','Arial','FontSize',9.5,'FontWeight','bold');
xlabel(ax,score,'FontName','Arial','FontSize',9,'FontWeight','bold');
ylabel(ax,y_label,'FontName','Arial','FontSize',9,'FontWeight','bold');
set(ax,'FontName','Arial','FontSize',8.5,'LineWidth',0.55,'Box','on','TickDir','out', ...
    'XLim',x_lim,'YLim',y_lim,'YTick',-100:50:y_lim(2),'Layer','top');
end

function fig = plot_auc(report, visible_state)
regimen_order = [1 3 2 4];
score_order = [1 3 2];
auc_values = report.auc(score_order,regimen_order)';
labels = {'Xel','Xel + aPD1','Xel + RT','Xel + RT + aPD1'};
palette = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
fig = figure('Visible',visible_state,'Color','white','Units','centimeters', ...
    'Position',[2 2 17.0011 7.5034]);
ax = axes(fig,'Position',[0.079 0.153 0.894 0.720]);
bars = bar(ax,auc_values,'grouped','BarWidth',0.8);
hold(ax,'on');
for j = 1:3; bars(j).FaceColor = palette(j,:); end
drawnow;
for j = 1:3
    estimable = isfinite(bars(j).YData);
    text(ax,bars(j).XEndPoints(estimable),bars(j).YEndPoints(estimable)+0.008, ...
        compose('%.3f',bars(j).YData(estimable)),'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontName','Arial','FontSize',8.5);
end
for k = find(all(~isfinite(auc_values),2))'
    response_rate = 100*report.regimens(regimen_order(k)).response_rate;
    text(ax,k,0.65,{'Not estimable',sprintf('%.0f%% response',response_rate)}, ...
        'HorizontalAlignment','center','FontName','Arial','FontSize',8.5,'Color',[0.35 0.35 0.35]);
end
y_lim = [0.50 1.16];
y_ticks = 0.5:0.1:1;
if any(auc_values(:)<y_lim(1))
    y_lim(1) = 0;
    y_ticks = 0:0.2:1;
end
set(ax,'FontName','Arial','FontSize',8.5,'LineWidth',0.55,'Box','off','TickDir','out', ...
    'XTick',1:4,'XLim',[0.5 4.5],'YLim',y_lim,'YTick',y_ticks);
xticklabels(ax,labels);
xtickangle(ax,0);
xlabel(ax,'Treatment group','FontSize',9);
ylabel(ax,'Area under the ROC curve (AUC)','FontSize',9);
key = legend(ax,bars,report.score_names(score_order),'Orientation','horizontal', ...
    'Location','north','Box','off','FontName','Arial','FontSize',8);
key.ItemTokenSize = [12 7];
end

function export_panel(fig, path, dimensions)
options = {'BackgroundColor','white','Units','centimeters','Width',dimensions(1), ...
    'Height',dimensions(2),'Padding','figure','PreserveAspectRatio','on'};
for extension = {'.pdf','.svg'}
    exportgraphics(fig,[path extension{1}],'ContentType','vector',options{:});
end
exportgraphics(fig,[path '.png'],'Resolution',600,options{:});
end
