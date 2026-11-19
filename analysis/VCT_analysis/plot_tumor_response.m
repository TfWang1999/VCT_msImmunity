function F = plot_tumor_response(export_files, visible_state)
% Plot the Fig. 5 panels from the shared virtual-cohort report.
arguments
    export_files (1,1) logical = false
    visible_state = 'on'
end
figure_dir = fullfile(fileparts(mfilename('fullpath')),'results');
loaded = load(fullfile(figure_dir,'vct_analysis_results.mat'),'report');
regimens = loaded.report.regimens([1 3 2 4]);
display_names = {'Xel','Xel + aPD1','Xel + RT','Xel + RT + aPD1'};
volume_values = cell(1,4);
proportions = zeros(4,3);
for k = 1:4
    volume_values{k} = 100*regimens(k).volume_change(:);
    classes = regimens(k).response_class;
    proportions(k,:) = 100*[mean(strcmp(classes,'PD')) ...
        mean(strcmp(classes,'SD')) mean(strcmp(classes,'PR'))];
end
pairs = nchoosek(1:4,2);
significance = cell(size(pairs,1),1);
for p = 1:size(pairs,1)
    p_value = ranksum(volume_values{pairs(p,1)},volume_values{pairs(p,2)});
    if p_value < 0.001
        significance{p} = '***';
    elseif p_value < 0.01
        significance{p} = '**';
    elseif p_value < 0.05
        significance{p} = '*';
    else
        significance{p} = 'ns';
    end
end
tumor_response_handle = gobjects(1,3);

%% 子图 A：小提琴图
k = 1;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.SizeCM = [9 7.8]; % 独立子图的 [宽 高]，单位 cm；同时控制图窗与导出图片尺寸
panel_options.AxesPosition = [0.205 0.240 0.755 0.670]; % 绘图区在本图窗内的 [左 下 宽 高]，归一化比例（无单位）；通常取 0–1
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 7.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 10; % 横轴和纵轴标题的字号，单位 pt
panel_options.AnnotationFontSize = 7; % 图内数值、显著性或分类标注的字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.Resolution = 600; % PNG 导出分辨率，单位 dpi；不改变矢量图的物理尺寸
panel_options.ViolinColors = [0.57 0.78 0.68;0.91 0.66 0.49;0.66 0.64 0.80;0.91 0.76 0.38]; % 四个小提琴的 RGB 配色（每项 0–1），行顺序与横轴治疗组一致
panel_options.TreatmentLabels = {'Xel','Xel + aPD1', ...
    'Xel + RT','Xel + RT + aPD1'}; % 横轴治疗组名称列表（无单位），顺序与各组数据一致
panel_options.BracketHeights = [205 230 255 280 305 330]; % 六条两两比较括号的纵轴位置，单位为肿瘤体积变化的百分点；按 pairs 顺序排列
panel_options.DensityPoints = 160; % 小提琴核密度曲线的采样点数，单位 个；增加后轮廓更细密
panel_options.ViolinWidth = 0.34; % 小提琴最大半宽，单位为横轴组间距（相邻组中心距离为 1）
panel_options.ViolinAlpha = 0.78; % 小提琴填充不透明度，无单位；0 完全透明，1 完全不透明
panel_options.XLimits = [0.5 4.5]; % 横轴 [最小值 最大值]，单位：治疗组序号；[] 表示沿用自动范围
panel_options.YLimits = [-105 350]; % 纵轴 [最小值 最大值]，单位：肿瘤体积变化百分数（%）；[] 表示沿用自动范围
panel_options.YTicks = -100:100:300; % 纵轴刻度位置，单位：肿瘤体积变化百分数（%）；冒号写法为 起点:间隔:终点

panel_options.XLabel = 'Treatment group'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = {'Tumor volume change (%)'}; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行

panel_options.XTicks = 1:4; % 横轴刻度位置，单位为治疗组序号；1–4 对应四组
panel_options.WhiskerLineWidth = 0.65; % 小提琴图内部上下须线的粗细，单位 pt
panel_options.MedianLineWidth = 0.9; % 小提琴图内部中位数横线的粗细，单位 pt
panel_options.BoxAlpha = 0.80; % 小提琴内部四分位箱的不透明度，无单位；0 完全透明，1 完全不透明
panel_options.BoxHalfWidth = 0.075; % 小提琴内部四分位箱的半宽，单位为横轴组间距（相邻组中心距离为 1）
panel_options.BracketHeight = 2.5; % 显著性括号两端短竖线的高度，单位为纵轴百分点
panel_options.SignificanceOffset = 9; % 显著性文字相对括号向上偏移的距离，单位为纵轴百分点

tumor_response_handle(1) = figure('Visible',visible_state,'Color','white', ...
        'Units','centimeters','Position',[2 2 panel_options.SizeCM]);
palette = panel_options.ViolinColors;
data_top=max(cellfun(@max,volume_values));
if data_top>min(panel_options.BracketHeights)-25
    spacing=max(25,0.06*(data_top+100));
    panel_options.BracketHeights=data_top+(1:6)*spacing;
    panel_options.YLimits(2)=panel_options.BracketHeights(end)+spacing;
    tick_step=100*max(1,ceil(panel_options.YLimits(2)/600));
    panel_options.YTicks=[-100,0:tick_step:panel_options.YLimits(2)];
end
ax_a = axes(tumor_response_handle(1),'Units','normalized','Position',panel_options.AxesPosition); %#ok<LAXES>
hold(ax_a,'on');
for k = 1:4
    ax = ax_a; y = volume_values{k}; x_center = k; color = palette(k,:);
    % Resolve the concentrated response peak as well as the long upper tail.
    [~,base_support,bandwidth] = ksdensity(y,'NumPoints',panel_options.DensityPoints);
    dense_support=quantile(y,linspace(0,1,4*panel_options.DensityPoints));
    support=unique([base_support(:);dense_support(:)])';
    [density,support] = ksdensity(y,support,'Bandwidth',bandwidth);
    density = panel_options.ViolinWidth*density/max(density);
    patch(ax,[x_center-density fliplr(x_center+density)], ...
        [support fliplr(support)],color,'FaceAlpha',panel_options.ViolinAlpha, ...
        'EdgeColor',color*0.72,'LineWidth',panel_options.LineWidth);
    q = prctile(y,[5 25 50 75 95]);
    plot(ax,[x_center x_center],[q(1) q(5)],'k-','LineWidth',panel_options.WhiskerLineWidth);
    patch(ax,[x_center-panel_options.BoxHalfWidth x_center+panel_options.BoxHalfWidth x_center+panel_options.BoxHalfWidth x_center-panel_options.BoxHalfWidth], ...
        [q(2) q(2) q(4) q(4)],[1 1 1],'FaceAlpha',panel_options.BoxAlpha, ...
        'EdgeColor','k','LineWidth',panel_options.LineWidth);
    plot(ax,[x_center-panel_options.BoxHalfWidth x_center+panel_options.BoxHalfWidth],[q(3) q(3)],'k-','LineWidth',panel_options.MedianLineWidth);
end
set(ax_a,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth, ...
    'Box','on','TickDir','out','XLim',panel_options.XLimits,'YLim',panel_options.YLimits, ...
    'XTick',panel_options.XTicks,'YTick',panel_options.YTicks);
xticklabels(ax_a,panel_options.TreatmentLabels);
xtickangle(ax_a,0);
ax_a.XAxis.FontSize = panel_options.TickFontSize;
ylabel(ax_a,panel_options.YLabel,'FontName',panel_options.FontName,'FontSize',panel_options.LabelFontSize);
xlabel(ax_a,panel_options.XLabel,'FontName',panel_options.FontName,'FontSize',panel_options.LabelFontSize);
bracket_heights = panel_options.BracketHeights;
for p = 1:size(pairs,1)
    ax = ax_a; x1 = pairs(p,1); x2 = pairs(p,2);
    y = bracket_heights(p); label = significance{p};
    h = panel_options.BracketHeight;
    plot(ax,[x1 x1 x2 x2],[y-h y y y-h],'k-','LineWidth',panel_options.LineWidth);
    text(ax,(x1+x2)/2,y+panel_options.SignificanceOffset,label,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontName',panel_options.FontName,'FontSize',panel_options.AnnotationFontSize,'FontWeight','bold');
end


if export_files
    output_stem = fullfile(figure_dir,'tumor_change_violin');
    for extension = {'.pdf','.svg'}
        exportgraphics(tumor_response_handle(1),[output_stem extension{1}],'ContentType','vector', ...
            'BackgroundColor','white','Units','centimeters', ...
            'Width',panel_options.SizeCM(1),'Height',panel_options.SizeCM(2), ...
            'Padding','figure','PreserveAspectRatio','on');
    end
    exportgraphics(tumor_response_handle(1),[output_stem '.png'],'Resolution',panel_options.Resolution, ...
        'BackgroundColor','white','Units','centimeters', ...
        'Width',panel_options.SizeCM(1),'Height',panel_options.SizeCM(2), ...
            'Padding','figure','PreserveAspectRatio','on');
end
if strcmp(visible_state,'off'); close(tumor_response_handle(1)); end

%% 子图 B：PR/SD/PD 比例图
k = 2;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.SizeCM = [9 7.8]; % 独立子图的 [宽 高]，单位 cm；同时控制图窗与导出图片尺寸
panel_options.AxesPosition = [0.145 0.240 0.815 0.670]; % 绘图区在本图窗内的 [左 下 宽 高]，归一化比例（无单位）；通常取 0–1
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.xTickFontSize = 8; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.yTickFontSize = 10; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.TickFontSize = 10; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 10; % 横轴和纵轴标题的字号，单位 pt
panel_options.LegendFontSize = 10; % 图例文字字号，单位 pt
panel_options.xLabelFontSize= 10; % 图例文字字号，单位 pt
panel_options.yLabelFontSize= 10; % 图例文字字号，单位 pt
panel_options.AnnotationFontSize = 8; % 图内数值、显著性或分类标注的字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.Resolution = 600; % PNG 导出分辨率，单位 dpi；不改变矢量图的物理尺寸
panel_options.OutcomeColors = [0.89 0.43 0.36;0.36 0.62 0.78;0.47 0.72 0.43]; % 堆叠柱的 RGB 配色（每项 0–1），行顺序为 PD、SD、PR
panel_options.TreatmentLabels = {'Xel','Xel + aPD1', ...
    'Xel + RT','Xel + RT + aPD1'}; % 横轴治疗组名称列表（无单位），顺序与各组数据一致
panel_options.BarWidth = 0.5; % 柱形宽度占可用组间距的比例，无单位；通常取 0–1
panel_options.GridAlpha = 0.18; % 网格线不透明度，无单位；0 完全透明，1 完全不透明
panel_options.XLimits = [0.45 4.55]; % 横轴 [最小值 最大值]，单位：治疗组序号；[] 表示沿用自动范围
panel_options.YLimits = [0 100]; % 纵轴 [最小值 最大值]，单位：比例百分数（%）；[] 表示沿用自动范围
panel_options.YTicks = 0:20:100; % 纵轴刻度位置，单位：比例百分数（%）；冒号写法为 起点:间隔:终点

panel_options.XLabel = 'Treatment group'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Proportion (%)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.LegendTokenSize = [10 6]; % 图例中线段或色块标记的 [宽 高]，单位 pt
panel_options.BarEdgeColor = [0.2 0.2 0.2]; % 柱形边框颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BarLineWidth = 0.6; % 柱形边框线粗细，单位 pt
panel_options.LegendPosition = [0.238 0.915 0.642 0.041]; % 图例在整幅图窗内的 [左 下 宽 高]，归一化比例（无单位）

panel_options.LegendLabels = {'PD','SD','PR'}; % 图例文字列表（无单位），按对应曲线或色块的绘制顺序排列
panel_options.XTicks = 1:4; % 横轴刻度位置，单位为治疗组序号；1–4 对应四组

tumor_response_handle(2) = figure('Visible',visible_state,'Color','white', ...
        'Units','centimeters','Position',[2 2 panel_options.SizeCM]);
ax_b = axes(tumor_response_handle(2),'Units','normalized','Position',panel_options.AxesPosition); %#ok<LAXES>
outcome_colors = panel_options.OutcomeColors;
b = bar(ax_b,proportions,'stacked','BarWidth',panel_options.BarWidth);
for c = 1:3
    b(c).FaceColor = outcome_colors(c,:);
    b(c).EdgeColor = panel_options.BarEdgeColor;
    b(c).LineWidth = panel_options.BarLineWidth;
end
hold(ax_b,'on');
for k = 1:4
    cumulative = 0;
    for c = 1:3
        value = proportions(k,c);
        if value > 5
            value_label=sprintf('%.1f%%',value);
            if value==100; value_label='100%'; end
            text(ax_b,k,cumulative+value/2,value_label, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontName',panel_options.FontName,'FontSize',panel_options.AnnotationFontSize,'Color','white','FontWeight','bold');
        end
        cumulative = cumulative+value;
    end
end
set(ax_b,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth,'Box','off', ...
    'TickDir','out','XLim',panel_options.XLimits,'YLim',panel_options.YLimits,'XTick',panel_options.XTicks, ...
    'YTick',panel_options.YTicks,'YGrid','on','GridAlpha',panel_options.GridAlpha,'Layer','top');
xticklabels(ax_b,panel_options.TreatmentLabels);
xtickangle(ax_b,0);
ax_b.XAxis.FontSize = panel_options.xTickFontSize;
ax_b.YAxis.FontSize = panel_options.yTickFontSize;
ylabel(ax_b,panel_options.YLabel,'FontName',panel_options.FontName,'FontSize',panel_options.yLabelFontSize,'FontWeight','normal');
xlabel(ax_b,panel_options.XLabel,'FontName',panel_options.FontName,'FontSize',panel_options.xLabelFontSize,'FontWeight','normal');

lgd = legend(ax_b,b,panel_options.LegendLabels,'Box','off', ...
    'FontName',panel_options.FontName,'FontSize',panel_options.LegendFontSize,'Orientation','horizontal','Units','normalized','Position',panel_options.LegendPosition);
lgd.ItemTokenSize = panel_options.LegendTokenSize;


if export_files
    output_stem = fullfile(figure_dir,'response_proportions');
    for extension = {'.pdf','.svg'}
        exportgraphics(tumor_response_handle(2),[output_stem extension{1}],'ContentType','vector', ...
            'BackgroundColor','white','Units','centimeters', ...
            'Width',panel_options.SizeCM(1),'Height',panel_options.SizeCM(2), ...
            'Padding','figure','PreserveAspectRatio','on');
    end
    exportgraphics(tumor_response_handle(2),[output_stem '.png'],'Resolution',panel_options.Resolution, ...
        'BackgroundColor','white','Units','centimeters', ...
        'Width',panel_options.SizeCM(1),'Height',panel_options.SizeCM(2), ...
            'Padding','figure','PreserveAspectRatio','on');
end
if strcmp(visible_state,'off'); close(tumor_response_handle(2)); end

%% 子图 C：排序分布
k = 3;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.SizeCM = [17.5 9]; % 独立子图的 [宽 高]，单位 cm；同时控制图窗与导出图片尺寸
panel_options.AxesPosition = [0.095 0.140 0.875 0.760]; % 绘图区在本图窗内的 [左 下 宽 高]，归一化比例（无单位）；通常取 0–1
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 10; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 10; % 横轴和纵轴标题的字号，单位 pt
panel_options.LegendFontSize = 10; % 图例文字字号，单位 pt
panel_options.AnnotationFontSize = 10; % 图内数值、显著性或分类标注的字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.Resolution = 600; % PNG 导出分辨率，单位 dpi；不改变矢量图的物理尺寸
panel_options.DistributionColors = [27 158 119;117 112 179;217 95 2;225 171 2]/255; % 排序分布的 RGB 配色；此处将 0–255 色值除以 255，行顺序与图例一致
panel_options.DisplayRanks = 100; % 每组显示的排序抽样点数，单位 个；仅控制显示，不改变统计所用队列
panel_options.DistributionAlpha = 0.6; % 排序分布填充不透明度，无单位；0 完全透明，1 完全不透明
panel_options.ThresholdColor = [0.2 0.2 0.2]; % 分类分界线颜色；RGB 三元组，每项 0–1（无单位）
panel_options.LegendLocation = [0.5,0.75]; % 图例位置（字符串，无单位）：north 上方、northeast 右上、northwest 左上
panel_options.XLimits = [0 panel_options.DisplayRanks+12]; % 横轴 [最小值 最大值]，单位：显示排序序号；[] 表示沿用自动范围
panel_options.YLimits = [-100 200]; % 纵轴 [最小值 最大值]，单位：肿瘤体积变化百分数（%）；[] 表示沿用自动范围
panel_options.YTicks = -100:50:200; % 纵轴刻度位置，单位：肿瘤体积变化百分数（%）；冒号写法为 起点:间隔:终点

panel_options.XLabel = 'Sample index'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = {'Tumor volume change (%)'}; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.LegendTokenSize = [12 6]; % 图例中线段或色块标记的 [宽 高]，单位 pt

panel_options.CurveLineWidth = 0.9; % 曲线粗细，单位 pt；散点图中控制分箱均值线
panel_options.ThresholdLineWidth = 1.0; % 疗效分类分界线粗细，单位 pt

panel_options.ResponseLabelY = [100 0 -82.5]; % PD、SD、PR 标注的纵坐标，单位为肿瘤体积变化百分点
panel_options.ResponseLabelXOffset = 4; % 分类文字相对最后一个显示点向右偏移的距离，单位为排序序号

tumor_response_handle(3) = figure('Visible',visible_state,'Color','white', ...
        'Units','centimeters','Position',[2 2 panel_options.SizeCM]);
data_top=max(cellfun(@max,volume_values));
if data_top>panel_options.YLimits(2)
    panel_options.YLimits(2)=100*ceil(1.05*data_top/100);
    panel_options.YTicks=-100:100:panel_options.YLimits(2);
end
ax_c = axes(tumor_response_handle(3),'Units','normalized','Position',panel_options.AxesPosition); %#ok<LAXES>
hold(ax_c,'on');
distribution_order = [1 3 2 4];
distribution_names = display_names(distribution_order);
distribution_colors = panel_options.DistributionColors;
n_display = min(panel_options.DisplayRanks,min(cellfun(@numel,volume_values)));
ranked_values = zeros(n_display,4);
for k = 1:4
    sorted_values = sort(volume_values{distribution_order(k)},'descend');
    selected_ranks = round(linspace(1,numel(sorted_values),n_display));
    ranked_values(:,k) = sorted_values(selected_ranks);
end
for i = 1:n_display
    row = ranked_values(i,:);
    positive = find(row>0);
    negative = find(row<0);
    [~,positive_order] = sort(row(positive),'descend');
    [~,negative_order] = sort(row(negative),'ascend');
    draw_order = [positive(positive_order),negative(negative_order)];
    for k = draw_order
        patch(ax_c,i+[-0.5 0.5 0.5 -0.5],[0 0 row(k) row(k)], ...
            distribution_colors(k,:),'FaceAlpha',panel_options.DistributionAlpha,'EdgeColor','none','HandleVisibility','off');
    end
end
curve_handles = gobjects(1,4);
for k = 1:4
    curve_handles(k) = plot(ax_c,1:n_display,ranked_values(:,k), ...
        'Color',distribution_colors(k,:),'LineWidth',panel_options.CurveLineWidth);
end
% 等效直径 +20% / -30% 对应的体积变化：+72.8% / -65.7%。
for threshold = 100*([1.2 0.7].^3-1)
    yline(ax_c,threshold,'--','Color',panel_options.ThresholdColor,'LineWidth',panel_options.ThresholdLineWidth,'HandleVisibility','off');
end
response_labels = {'PD','SD','PR'};
for label_index = 1:3
    text(ax_c,n_display+panel_options.ResponseLabelXOffset,panel_options.ResponseLabelY(label_index),response_labels{label_index},'FontName',panel_options.FontName,'FontSize',panel_options.AnnotationFontSize);
end
set(ax_c,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth,'Box','off', ...
    'TickDir','out','XLim',panel_options.XLimits,'YLim',panel_options.YLimits, ...
    'XTick',[],'YTick',panel_options.YTicks,'Layer','top');
xlabel(ax_c,panel_options.XLabel,'FontName',panel_options.FontName,'FontSize',panel_options.LabelFontSize,'FontWeight','normal');
ylabel(ax_c,panel_options.YLabel,'FontName',panel_options.FontName,'FontSize',panel_options.LabelFontSize,'FontWeight','normal');

lgd = legend(ax_c,curve_handles,distribution_names,'Location','none','Box','off', ...
    'NumColumns', 2,'FontName',panel_options.FontName,'FontSize',panel_options.LegendFontSize);
lgd.ItemTokenSize = panel_options.LegendTokenSize;

lgd.Units = 'normalized';
pos = lgd.Position;
pos(1:2) = panel_options.LegendLocation;
lgd.Position = pos;

if export_files
    output_stem = fullfile(figure_dir,'tumor_change_distribution');
    for extension = {'.pdf','.svg'}
        exportgraphics(tumor_response_handle(3),[output_stem extension{1}],'ContentType','vector', ...
            'BackgroundColor','white','Units','centimeters', ...
            'Width',panel_options.SizeCM(1),'Height',panel_options.SizeCM(2), ...
            'Padding','figure','PreserveAspectRatio','on');
    end
    exportgraphics(tumor_response_handle(3),[output_stem '.png'],'Resolution',panel_options.Resolution, ...
        'BackgroundColor','white','Units','centimeters', ...
        'Width',panel_options.SizeCM(1),'Height',panel_options.SizeCM(2), ...
            'Padding','figure','PreserveAspectRatio','on');
end
if strcmp(visible_state,'off'); close(tumor_response_handle(3)); end

F = tumor_response_handle;
end
