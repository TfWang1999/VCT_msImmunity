%% 小鼠拟合：参数设置
% 本脚本只绘图；数值由同目录的计算脚本生成。
% 首次先运行“参数设置”和“初始化数据和画布”，随后任意子图节可单独 Ctrl+Enter。
% 每个子图节开头列出全部实际绘图参数；修改后 Ctrl+Enter 仅更新该子图。
% SizeCM 为独立图窗宽高；PositionCM 为组图内绘图区 [左 下 宽 高]，单位 cm；字号 pt。
save_mouse_fit = false;         % true：导出图片；false：只预览
visible_state = 'on';           % 'on'：显示图窗；'off'：后台绘制并关闭
mouse_fit_canvas_cm = [17.5 14]; % 整幅组图画布；各绘图区尺寸在对应分节设置

export_files = save_mouse_fit;
visible_state = validatestring(visible_state,{'on','off'});

mouse_fit_visible = visible_state;


%% 读取数据并初始化画布
entry_file = mfilename('fullpath');
if isempty(entry_file) && usejava('desktop')
    entry_file = matlab.desktop.editor.getActiveFilename;
end
script_dir = fileparts(entry_file);
addpath(genpath(fullfile(fileparts(script_dir),'src')));
loaded = load(fullfile(script_dir,'results','mouse_results.mat'),'results');
results = loaded.results;
cases = mouse_cases();
figure_dir = fullfile(script_dir,'results');
mouse_fit_output_dir = figure_dir;
if exist('mouse_fit_handle','var')
    delete(mouse_fit_handle(isgraphics(mouse_fit_handle,'figure')));
end
mouse_fit_results = results;
mouse_fit_cases = cases;
mouse_fit_handle = figure('Visible',visible_state,'Units','centimeters','Position',[5 5 mouse_fit_canvas_cm]);
mouse_fit_axes = gobjects(1,4);

%% 子图 A：LLC
i = 1;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [1.29600 7.81008 6.92431 5.69138]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 7.2; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 8.5; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 9.2; % 子图标题字号，单位 pt
panel_options.LegendFontSize = 7.2; % 图例文字字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.CurveLineWidth = 1; % 曲线粗细，单位 pt；散点图中控制分箱均值线
panel_options.MarkerSize = 35; % 散点标记面积，单位 pt²（scatter 的 SizeData）；数值越大，点越大
panel_options.CapSize = 3; % 误差棒末端横帽的长度，单位 pt
panel_options.Colors = {'#0072BD','#77AC30','#D95319','#7E2F8E'}; % 各实验组颜色列表（十六进制字符串，无单位）；组别顺序见本节图例与绘图代码
panel_options.Markers = {'s','d','^','o'}; % 各实验组的点形状（无单位）：s 方形、d 菱形、^ 三角形、o 圆形
panel_options.LegendLocation = 'northwest'; % 图例位置（字符串，无单位）：north 上方、northeast 右上、northwest 左上
panel_options.Title = 'LLC cells injection-induced mouse model'; % 当前子图标题文字（无单位）；可替换为自定义字符串
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：mm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (mm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.LegendTokenSize = [13 8]; % 图例中线段或色块标记的 [宽 高]，单位 pt

panel_options.LegendLabels = {'Control','Anti-PD-1 mAb'}; % 图例文字列表（无单位），按对应曲线或色块的绘制顺序排列

visible_state = mouse_fit_visible;
export_files = save_mouse_fit;
figure_dir = mouse_fit_output_dir;
assert(isgraphics(mouse_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(mouse_fit_axes(i)); delete(mouse_fit_axes(i)); end
ax = axes(mouse_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
mouse_fit_axes(i) = ax;
control_model = mouse_fit_results(i).evaluation.groups(1);
anti_model = mouse_fit_results(i).evaluation.groups(2);
plot(anti_model.time,anti_model.volume_mm3,'Color',panel_options.Colors{3},'LineWidth',panel_options.CurveLineWidth);
hold on
plot(control_model.time,control_model.volume_mm3,'Color',panel_options.Colors{1},'LineWidth',panel_options.CurveLineWidth);
hold on
control = mouse_fit_cases(i).groups(1);
anti = mouse_fit_cases(i).groups(2);
d1 = scatter(control.time,control.mean,panel_options.MarkerSize,panel_options.Markers{1},'filled', ...
    'MarkerEdgeColor',panel_options.Colors{1},'MarkerFaceColor',panel_options.Colors{1});
errorbar(control.time,control.mean,control.sd,'Color',panel_options.Colors{1}, ...
    'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
d2 = scatter(anti.time,anti.mean,panel_options.MarkerSize,panel_options.Markers{3},'filled', ...
    'MarkerEdgeColor',panel_options.Colors{3},'MarkerFaceColor',panel_options.Colors{3});
errorbar(anti.time,anti.mean,anti.sd,'Color',panel_options.Colors{3}, ...
    'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
lgd = legend(ax,[d1,d2],panel_options.LegendLabels, ...
    'Location',panel_options.LegendLocation,'FontName',panel_options.FontName,'FontSize',panel_options.LegendFontSize,'Box','on');
lgd.ItemTokenSize = panel_options.LegendTokenSize;
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize);
title(ax,panel_options.Title,'FontSize',panel_options.TitleFontSize,'FontWeight','normal');
set(ax,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth, ...
    'Box','on','TickDir','out','LooseInset',get(ax,'TightInset'));
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end

%% 子图 B：MC38
i = 2;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [9.63965 7.81008 6.92431 5.69138]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 7.2; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 8.5; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 9.2; % 子图标题字号，单位 pt
panel_options.LegendFontSize = 7.2; % 图例文字字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.CurveLineWidth = 1; % 曲线粗细，单位 pt；散点图中控制分箱均值线
panel_options.MarkerSize = 35; % 散点标记面积，单位 pt²（scatter 的 SizeData）；数值越大，点越大
panel_options.CapSize = 3; % 误差棒末端横帽的长度，单位 pt
panel_options.Colors = {'#0072BD','#77AC30','#D95319','#7E2F8E'}; % 各实验组颜色列表（十六进制字符串，无单位）；组别顺序见本节图例与绘图代码
panel_options.Markers = {'s','d','^','o'}; % 各实验组的点形状（无单位）：s 方形、d 菱形、^ 三角形、o 圆形
panel_options.LegendLocation = 'northwest'; % 图例位置（字符串，无单位）：north 上方、northeast 右上、northwest 左上
panel_options.Title = 'MC38 cells injection-induced mouse model'; % 当前子图标题文字（无单位）；可替换为自定义字符串
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：mm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (mm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.LegendTokenSize = [13 8]; % 图例中线段或色块标记的 [宽 高]，单位 pt

panel_options.LegendLabels = {'Control','Anti-PD-1 mAb'}; % 图例文字列表（无单位），按对应曲线或色块的绘制顺序排列

visible_state = mouse_fit_visible;
export_files = save_mouse_fit;
figure_dir = mouse_fit_output_dir;
assert(isgraphics(mouse_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(mouse_fit_axes(i)); delete(mouse_fit_axes(i)); end
ax = axes(mouse_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
mouse_fit_axes(i) = ax;
control_model = mouse_fit_results(i).evaluation.groups(1);
anti_model = mouse_fit_results(i).evaluation.groups(2);
plot(anti_model.time,anti_model.volume_mm3,'Color',panel_options.Colors{3},'LineWidth',panel_options.CurveLineWidth);
hold on
plot(control_model.time,control_model.volume_mm3,'Color',panel_options.Colors{1},'LineWidth',panel_options.CurveLineWidth);
hold on
control = mouse_fit_cases(i).groups(1);
anti = mouse_fit_cases(i).groups(2);
d1 = scatter(control.time,control.mean,panel_options.MarkerSize,panel_options.Markers{1},'filled', ...
    'MarkerEdgeColor',panel_options.Colors{1},'MarkerFaceColor',panel_options.Colors{1});
errorbar(control.time,control.mean,control.sd,'Color',panel_options.Colors{1}, ...
    'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
d2 = scatter(anti.time,anti.mean,panel_options.MarkerSize,panel_options.Markers{3},'filled', ...
    'MarkerEdgeColor',panel_options.Colors{3},'MarkerFaceColor',panel_options.Colors{3});
errorbar(anti.time,anti.mean,anti.sd,'Color',panel_options.Colors{3}, ...
    'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
lgd = legend(ax,[d1,d2],panel_options.LegendLabels, ...
    'Location',panel_options.LegendLocation,'FontName',panel_options.FontName,'FontSize',panel_options.LegendFontSize,'Box','on');
lgd.ItemTokenSize = panel_options.LegendTokenSize;
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize);
title(ax,panel_options.Title,'FontSize',panel_options.TitleFontSize,'FontWeight','normal');
set(ax,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth, ...
    'Box','on','TickDir','out','LooseInset',get(ax,'TightInset'));
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end

%% 子图 C：A375
i = 3;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [1.29600 0.88779 6.92431 5.69138]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 7.2; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 8.5; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 9.2; % 子图标题字号，单位 pt
panel_options.LegendFontSize = 7.2; % 图例文字字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.CurveLineWidth = 1; % 曲线粗细，单位 pt；散点图中控制分箱均值线
panel_options.MarkerSize = 35; % 散点标记面积，单位 pt²（scatter 的 SizeData）；数值越大，点越大
panel_options.CapSize = 3; % 误差棒末端横帽的长度，单位 pt
panel_options.Colors = {'#0072BD','#77AC30','#D95319','#7E2F8E'}; % 各实验组颜色列表（十六进制字符串，无单位）；组别顺序见本节图例与绘图代码
panel_options.Markers = {'s','d','^','o'}; % 各实验组的点形状（无单位）：s 方形、d 菱形、^ 三角形、o 圆形
panel_options.LegendLocation = 'northwest'; % 图例位置（字符串，无单位）：north 上方、northeast 右上、northwest 左上
panel_options.Title = 'A375 cells injection-induced mouse model'; % 当前子图标题文字（无单位）；可替换为自定义字符串
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：mm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (mm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.LegendTokenSize = [13 8]; % 图例中线段或色块标记的 [宽 高]，单位 pt

panel_options.LegendLabels = {'Control','Anti-PD-1 mAb'}; % 图例文字列表（无单位），按对应曲线或色块的绘制顺序排列

visible_state = mouse_fit_visible;
export_files = save_mouse_fit;
figure_dir = mouse_fit_output_dir;
assert(isgraphics(mouse_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(mouse_fit_axes(i)); delete(mouse_fit_axes(i)); end
ax = axes(mouse_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
mouse_fit_axes(i) = ax;
control_model = mouse_fit_results(i).evaluation.groups(1);
anti_model = mouse_fit_results(i).evaluation.groups(2);
plot(anti_model.time,anti_model.volume_mm3,'Color',panel_options.Colors{3},'LineWidth',panel_options.CurveLineWidth);
hold on
plot(control_model.time,control_model.volume_mm3,'Color',panel_options.Colors{1},'LineWidth',panel_options.CurveLineWidth);
hold on
control = mouse_fit_cases(i).groups(1);
anti = mouse_fit_cases(i).groups(2);
d1 = scatter(control.time,control.mean,panel_options.MarkerSize,panel_options.Markers{1},'filled', ...
    'MarkerEdgeColor',panel_options.Colors{1},'MarkerFaceColor',panel_options.Colors{1});
errorbar(control.time,control.mean,control.sd,'Color',panel_options.Colors{1}, ...
    'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
d2 = scatter(anti.time,anti.mean,panel_options.MarkerSize,panel_options.Markers{3},'filled', ...
    'MarkerEdgeColor',panel_options.Colors{3},'MarkerFaceColor',panel_options.Colors{3});
errorbar(anti.time,anti.mean,anti.sd,'Color',panel_options.Colors{3}, ...
    'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
lgd = legend(ax,[d1,d2],panel_options.LegendLabels, ...
    'Location',panel_options.LegendLocation,'FontName',panel_options.FontName,'FontSize',panel_options.LegendFontSize,'Box','on');
lgd.ItemTokenSize = panel_options.LegendTokenSize;
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize);
title(ax,panel_options.Title,'FontSize',panel_options.TitleFontSize,'FontWeight','normal');
set(ax,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth, ...
    'Box','on','TickDir','out','LooseInset',get(ax,'TightInset'));
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end

%% 子图 D：CT26
i = 4;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [9.63965 0.88779 6.92431 5.69138]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 7.2; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 8.5; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 9.2; % 子图标题字号，单位 pt
panel_options.CT26LegendFontSize = 6.5; % CT26 面板图例文字字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.CurveLineWidth = 1; % 曲线粗细，单位 pt；散点图中控制分箱均值线
panel_options.MarkerSize = 35; % 散点标记面积，单位 pt²（scatter 的 SizeData）；数值越大，点越大
panel_options.CapSize = 3; % 误差棒末端横帽的长度，单位 pt
panel_options.Colors = {'#0072BD','#77AC30','#D95319','#7E2F8E'}; % 各实验组颜色列表（十六进制字符串，无单位）；组别顺序见本节图例与绘图代码
panel_options.Markers = {'s','d','^','o'}; % 各实验组的点形状（无单位）：s 方形、d 菱形、^ 三角形、o 圆形
panel_options.LegendLocation = 'northwest'; % 图例位置（字符串，无单位）：north 上方、northeast 右上、northwest 左上
panel_options.Title = 'CT26 cells injection-induced mouse model'; % 当前子图标题文字（无单位）；可替换为自定义字符串
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：mm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (mm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.LegendTokenSize = [12 7]; % 图例中线段或色块标记的 [宽 高]，单位 pt

panel_options.LegendLabels = {'Control','Low PD-L1 expression', ...
    'Control + Anti-PD-1 mAb',sprintf('Low PD-L1 expression\n+ Anti-PD-1 mAb')}; % 图例文字列表（无单位），按对应曲线或色块的绘制顺序排列

visible_state = mouse_fit_visible;
export_files = save_mouse_fit;
figure_dir = mouse_fit_output_dir;
assert(isgraphics(mouse_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(mouse_fit_axes(i)); delete(mouse_fit_axes(i)); end
ax = axes(mouse_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
mouse_fit_axes(i) = ax;
colors = panel_options.Colors;
markers = panel_options.Markers;
legend_handles = gobjects(1,4);
for g = 1:4
    model_group = mouse_fit_results(4).evaluation.groups(g);
    plot(model_group.time,model_group.volume_mm3,'Color',colors{g},'LineWidth',panel_options.CurveLineWidth);
    hold on
end
for g = 1:4
    source_group = mouse_fit_cases(4).groups(g);
    legend_handles(g) = scatter(source_group.time,source_group.mean,panel_options.MarkerSize,markers{g},'filled', ...
        'MarkerEdgeColor',colors{g},'MarkerFaceColor',colors{g});
    errorbar(source_group.time,source_group.mean,source_group.sd,'Color',colors{g}, ...
        'LineStyle','none','LineWidth',panel_options.CurveLineWidth,'CapSize',panel_options.CapSize);
end
lgd = legend(ax,legend_handles,panel_options.LegendLabels, ...
    'Location',panel_options.LegendLocation,'FontName',panel_options.FontName,'FontSize',panel_options.CT26LegendFontSize,'Box','on');
lgd.ItemTokenSize = panel_options.LegendTokenSize;
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize);
title(ax,panel_options.Title, ...
    'FontSize',panel_options.TitleFontSize,'FontWeight','normal');
set(ax,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize,'LineWidth',panel_options.LineWidth, ...
    'Box','on','TickDir','out','LooseInset',get(ax,'TightInset'));

if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end


%% 导出当前组图（只保存，不重画）
F = mouse_fit_handle;
export_resolution = 600; % 本次组图导出分辨率（dpi）
export_files = save_mouse_fit;
visible_state = mouse_fit_visible;
figure_dir = mouse_fit_output_dir;
output_stem = fullfile(figure_dir,'mouse_fit');
if export_files && ~exist(figure_dir,'dir'); mkdir(figure_dir); end
if export_files
    exportgraphics(F,[output_stem '.pdf'],'ContentType','vector','BackgroundColor','none');
    exportgraphics(F,[output_stem '.png'],'Resolution',export_resolution,'BackgroundColor','white');
end

if strcmp(visible_state,'off'); close(mouse_fit_handle); end
