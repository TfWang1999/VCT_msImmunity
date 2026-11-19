%% CRC 拟合：参数设置
% 本脚本只绘图；数值由同目录的计算脚本生成。
% 首次先运行“参数设置”和“初始化数据和画布”，随后任意子图节可单独 Ctrl+Enter。
% 每个子图节开头列出全部实际绘图参数；修改后 Ctrl+Enter 仅更新该子图。
% SizeCM 为独立图窗宽高；PositionCM 为组图内绘图区 [左 下 宽 高]，单位 cm；字号 pt。
save_crc_fit = false;         % true：导出图片；false：只预览
visible_state = 'on';           % 'on'：显示图窗；'off'：后台绘制并关闭
crc_fit_canvas_cm = [18.3 10.2]; % 整幅组图画布；各绘图区尺寸在对应分节设置

export_files = save_crc_fit;
visible_state = validatestring(visible_state,{'on','off'});

crc_fit_visible = visible_state;


%% 读取数据并初始化画布
entry_file = mfilename('fullpath');
if isempty(entry_file) && usejava('desktop')
    entry_file = matlab.desktop.editor.getActiveFilename;
end
script_dir = fileparts(entry_file);
addpath(genpath(fullfile(fileparts(script_dir),'src')));
loaded = load(fullfile(script_dir,'results','crc_results.mat'),'results');
results = loaded.results;
cases = crc_cases();
figure_dir = fullfile(script_dir,'results');
crc_fit_output_dir = figure_dir;
if exist('crc_fit_handle','var')
    delete(crc_fit_handle(isgraphics(crc_fit_handle,'figure')));
end
crc_fit_results = results;
crc_fit_cases = cases;
crc_fit_handle = figure('Visible',visible_state,'Units','centimeters','Position',[2 2 crc_fit_canvas_cm]);
crc_fit_axes = gobjects(1,8);

%% 患者 1：CRC 拟合
i = 1;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [1.21596 5.37309 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.LegendFontSize = 6.7; % 图例文字字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.LegendPosition = [0.15 0.925 0.80 0.045]; % 图例在整幅图窗内的 [左 下 宽 高]，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串

panel_options.LegendLabels = {'Radiation','Chemotherapy','Immunotherapy','Simulation','±15% band','Real data'}; % 图例文字列表（无单位），按对应曲线或色块的绘制顺序排列

visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end
legend_handles = [handles.radiation,handles.chemotherapy,handles.immunotherapy,handles.model,handles.ci,handles.data];
lgd = legend(ax,legend_handles, ...
    panel_options.LegendLabels, ...
    'Orientation','horizontal','FontName',panel_options.FontName,'FontSize',panel_options.LegendFontSize, ...
    'Box','off','NumColumns',6);
set(lgd,'Units','normalized','Position',panel_options.LegendPosition);

%% 患者 2：CRC 拟合
i = 2;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [5.39843 5.37309 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 患者 3：CRC 拟合
i = 3;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [9.58089 5.37309 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 患者 4：CRC 拟合
i = 4;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [13.76336 5.37309 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 患者 5：CRC 拟合
i = 5;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [1.21596 0.81040 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 患者 6：CRC 拟合
i = 6;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [5.39843 0.81040 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 患者 7：CRC 拟合
i = 7;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [9.58089 0.81040 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 患者 8：CRC 拟合
i = 8;
panel_options = struct(); % 以下数值仅控制本子图，直接修改后运行本节
panel_options.PositionCM = [13.76336 0.81040 3.55660 3.69067]; % 组图内本绘图区的 [左 下 宽 高]，单位 cm；原点在整图左下角
panel_options.FontName = 'Arial'; % 字体名称（字符串，无单位），如 Arial；需为本机已安装字体
panel_options.TickFontSize = 6.5; % 坐标轴刻度文字大小，单位 pt（1 pt = 1/72 英寸）
panel_options.LabelFontSize = 7.2; % 横轴和纵轴标题的字号，单位 pt
panel_options.TitleFontSize = 8; % 子图标题字号，单位 pt
panel_options.AnnotationFontSize = 6.8; % R² 标注字号，单位 pt
panel_options.LineWidth = 0.55; % 坐标轴边框及引用此参数的图形线条粗细，单位 pt
panel_options.ModelLineWidth = 0.85; % 模型模拟曲线粗细，单位 pt
panel_options.MarkerSize = 2.4; % 观测点标记的线性大小，单位 pt（errorbar 的 MarkerSize）
panel_options.CapSize = 1.5; % 误差棒末端横帽的长度，单位 pt
panel_options.ModelColor = '#0072BD'; % 模型曲线颜色；十六进制颜色字符串（无单位）
panel_options.DataColor = [0.9 0.4 0.1]; % 观测数据标记边缘及误差棒颜色；RGB 三元组，每项 0–1（无单位）
panel_options.MarkerFaceColor = [1 0.6 0.2]; % 观测点内部填充颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandColor = [0.3010 0.7450 0.9330]; % 模型曲线阴影带颜色；RGB 三元组，每项 0–1（无单位）
panel_options.BandAlpha = 0.15; % 阴影带不透明度，无单位；0 完全透明，1 完全不透明
panel_options.TherapyColors = [1 0 0;0 0 1;0 1 0]; % 治疗条颜色，每行一个 RGB（0–1）；依次为放疗、化疗、免疫治疗
panel_options.TherapyBarHeight = 0.035; % 单条治疗标记高度占本例最大肿瘤体积的比例，无单位
panel_options.TherapyBarTop = 1.18; % 治疗条基准位置相对于本例最大肿瘤体积的倍数，无单位
panel_options.TherapyBarInset = 0.10; % 每个治疗日色条左右两端向内缩进的宽度，单位 天
panel_options.XTickStep = 5; % 时间轴相邻刻度的间隔，单位 天
panel_options.YPadding = 1.36; % 自动纵轴上限相对于本例最大肿瘤体积的倍数，无单位
panel_options.R2Position = [0.52 0.53]; % R² 标注在绘图区内的 [横向 纵向] 位置，归一化比例（无单位）
panel_options.XLimits = []; % 横轴 [最小值 最大值]，单位：天；[] 表示沿用自动范围
panel_options.YLimits = []; % 纵轴 [最小值 最大值]，单位：cm³；[] 表示沿用自动范围

panel_options.XLabel = 'Time (days)'; % 横轴标题文字（无单位）；若修改数据口径，请同步修改文字中的单位
panel_options.YLabel = 'Tumor volume (cm^3)'; % 纵轴标题文字（无单位）；元胞数组中的不同字符串显示为不同行
panel_options.Title = crc_fit_cases(i).id; % 当前子图标题文字（无单位）；可替换为自定义字符串


visible_state = crc_fit_visible;
export_files = save_crc_fit;
figure_dir = crc_fit_output_dir;
assert(isgraphics(crc_fit_handle,'figure'),'请先运行本脚本的初始化数据和画布分节。');
if isgraphics(crc_fit_axes(i)); delete(crc_fit_axes(i)); end
ax = axes(crc_fit_handle,'Units','centimeters','Position',panel_options.PositionCM);
crc_fit_axes(i) = ax;
exp_data = crc_fit_cases(i).experimental_data;
exp_data_name = panel_options.Title;
saved_trajectory = crc_fit_results(i).evaluation;
t_sim = saved_trajectory.time;
V_sim = saved_trajectory.volume;
t_sim = t_sim + 1;
hold(ax,'on');

handles.ci = fill(ax,[t_sim;flipud(t_sim)], ...
    [V_sim*0.85;flipud(V_sim*1.15)],panel_options.BandColor, ...
    'FaceAlpha',panel_options.BandAlpha,'EdgeColor','none');
handles.data = errorbar(ax,exp_data.Volume(:,1),exp_data.Volume(:,2), ...
    exp_data.Volume(:,2)-exp_data.Volume(:,3), ...
    exp_data.Volume(:,4)-exp_data.Volume(:,2),'o', ...
    'Color',panel_options.DataColor,'MarkerFaceColor',panel_options.MarkerFaceColor, ...
    'MarkerSize',panel_options.MarkerSize,'LineWidth',panel_options.LineWidth,'CapSize',panel_options.CapSize);
handles.model = plot(ax,t_sim,V_sim,'Color',panel_options.ModelColor,'LineWidth',panel_options.ModelLineWidth);

max_volume = max([V_sim(:); exp_data.Volume(:,4)]);
bar_height = max_volume * panel_options.TherapyBarHeight;
bar_top = max_volume * panel_options.TherapyBarTop;
bar_width = panel_options.TherapyBarInset;
for therapy_i = 1:size(exp_data.Therapy,1)
    t = exp_data.t(therapy_i)+1;
    therapy = exp_data.Therapy(therapy_i,:);
    if therapy(1)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top,bar_top,bar_top+bar_height,bar_top+bar_height], ...
            panel_options.TherapyColors(1,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(2)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*1.25,bar_top+bar_height*1.25, ...
            bar_top+bar_height*2.25,bar_top+bar_height*2.25], ...
            panel_options.TherapyColors(2,:),'EdgeColor','none','HandleVisibility','off');
    end
    if therapy(4)>0
        fill(ax,[t+bar_width,t+1-bar_width,t+1-bar_width,t+bar_width], ...
            [bar_top+bar_height*2.5,bar_top+bar_height*2.5, ...
            bar_top+bar_height*3.5,bar_top+bar_height*3.5], ...
            panel_options.TherapyColors(3,:),'EdgeColor','none','HandleVisibility','off');
    end
end

handles.radiation = patch(ax,nan,nan,panel_options.TherapyColors(1,:),'EdgeColor','none');
handles.chemotherapy = patch(ax,nan,nan,panel_options.TherapyColors(2,:),'EdgeColor','none');
handles.immunotherapy = patch(ax,nan,nan,panel_options.TherapyColors(3,:),'EdgeColor','none');

xlim(ax,[1 max(exp_data.t)+2]);
ylim(ax,[0 max_volume*panel_options.YPadding]);
xticks(ax,0:panel_options.XTickStep:max(exp_data.t)+1);

obs_days = exp_data.Volume(:,1);
y_predicted = interp1(t_sim,V_sim,obs_days,'linear');
y_actual = exp_data.Volume(:,2);
ss_total = sum((y_actual-mean(y_actual)).^2);
ss_residual = sum((y_actual-y_predicted).^2);
r_squared = 1-(ss_residual/max(ss_total,eps));

set(ax,'LineWidth',panel_options.LineWidth,'FontName',panel_options.FontName,'FontSize',panel_options.TickFontSize, ...
    'Box','on','TickDir','out','Layer','top');
xlabel(ax,panel_options.XLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
ylabel(ax,panel_options.YLabel,'FontSize',panel_options.LabelFontSize,'FontName',panel_options.FontName);
title(ax,exp_data_name,'FontSize',panel_options.TitleFontSize,'FontName',panel_options.FontName, ...
    'FontWeight','normal');
text(ax,panel_options.R2Position(1),panel_options.R2Position(2),sprintf('R^2 = %.3f',r_squared),'Units','normalized', ...
    'FontSize',panel_options.AnnotationFontSize,'FontName',panel_options.FontName,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');
if ~isempty(panel_options.XLimits); xlim(ax,panel_options.XLimits); end
if ~isempty(panel_options.YLimits); ylim(ax,panel_options.YLimits); end
if i <= 4
    xlabel(ax,'');
end
if ~ismember(i,[1 5])
    ylabel(ax,'');
end

%% 导出当前组图（只保存，不重画）
F = crc_fit_handle;
export_resolution = 600; % 本次组图导出分辨率（dpi）
export_files = save_crc_fit;
visible_state = crc_fit_visible;
figure_dir = crc_fit_output_dir;
output_stem = fullfile(figure_dir,'crc_fit');
if export_files && ~exist(figure_dir,'dir'); mkdir(figure_dir); end
if export_files
    exportgraphics(F,[output_stem '.pdf'],'ContentType','vector','BackgroundColor','none');
    exportgraphics(F,[output_stem '.png'],'Resolution',export_resolution,'BackgroundColor','white');
end

if strcmp(visible_state,'off'); close(crc_fit_handle); end
