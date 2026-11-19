function plot_prcc_bars(export_files, visible_state)
% Supplementary Fig. S7 uses the union of the top five parameters per scenario.
arguments
    export_files (1,1) logical = false
    visible_state = 'on'
end
script_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(script_dir,'results');
load(fullfile(results_dir,'p01_prcc_results.mat'),'prcc_values');
load(fullfile(script_dir,'input_data','p01_prcc_design.mat'),'parameter_symbols','scenario_names');
draw_prcc_bars(prcc_values,parameter_symbols,scenario_names,5,[-12,40], ...
    fullfile(results_dir,'terminal_volume_prcc'),visible_state,export_files);
end

function selected = draw_prcc_bars(values,parameter_symbols,scenario_names,top_n,view_angle,output_stem,visible_state,export_files)
% Select the union by absolute PRCC within EACH scenario and retain the
% original best-absolute-rank order. N is not the final column count.
validateattributes(top_n,{'numeric'},{'scalar','integer','>=',1,'<=',numel(parameter_symbols)});
n_parameters=numel(parameter_symbols);
selected_mask=false(n_parameters,1);
best_rank=inf(n_parameters,1);
for k=1:5
    [~,order]=sort(abs(values(:,k)),'descend');
    selected_mask(order(1:top_n))=true;
    positions=zeros(n_parameters,1); positions(order)=1:n_parameters;
    best_rank=min(best_rank,positions);
end
selected=find(selected_mask);
[~,order]=sortrows([best_rank(selected),selected],[1,2]);
selected=selected(order);
% bar3 rows are treatment strategies; columns are the selected parameters.
display_order=[2,3,4,5,1]; % Place no treatment after the triple-combination row.
shown=values(selected,display_order)';
heights=abs(shown);
labels=parameter_symbols(selected);
% Match Tumor_response/plot_tumor_response.m (tumor_change_violin),
% keeping each label paired with its original scenario and data row.
treatment_labels=replace(scenario_names(display_order),["XELOX","anti-PD-1"],["Xel","aPD1"]);
n_parameters_shown=numel(selected);
figure_width=max(28,1.2*n_parameters_shown+14);
figure_height=16;
fig=figure('Color','w','Visible',visible_state,'Units','centimeters', ...
    'Position',[2,2,figure_width,figure_height],'Name','P01 PRCC: 3-D bars');
ax=axes(fig,'Position',[0.09,0.14,0.80,0.82]);
bars=bar3(ax,heights,0.44);
hold(ax,'on');
for p=1:numel(bars)
    % Heights are absolute; signed jet gradients run from green bases to
    % each original PRCC value at the tip, including blue negative tips.
    z=bars(p).ZData;
    y=bars(p).YData;
    colors=nan(size(z));
    for k=1:size(shown,1)
        row_mask=abs(y-k)<0.5;
        colors(row_mask)=sign(shown(k,p))*z(row_mask);
    end
    set(bars(p),'CData',colors,'FaceColor','interp', ...
        'EdgeColor',[0.25,0.25,0.25],'LineWidth',0.35);
end
% A light wire grid marks the zero plane.
[xzero,yzero]=meshgrid(0.5:1:n_parameters_shown+0.5,0.5:1:5.5);
mesh(ax,xzero,yzero,zeros(size(xzero)),'FaceColor','none', ...
    'EdgeColor',[0.60,0.60,0.60],'EdgeAlpha',0.40,'LineWidth',0.4);
% Tick coordinates are the centers of the bar columns and treatment rows.
colormap(ax,jet(256));
clim(ax,[-1,1]);
set(ax,'XLim',[0.4,n_parameters_shown+0.6],'YLim',[0.4,5.6], ...
    'ZLim',[0,1],'XTick',1:n_parameters_shown,'XTickLabel',cellstr(labels), ...
    'YTick',1:size(shown,1),'YTickLabel',cellstr(treatment_labels),'ZTick',0:0.25:1, ...
    'TickLabelInterpreter','latex','FontSize',11,'TickLength',[0.006,0.006], ...
    'LineWidth',0.7,'Box','off','Color','w', ...
    'GridColor',[0.65,0.65,0.65],'GridAlpha',0.22,'Projection','orthographic');
xlabel(ax,'Parameters','FontName','Arial','FontSize',12);
strategy_label=ylabel(ax,'Treatment strategy','FontName','Arial','FontSize',12);
zlabel(ax,'|PRCC|','FontName','Arial','FontSize',12);
grid(ax,'on');
view(ax,view_angle);
pbaspect(ax,[max(3.75,0.32*n_parameters_shown),2.2,1.05]);
drawnow;
% Manual treatment-axis title position, shifted left as requested.
set(strategy_label,'Units','normalized','Position',[0.02,-0.02,0], ...
    'HorizontalAlignment','center','Rotation',0);
cb=colorbar(ax,'Position',[0.92,0.31,0.018,0.43]);
cb.Ticks=-1:0.5:1; cb.Label.String='PRCC';
cb.FontName='Arial'; cb.FontSize=10;
if export_files
    resolution = 300;
    png_file = [output_stem '.png'];
    exportgraphics(fig,png_file,'Resolution',resolution,'BackgroundColor','white', ...
        'Units','centimeters','Padding',0.12);

    % Use the same tight image bounds in the PDF.
    rgb = imread(png_file);
    size_inches = [size(rgb,2),size(rgb,1)]/resolution;
    pdf_figure = figure('Color','w','Visible','off','Units','inches', ...
        'Position',[1,1,size_inches]);
    pdf_axes = axes(pdf_figure,'Position',[0,0,1,1]);
    image(pdf_axes,rgb);
    axis(pdf_axes,'image','off');
    exportgraphics(pdf_figure,[output_stem '.pdf'],'ContentType','image', ...
        'Resolution',resolution,'Padding','figure');
    close(pdf_figure);
end
if strcmp(visible_state,'off'); close(fig); end
end

