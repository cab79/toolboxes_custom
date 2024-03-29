function g = gplot_timeseries(g,P)
% requires
    % P.x
    % P.y
    % P.cond
    % P.xy
    % P.fact_names
    % P.colours
    % P.poly
    % P.xlinesolid
    % P.xlinedashed
    % P.ptitle - individual plots
    % P.gtitle - whole figure

%construct gramm plot
%cond = cond(end:-1:1);
%y = y(end:-1:1);

if P.color_order==-1
    P.colours = fliplr(P.colours')';
    P.y = fliplr(P.y')';
    P.cond = fliplr(P.cond);
end

g(P.xy(1),P.xy(2))=gramm('x',P.x,'y',P.y,'color',P.cond,'size',P.condsize);
g(P.xy(1),P.xy(2)).set_names('x',P.xaxisname,'y',P.yaxisname,'color',P.fact_names{end},'column','','row','');
g(P.xy(1),P.xy(2)).set_order_options('color',P.color_order);
g(P.xy(1),P.xy(2)).set_color_options('map',P.colours,'chroma_range',[30 90]);
g(P.xy(1),P.xy(2)).set_point_options('base_size',3);
g(P.xy(1),P.xy(2)).set_title(P.ptitle);
if ~P.legend
    g(P.xy(1),P.xy(2)).no_legend();
end
if ~isempty(P.xlinedashed)
    g(P.xy(1),P.xy(2)).geom_vline('xintercept',P.xlinedashed,'style','k--','extent',4);
end
if ~isempty(P.xlinesolid)
    g(P.xy(1),P.xy(2)).geom_vline('xintercept',P.xlinesolid,'style','k');
end
if P.poly
    g(P.xy(1),P.xy(2)).geom_polygon('x',{[min(P.poly) max(P.poly)]},'color',[0.6 0.6 0.6]);
end
%range(cl,:) = [min(P.E_valz) max(P.E_valz)];

switch P.plottype 
    case 'stat_smooth'
        g(P.xy(1),P.xy(2)).stat_smooth();
        g(P.xy(1),P.xy(2)).set_line_options('base_size',0.5);
    case 'stat_summary'
        g(P.xy(1),P.xy(2)).stat_summary('type','ci','geom','area','setylim',true);
    case 'geom_point'
        g(P.xy(1),P.xy(2)).geom_point();

end
end
%tix=get(g(P.xy(1),P.xy(2)).facet_axes_handles,'ytick')
%g(P.xy(1),P.xy(2)).axe_property('yticklabel',num2str(tix,'%.1f'));

%Possibility to set color and fill by indices (using a column vector of
%integers. Colormap generated between 1 and max(vector))
%   g(2,1).geom_polygon('y',{[5 20];  [20 30];  [30 50]},'color',[1 ; 3;  2]);
