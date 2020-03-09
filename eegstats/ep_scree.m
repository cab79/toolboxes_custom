function ep_scree(FactorResults,noiseResults,randResults)
% ep_scree(FactorResults,noiseResults,randResults) -
% Provides tests for deciding how many factors to retain.
%
%Input:
%  FactorResults    : Structured array with the PCA results.  See doPCA.
%  noiseResults     : Results of PCA of noise average.
%  randResults      : Results of PCA of random normally distributed noise.
%

%History
%  by Joseph Dien (6/8/09)
%  jdien07@mac.com
%
%  bugfix 9/5/09 JD
%  Location of windows appearing partly off screen on some systems due to offset from having Dock (on Mac) or Task Bar
%  (on PC) at the bottom of the screen.
%
%  bugfix 2/26/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 1999-2020  Joseph Dien
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windowHeight=500;

global EPscree
    
    screeFigure=findobj('Name', 'ScreeWindow');
    if ~isempty(screeFigure)
        close(screeFigure)
    end
    EPscree.handles.window = figure('Name', 'ScreeWindow', 'NumberTitle', 'off');
    colormap jet;
    
    EPscree.FactorResults=FactorResults;
    EPscree.noiseResults=noiseResults;
    EPscree.randResults=randResults;
    
    if isfield(randResults,'screeST')
        EPscree.theDataScree = mean(FactorResults.screeST,2);
        EPscree.theRandScree = mean(randResults.screeST,2);
        if ~isempty(noiseResults)
            EPscree.theNoiseScree = mean(noiseResults.screeST,2);
        end
    elseif isfield(randResults,'scree')
        EPscree.theDataScree = FactorResults.scree;
        EPscree.theRandScree = randResults.scree;
        if ~isempty(noiseResults)
            EPscree.theNoiseScree = noiseResults.scree;
        end
    else
        error('No scree for data.');
    end
    
    EPscree.theRandScree=EPscree.theRandScree*(sum(EPscree.theDataScree)/sum(EPscree.theRandScree)); %scale to same size as real data
    
    EPscree.screeLines(1,:)=EPscree.theDataScree;
    EPscree.screeLines(2,:)=EPscree.theRandScree;
    EPscree.legend{1}='data';
    EPscree.legend{2}='random';
    if ~isempty(noiseResults)
        EPscree.screeLines(3,:)=EPscree.theNoiseScree;
        EPscree.legend{3}='noise';
    end
    
    EPscree.cumPercent=zeros(length(EPscree.theDataScree),1);
    for i=1:length(EPscree.theDataScree)
        EPscree.cumPercent(i)=sum(EPscree.theDataScree(1:i))/sum(EPscree.theDataScree);
    end
    
    EPscree.scree.ymax=max(max(EPscree.screeLines));
    EPscree.scree.ymin=min(min(EPscree.screeLines));
    EPscree.scree.xmax=min(20,length(EPscree.theDataScree));
    EPscree.scree.xmin=1;
    
    EPscree.percent=95;
    

EPscree.percentScree=max(find((EPscree.percent/100) >= EPscree.cumPercent));

%scree test

%EPscree.handles.scree.frame=uicontrol('Style','frame','Position',[45 windowHeight-355 555 260]);

EPscree.handles.screeChart = axes('units','pixels','position',[150 windowHeight-300 400 200]);

EPscree.handles.scree.Lines.data = plot([EPscree.scree.xmin:EPscree.scree.xmax],EPscree.screeLines(:,EPscree.scree.xmin:EPscree.scree.xmax),'Marker','o');
axis([EPscree.scree.xmin EPscree.scree.xmax EPscree.scree.ymin EPscree.scree.ymax]);
legend(EPscree.legend);

EPscree.handles.scree.ymax = uicontrol('Style','edit','HorizontalAlignment','left',...
    'String',EPscree.scree.ymax,'Position',[50 windowHeight-120 60 20],...
    'Callback', 'ep_scree;');
EPscree.handles.scree.ymin = uicontrol('Style','edit','HorizontalAlignment','left',...
    'String',EPscree.scree.ymin,'Position',[50 windowHeight-300 60 20],...
    'Callback', 'ep_scree;');
EPscree.handles.scree.xmin = uicontrol('Style','edit','HorizontalAlignment','left',...
    'String',EPscree.scree.xmin,'Position',[150 windowHeight-350 60 20],...
    'Callback', 'ep_scree;');
EPscree.handles.scree.xmax = uicontrol('Style','edit','HorizontalAlignment','left',...
    'String',EPscree.scree.xmax,'Position',[480 windowHeight-350 60 20],...
    'Callback', 'ep_scree;');


uicontrol('Style','text','HorizontalAlignment','left',...
    'String','Minimum %age accounted for criterion (nothing to do with scree test above)','Position',[50 windowHeight-380 500 20]);

EPscree.handles.percent = uicontrol('Style','edit','HorizontalAlignment','left',...
    'String',EPscree.percent,'Position',[50 windowHeight-400 60 20],...
    'Callback', 'ep_scree;');

uicontrol('Style','text','HorizontalAlignment','left',...
    'String',sprintf('Factors to retain: %d',EPscree.percentScree),'Position',[150 windowHeight-400 200 20]);

EPscree.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done',...
    'Position', [20 windowHeight-450 80 35], 'Callback', ['global EPscree; close(EPscree.handles.window); clear EPscree; ep(''start''); ']);



