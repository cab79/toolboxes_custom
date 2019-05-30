% analyzeBAT  analyze BAT experiment data from Presentation logfiles
%

%experimental params
tBeepStart = 5000; % approximate time of beep starting (we evaluate taps after this point only)

forceRecomputation = false;
doPlotIndividual = false; %plot tapping data for individuals

close all

%across subject vars
clear sub*

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up experimental subject information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%repeat tapping twice for each excerpt.
nTapRep = 2;
subjects = {'JRI2','DBE','JW','RA','DN',...
    'BF','JAG','WZ','RM','ATT',...
    'TA','DH','BVS','JLK','EW',...
    'SW','kd','nsd','slc','GDB',...
    'JCR','CRC','SRV','cv','raj',...
    'JRC','emi','dsh','fsj','srh'};

%how many metronome and free tapping trials run per subject--initially I did not
% use these additional tests.
nMetronome = zeros(size(subjects));
nFreetap = zeros(size(subjects));

%Later added metronome taps (3 tempo) and 2 freetapping runs
%-keyboard used...keyboard
if 0,
  subjects = {subjects{:}, 'LV-keyboard','KW-keyboard','JAMF','lf','ht' };
  nMetronome = [nMetronome 3 3 3 3 3];
  nFreetap = [nFreetap 2 2 2 2 2];
end

subsToDo = 1:length(subjects);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up excerpt lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excerptList = {'hth','owa','hsg','pan','ocj','sas','tju','kps','nyn','acl','sma','rrw'}';
%excerptList = {'kps'}
excerptGenre = ['rrrrjjjjoooo'];

% timing information for each excerpt
%   these are derived from tapping data
% CURRENTLY NOT USED
excerptIBI = [];
for iE = 1:size(excerptList,1),
  timing = batBeatTimes(excerptList{iE});
  excerptIBI(iE) = timing.IBI;
  excerptBeats{iE} = timing.t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process log file for each subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logdir = fullfile(G.paths.root, 'projects','bat','logs','');
savedir = fullfile(G.paths.root, 'projects','bat','results','');

%across subject variables
meanCV=[];
latency=[];
pctCorrect=[];
cor = [];
pp = [];

for iS = subsToDo,

  %historical artifact: Handle cases where we did not collect the metronome
  if nMetronome(iS) > 0,
    hasMetronomeTap = true;
  else
    hasMetronomeTap = false;
  end
  subject = subjects{iS};
  outfile = parseBATlogfile(subject, logdir, savedir, hasMetronomeTap);
  load(outfile); %loads three variables 'subject','trials','expDate'

  if doPlotIndividual,
    figure
    jisubplot(13,3,0,'tall',[],'fontsize',8)
    sts = ['\bf' subject repStr '\rm - ' expDate];
    jisuptitle(sts)
  end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% analyze excerpts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for iE = 1:size(excerptList,1),

        thisExcerpt = excerptList{iE};

        timing = batBeatTimes(thisExcerpt);

        %get all trials with this excerpt
        idx = strmatch(thisExcerpt,cellstr(char(trials(:).excerpt)),'exact');

        nITI = length(idx) - 3;
        require(nITI == nTapRep,'incorrect number of ITI trials')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ITI, Async timeseries
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        style = {'b:','r-','g-'};

        tstr = '';
        astr = '';
        meanITIs = [];
        CVs = [];
        latencys = [];
        sdAsyncs = [];
        
        
        if doPlotIndividual,
            nextplot
        end
        
        for i = 1:nITI,
            
            %no taps?
            if isempty(trials(idx(i)).data.tTap),
                meanITIs(i) = nan;
                CVs(i) = nan;
                latencys(i) = nan;
                sdAsyncs(i) = nan;
                subIBI(iS,iE,i) = nan;
                subTactusIBI(iS,iE,i) = nan;
                subITI(iS,iE,i) = nan;
                timing.tactus_IBI = timing.IBI;
                tstr = [tstr 'xx'];
                astr = [astr 'xx'];
                continue
            end
            
            tap = calc_tap(trials(idx(i)).data.tTap, timing.t, true); %adjust_tactus
            %% special case: offbeat taps
            % if a lot of relphases are near 0.5, need to double the target
            % rate to fill in offbeats and use these to calc async. Leave
            % rp as is in order to signal offbeats.
            tmprp = abs(tap.rp) / 0.5;
            ioff = (tmprp>0.9 & tmprp<1.1); %20% window around 0.5
            pctoff = sum(ioff)/length(ioff);
            tactus_multiplier = 1;
            if pctoff > 0.10, %arb: more than 10% are offbeats
                tactus_multiplier = 2;
                tap2 = calc_tap(trials(idx(i)).data.tTap, timing.t, true, tactus_multiplier); %doubled tactus
                tap.async = tap2.async;
            end
            
            timing.tactus_IBI = mean(diff(tap.tactus_target_times));
            
            if 0 && (tap.pct_missed > 5 || tactus_multiplier ~= 1),
               figure
               title([subject ', ' excerptList{iE} ', ' num2str(tap.pct_missed) '%, ' ...
                   num2str(tap.tapTargRatio) ', (' num2str(tactus_multiplier) ')'])
               gridx(timing.t,'k-',2)
               gridx(tap.tactus_target_times)
               gridx(trials(idx(i)).data.tTap,'r--')
               pause      
               close
            end
            
            %temp plot

%             z=exp(sqrt(-1)*2*pi*tap.rp);
%             figure
%             compass(z)
%             title([subjects{iS} ', ' excerptList{iE} '# ' num2str(i)]);
%             pause
%             close

            %collect results
            mITI = mean(tap.iti);
            sdITI = std(tap.iti);
            
            alltap{iS,iE,i} = trials(idx(i)).data.tTap;

            mAsync = mean(tap.async);
            sdAsync = std(tap.async);

            mRp = mean(tap.rp);
            
            meanITIs(i) = mITI;
            CVs(i) = sdITI / mITI;
            if ~isnan(meanITIs(i)),
                latencys(i) = tap.valid_t_down(1); %time to first tap
            else
                latencys(i) = nan;
            end
            sdAsyncs(i) = sdAsync;
            
            subTactusIBI(iS,iE,i) = timing.tactus_IBI;
            subITI(iS,iE,i) = mITI;
            subIBI(iS,iE,i) = timing.IBI;
            
            %plot
            if doPlotIndividual,
                %plot(trials(idx(i)).data.tITI/1000, trials(idx(i)).data.ITI, style{i})
                %hold on
                tstr = sprintf('%s%.1f (%.1f), ',...
                    tstr,mITI,sdITI);

                %plot iti
                plot(tap.t_iti/1000, tap.iti, style{i})
                hold on

                nextplot
                %asynchrony
                plot(tap.t_async/1000, tap.async, style{i})
                %plot(tap.valid_t_down/1000, tap.async_wrapped, style{i})
                %plot(tap.valid_t_down/1000, tap.rp, style{i})
                hold on

                astr = sprintf('%s%.1f, ',astr,sdAsync);

                nextplot('delta',[0 -1]); %back up to previous panel
            end
        end


        %finalize plots
        if doPlotIndividual,
            %ITI finalize
            xlim([0 19])
            if strcmp(thisExcerpt,'sma'),
                xlim([0 25])
            end

            %take overall (for plotting) as mean of ITIs
            if 0
                mITI = nanmean(meanITIs);
                boundyl = timing.tactus_IBI + (timing.tactus_IBI*[-.1 .1]);
                boundyl = max(0,boundyl);
                %keep plotted limits, unless they're extreme
                yl = ylim;
                yl(1) = min(boundyl(1),yl(1));
                yl(2) = min(boundyl(2),yl(2));
                if sum(yl)>0,
                    ylim(yl)
                end
            end

            gridy(timing.tactus_IBI)

            ullabel(thisExcerpt,'bold')
            tstr(end-1:end)=[]; %last ,
            tstr = sprintf('\\bf%.1f:\\rm %s',timing.tactus_IBI,tstr);
            title(tstr)

            %    if iE==1,
            %        hleg = legend('1st','2nd');
            %        set(hleg,'fontsize',8)
            %    end

            if ~currentplotis('atColumnEnd'),hideAxisLabels('x'),end

            %Async finalize
            nextplot
            %ITI finalize
            xlim([0 19])
            if strcmp(thisExcerpt,'sma'),
                xlim([0 25])
            end

            gridy
            if ~currentplotis('atColumnEnd'),hideAxisLabels('x'),end
            astr(end-1:end)=[]; %last ,
            title(astr)

        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% BAT test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %rt, correctness & confidence for answers
        %    rt = height of bar. green = correct, red = incorrect, alpha = confidence

        clear target answer conf
        rt = nan(1,5);
        correct = nan(1,5);
        conf = nan(1,5);

        for i = nITI+1:length(idx),
            d = trials(idx(i)).data;
            m = trials(idx(i)).manipulation;
            a = trials(idx(i)).amount;

            %generate column number: b0 b-10 b20 p-25 p25
            if (m=='b'), col = 1; else col = 3; end
            if (a<0), col = col + 1; elseif (a>0), col = col + 2; end

            target(col) = trials(idx(i)).target;
            answer(col) = trials(idx(i)).data.onBeatCode;
            conf(col) = trials(idx(i)).data.confidence;
            correct(col) = target(col)==answer(col);
            rt(col) = trials(idx(i)).data.rt;

        end %loop on perception trials

        %quantify this excerpt
        subMeanCV(iS,iE) = nanmean(CVs);
        subCorrect(iS,iE,:) = correct; %store all five columns
        subLatency(iS,iE) = nanmean(latencys);
        subSdAsync(iS,iE) = nanmean(sdAsyncs);
        subConfidence(iS,iE,:) = conf;


        %plot barplot
        if doPlotIndividual

            ymax = 10; %sec
            nextplot
            xlim([0 6])
            hold on

            for col = 1:5,
                if ~isnan(rt(col)),
                    if correct(col) == 1,
                        color = 'g';
                    else
                        color = 'r';
                    end
                    thisrt = rt(col)/1000;
                    if thisrt<0, %indicates no rt response
                        thisrt = -1;
                    end

                    h = bar(col,thisrt, 0.7, color);
                    fc = get(h,'facecolor');
                    fc = fc * (1 - (3-conf(col))/3);
                    set(h,'facecolor',fc)

                    if thisrt > ymax,
                        plot(col, ymax, 'k^','markersize',8); %indicate it's off scale
                    end

                end
                set(gca,'xtick',1:5,'xticklabel',{'b0','b-10','b10','p-25','p25'})
                ylabel('rt [s]')
                ylim([-1.5 ymax])
            end


        end %loop on seq type

        drawnow

    end %loop on excerpts

    %subject  mean (across excerpts)
    meanSdAsync(iS) = nanmean(subSdAsync(iS,:));
    meanCV(iS) = nanmean(subMeanCV(iS,:));
    meanLatency(iS) = nanmean(subLatency(iS,:));
    tmp = squeeze(subCorrect(iS,:,:));
    pctCorrect(iS) = nanmean(tmp(:)) * 100;
    
    pctCorrect_on(iS) = nanmean(tmp(:,1)) * 100;
    pctCorrect_phase(iS) = nanmean(nanmean(tmp(:,4:5))) * 100;
    pctCorrect_tempo(iS) = nanmean(nanmean(tmp(:,2:3))) * 100;
    
    xx = subTactusIBI(iS,:,:);
    yy = subITI(iS,:,:);
    xx=xx(:); yy=yy(:);
    [r,p]=corrcoef(xx(~isnan(xx+yy)),yy(~isnan(xx+yy)));
    r = r(1,2);
    p = p(1,2);
    cor(iS) = r;
    pp(iS) = p;
    
    if doPlotIndividual
        nextplot


        plot(xx,yy,'+')
        title('ITI vs IBI')

        nextplot
        hist(subMeanCV(iS,:))
        title('ITI CV');
        %xlim([0 0.1])

        nextplot
        hist(subSdAsync(iS,:))
        ullabel('sd Async')

        %printc

        sts = [sts ', ' num2str(pctCorrect(iS),2) '% cor., ITI rho=' num2str(cor(iS),2)];
        jisuptitle(sts)

        drawnow
        pause

    end

end %loop on subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tap examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSubj = length(subjects);
for iE = 1:size(excerptList,1),
    figure
    plotgrid(excerptBeats{iE}/1000,[0 nSubj+1]-0.5,[],[])
    jititle(excerptList{iE})
    %sort subjs for this excerpt
    %mITI = nanmean(subITI(:,iE,:),3);
    mITI = subITI(:,iE,1); %only first trial
    [jnk,order] = sort(mITI);
    %put in subject order
    %order = 1:nSubj;
    for i = 1:nSubj,
        iS = order(i);
        yl = i + 0.5 *[-1 1];
       plotgrid(alltap{iS,iE,1 }/1000,yl,[],[],'r',2,'-')
       %plotgrid(alltap{iS,iE,2 },yl,[],[],'b',1,'-')
       text(-1500/1000,i,num2str(mITI(iS),4),'fontsize',10)
    end
    set(gca,'ydir','reverse')
    set(gca,'ytick',[1:nSubj],'yticklabel',subjects(order))
    xlim([-2000 22000]/1000)
end
    
%note: 8/21, wow--tremendous variation in tapping styles, some change midway
%   several tapped waltz in 4? or slower, some tapped patterns
%   clearly current look at mean ITI is flawed, given the possibility of change
%   how to handle? need something on a tap by tap basis, to track changes in
%   tactus. Also, try relative phase analysis w/ varying multipliers (instead
%   of current ratio of ITIs
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create jmp data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jmptitle = 'BAT Test results';
% jmpcol = {'subject', 'pctCorrect','pctCorrect_on','pctCorrect_phase','pctCorrect_tempo',...
%     'tempo corr','mean CV','mean sd Async','mean Latency'};
% jmpdata = {subjects', pctCorrect', pctCorrect_on', pctCorrect_phase', pctCorrect_tempo',...
%     cor',meanCV', meanSdAsync', meanLatency'};
% openInJmp(jmpdata, jmptitle, jmpcol)

%% to do
%a plot of an excerpt w/ its beats w/ all tap times underneath
%three individual ITI plots: JRI, Joe, one with large slow spread



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot across subject results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
jisubplot(2,2,1,'portrait',[],'fontsize',9)

plot(pctCorrect,meanCV,'bo','markerfacecolor','b')
title(['BAT n=' num2str(length(subsToDo))])
xlabel('perception score [%]')
ylabel('synchronization CV')
ylim([0 .13])

% least squares fit


nextplot
xx=cat(2,subTactusIBI(:,:,1),subTactusIBI(:,:,2));
%xx=cat(2,subIBI(:,:,1),subIBI(:,:,2));
yy=cat(2,subITI(:,:,1),subITI(:,:,2));
plot(xx',yy','o','markerfacecolor','auto');
hold on
plot([0 2200],[0 2200],'k-')
%plot([0 2200],[0 4400],'k-')
%plot([0 2200],[0 1100],'k-')

%title('iti vs ibi')
jixlabel('Stimulus Tactus IBI [ms]')
jiylabel('Tapping ITI [ms]')
axis([0 2500 0 2500])

%overall correlation coeficient

xx = subTactusIBI(:);
yy = subITI(:);
[rAll,pAll]=corrcoef(xx(~isnan(xx+yy)),yy(~isnan(xx+yy)));
    rAll = rAll(1,2);
    pAll = pAll(1,2);

%figure out deviation from perfect tempo tracking
ssreg = nansum((xx-yy).^2,2);
ssdat = nansum(xx.^2,2);
r2 = 1 - ssreg./ssdat;
%d = nanmean(d,2);
%d = log(d)

nextplot

xlabel('perception score [%]')
ylabel('tempo tracking correlation')
axis([40 102 0.6 1.05])
gridy(1)
hold on
plot(pctCorrect,cor,'bo','markerfacecolor','b')
box on

nextplot
hold on
plot([0 1200],[0 1200],'k-')
plot(xx(1,:),yy(1,:),'bo','markerfacecolor','b')
plot(xx(7,:),yy(7,:),'rs')
axis([200 1200 200 1200])
box on


xlabel('Stimulus IBI [ms]')
ylabel('Tapping ITI [ms]')

%pct correct
nextplot
boxplot([pctCorrect; pctCorrect_on; pctCorrect_tempo; pctCorrect_phase]','notch','on');
ylim([0 109])
set(gca,'xtick',1:4,'xticklabel',{'All','On Beat', 'Tempo Error', 'Phase Error'})
xlabel('Stimulus Condition')
ylabel('Percent Correct')

%num correct
nextplot
boxplot((12/100)*[pctCorrect; pctCorrect_on; pctCorrect_tempo; pctCorrect_phase]','notch','on');
ylim([0 15])
set(gca,'xtick',1:4,'xticklabel',{'All/3','On Beat', 'Tempo Error', 'Phase Error'})
xlabel('Stimulus Condition')
ylabel('Number Correct')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% additional breakdowns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% per song
figure
jisubplot(6,2,0,'tall')

xx=cat(2,subTactusIBI(:,:,1),subTactusIBI(:,:,2));
yy=cat(2,subITI(:,:,1),subITI(:,:,2));

for i = 1:12,
    nextplot
    plot(xx(:,[i i+12]), yy(:,[i i+12]), 'bo','markerfacecolor','b')
    title(excerptList{i})
    hold on
    plot([0 2200],[0 2200],'k-')
    axis([0 2500 0 2500])
end

%% per person
figure
jisubplot(10,3,0,'tall')

xx=cat(2,subTactusIBI(:,:,1),subTactusIBI(:,:,2));
yy=cat(2,subITI(:,:,1),subITI(:,:,2));

for i = 1:size(xx,1),
    nextplot
    plot(xx(i,:), yy(i,:), 'bo','markerfacecolor','b')
    title(['subj ' num2str(i)])
    hold on
    plot([0 2200],[0 2200],'k-')
    axis([0 2500 0 2500])
end

%% distribution of CV by excerpt

grp = repmat(excerptList', nSubj,1);
grp = grp(:);

%order by increasing mean CV
[meanCV, rmeanCV] = grpstats(subMeanCV(:),grp,{'mean','robust_mean'});
[rmeanCV,order] = sort(rmeanCV); 

%reform the data in this order
xx=subMeanCV(:,order);
xx=xx(:);
grp = repmat(excerptList(order)', nSubj,1);
grp = grp(:);

%cluster into three categories
k = kmeans(rmeanCV, 3);
borders = find(diff(k)~=0); %assumes categories are contiguous
%however, this does not agree with visal inspection, which looks more like
borders = [5 10];

borderCV = [];
for i = 1:length(borders),
  borderCV(i) = mean(rmeanCV(borders(i)+[0 1]));
end

figure
jisuptitle('BAT 2.0 excerpt synchronization accuracy')
subplot(2,1,1)
grpplot(xx,grp,'se','robust_mean')
gridx(borders + 0.5)
ylabel('Accuracy (CV)')

% subplot(2,1,2)
% boxplot(xx,grp)
% gridx(borders + 0.5)
% ylabel('Accuracy (CV)')
% xlabel('Excerpt')

%% distribution of pctCorrect by excerpt

pctCorrectExcerpt = 100 * nanmean(subCorrect,3);
grp = repmat(excerptList', nSubj,1);
grp = grp(:);

%order by decreasing correct (easy to hard)
[meanPC, rmeanPC] = grpstats(pctCorrectExcerpt(:),grp,{'mean','robust_mean'});
[~,order] = sort(-rmeanPC); 
rmeanPC = rmeanPC(order);

%reform the data in this order
xx=pctCorrectExcerpt(:,order);
xx=xx(:);
grp = repmat(excerptList(order)', nSubj,1);
grp = grp(:);

%cluster into three categories
k = kmeans(rmeanPC, 3);
borders = find(diff(k)~=0); %assumes categories are contiguous
%however, this does not agree with visal inspection, which looks more like
%borders = [5 10];

borderPC = [];
for i = 1:length(borders),
  borderPC(i) = mean(meanPC(borders(i)+[0 1]));
end

figure
jisuptitle('BAT 2.0 beat perception accuracy')
subplot(2,1,1)
grpplot(xx,grp,'se','robust_mean')
gridx(borders + 0.5)
ylabel('Accuracy (% correct)')

% subplot(2,1,2)
% boxplot(xx,grp)
% gridx(borders + 0.5)
% ylabel('Accuracy (% correct)')
% xlabel('Excerpt')


%% finally, do kmeans to partition based on pctCorrect and CV

%recalc unsorted data
[meanCV, rmeanCV] = grpstats(subMeanCV(:),grp,{'mean','robust_mean'});
[meanPC, rmeanPC] = grpstats(pctCorrectExcerpt(:),grp,{'mean','robust_mean'});
data = [rmeanCV meanPC];

nClust = 3;
k = kmeans(standardize(data),nClust);
color = 'rkbm';
figure
for i = 1:nClust,
  idx = find(k==i);
  plot(data(idx,1), data(idx,2),'o','color',color(i),'markerfacecolor',color(i))
  text(data(idx,1), data(idx,2), excerptList(idx),...
    'VerticalAlignment','bottom','HorizontalAlignment','left','color',color(i),'fontweight','bold')
  hold on
end
egend({'hard','medium','easy'}) %this will change on each run
ylabel('BAT perception (% correct)')
xlabel('BAT synchronization accuracy (CV)')
title('BAT 2.0 excerpt clustered by perception & production accuracy')
xlim([.04 .09])
ylim([60 90])
gridx(borderCV)
gridy(abs(borderPC))





