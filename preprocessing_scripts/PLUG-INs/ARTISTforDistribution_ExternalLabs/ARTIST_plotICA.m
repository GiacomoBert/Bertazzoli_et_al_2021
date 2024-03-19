function [] = ARTIST_plotICA(EEGepoch,cfg,artcomp,runNum)
% FUNCTION TO PLOT ICA MAPS FOR EACH ICA RUN FROM ARTIST

cd([cfg.folderPath filesep 'ICA' num2str(runNum)]);

W=EEGepoch.icaweights*EEGepoch.icasphere; Wm=pinv(W); sourcesig=reshape(W*EEGepoch.data(:,:), size(W,1), size(EEGepoch.data,2), size(EEGepoch.data,3));

% PLOT ALL COMPONENTS IN EACH ROUND OF ICA
figure();
for jjj=1:size(EEGepoch.icaweights, 1)
    clf; colormap(jet)
    xlimm = [-150 350]; dat = find(EEGepoch.times>xlimm(1) & EEGepoch.times<xlimm(2));
    if ~all(all(isnan(squeeze(W(jjj, :, :)))))
        subplot(3, 1, 1); xx = get(gca, 'Position'); set(gca, 'Position', [xx(1) xx(2)+.05 xx(3) xx(4)]);
        erpimage(squeeze(sourcesig(jjj,dat,:)), [], EEGepoch.times(dat),'', 1, 0, 'erp', 1, 'cbar'); % erpstd
        title(['Comp.' int2str(jjj)]);
        % CAXIS BASED ON 0-500MS (BC RTMS AND LICI HAVE
        % PRE-STIM ARTIFACT THAT PUSHES RANGE
        caxis([min(squeeze(nanmean(sourcesig(jjj, EEGepoch.times>0, :), 3))) max(squeeze(nanmean(sourcesig(jjj, EEGepoch.times<400, :), 3)))])
        
        subplot(3, 1, 2);xx = get(gca, 'Position');set(gca, 'Position', [xx(1)-.05 xx(2)-.12 xx(3)*.5 xx(4)*2]);
        topoplot(Wm(:,jjj), EEGepoch.chanlocs, 'shrink', 'on', 'plotrad', 0.60, 'electrodes', 'numbers');
        
        subplot('Position', [xx(1)-.04 xx(2)-.35 xx(3)+.08 xx(4)]);
        vv = squeeze(sourcesig(jjj, :, :));vvmax = max(abs(vv(:, :)));
        vvmaxperc = find(vvmax>prctile(vvmax, 90));
        plot(vvmax); xl = get(gca, 'ylim'); xl(1) = [];
        for i = 1:floor(size(EEGepoch.data,3)/10),line([i*10 i*10], [0 ceil(max(vvmax)*1.25)], 'Color', [.6 .6 .6], 'LineStyle','--'),end
        box off, xlabel('Trial'), ylabel('Max uV in trial (0-200ms)'), xlim([0 size(sourcesig,3)+1]);
        
        % PLOT TRIALS
        subplot('Position', [xx(1)+.5 xx(2)-.07 xx(3)*.4 xx(4)*1.6]);
        mx = 0;
        for vv = 1:length(vvmaxperc)
            plot(EEGepoch.times, squeeze(sourcesig(jjj, :, vvmaxperc(vv)))+sum(mx)); hold all; box off;
            mx(vv) = max(squeeze(sourcesig(jjj, :, vvmaxperc(vv))))-min(squeeze(sourcesig(jjj, :, vvmaxperc(vv))))+5;
        end
        if mx(1)~=sum(mx)
            ylim([-1*mx(1) sum(mx)])
        end
        line([0 0], [-1*mx(1) sum(mx)], 'Color', 'k')
        xlim([EEGepoch.times(1) EEGepoch.times(end)+5])
        xlabel('Times (ms)');ylabel('Amp of large trials (uV)')        
        saveas(gcf, [cfg.fullCondName '_'  '_Comp' int2str(jjj)], 'jpg');
    end
end;close all

% PLOT ALL COMPONENTS FOR THIS ITERATION AND WHICH ONES WERE CONSIDERED BAD
close all, figure('units', 'normalized', 'outerposition', [0 0 .8 .8]); colormap(jet);
for jjj = 1:size(EEGepoch.icaweights, 1)
    if size(EEGepoch.data, 1)>40; scalar = 7; else scalar = 4;end
    subplot(ceil(size(EEGepoch.data, 1)/scalar),ceil(size(EEGepoch.data,1)/scalar), jjj);
    xx = get(gca, 'Position'); set(gca, 'Position', [xx(1) xx(2) xx(3)*1.2 xx(4)*1.2]);
    xlimm = [-200 500];
    dat = find(EEGepoch.times>xlimm(1) & EEGepoch.times<xlimm(2));
    if isreal(Wm(:, jjj))
        topoplot(Wm(:,jjj),EEGepoch.chanlocs, 'shrink', 'on');
        if sum(ismember((artcomp), jjj))>0
            titl = 'BAD';
        else
            titl = '';
        end
        title(titl);
    end
end
savefig_ARTIST([cfg.fullCondName  '_' 'Comp_ALL' num2str(size(EEGepoch.icaweights,1)) 'comp'], 15, 15, 250, '', 4, [10 8]*3.5);
close all

