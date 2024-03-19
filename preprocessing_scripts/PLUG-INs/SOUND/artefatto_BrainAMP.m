function [datatmp]=artefatto_BrainAMP(EEG,triggerstmp,janelacut,swindow,metodo,n)
samples1=round((EEG.srate*janelacut(1))/1000);
samples2=round(EEG.srate*janelacut(2)/1000);
samplescut=samples1:1:samples2;
samplescopy=samples1-1:-1:samples1-length(samplescut);

trialscut=repmat(triggerstmp',1,size(samplescut,2))+repmat(samplescut,size(triggerstmp,2),1);
trialscopy=repmat(triggerstmp',1,size(samplescut,2))+repmat(samplescopy,size(triggerstmp,2),1);
datatmp=EEG.data;
datatmp(:,trialscut)=datatmp(:,trialscopy);

samples11=round((EEG.srate*swindow(1))/1000);
samples22=round(EEG.srate*swindow(2)/1000);
samplessmooth=samples11:1:samples22;
% borda1;
samplessmooth1=samples11:1:samples22+samples1;
% borda2;
samplessmooth2=samples11:1:samples22+samples2;
trialssmooth1=repmat(triggerstmp',1,size(samplessmooth1,2))+repmat(samplessmooth1,size(triggerstmp,2),1);
trialssmooth2=repmat(triggerstmp',1,size(samplessmooth2,2))+repmat(samplessmooth2,size(triggerstmp,2),1);

h=waitbar(0,'Smoothing Data...');
for i=1:1:size(datatmp,1)
    for j=1:1:size(trialssmooth1,1)
        datatmp(i,trialssmooth2(j,:))=smooth(datatmp(i,trialssmooth2(j,:)),n,metodo);
    end
    waitbar(i/size(datatmp,1),h)
end
delete(h)

clear i j h janelacut metodo n samples1 samples11 samples2 samples22...
    samplescopy samplescut samplessmooth samplessmooth1 samplessmooth2...
    swindow trialscopy trialscut trialssmooth1 trialssmooth2 triggerstmp