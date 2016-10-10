function letter_discrim_behav_perf
allLetters=['EIXT';'OUKV';'SDAZ';'LNYH'];
% dirName=cd;
% load([dirName,'\test\',date,'_perf.mat'])
load('C:\Users\Xing\Lick\test\05-Oct-2016_perf - Copy.mat')
indCorr=find(performance==1);
indIncorr=find(performance==-1);
numCorrTrials=length(indCorr);
numIncorrTrials=length(indIncorr);
meanPerfAll=numCorrTrials/(numCorrTrials+numIncorrTrials);%average performance across conditions

for i=1:size(allLetters,1)
    for j=1:size(allLetters,2)
        trialsMatchingCondition=find(char(allTargetLetters)==allLetters(i,j));
        trialsMatchingCondition=trialsMatchingCondition(trialsMatchingCondition<=length(performance));
        trialCorrCond{i,j}=trialsMatchingCondition(performance(trialsMatchingCondition)==1);%trial numbers where correct response was made, for each condition
        trialIncorrCond{i,j}=trialsMatchingCondition(performance(trialsMatchingCondition)==-1);%trial numbers where incorrect response was made, for each condition
        numCorrTrialsCond(i,j)=length(trialCorrCond{i,j});%number of trials where correct response was made, for each condition
        numIncorrTrialsCond(i,j)=length(trialIncorrCond{i,j});%number of trials where incorrect response was made, for each condition
        meanPerfCond(i,j)=numCorrTrialsCond(i,j)/(numCorrTrialsCond(i,j)+numIncorrTrialsCond(i,j));%performance for each condition
        incorrResponses{i,j}=behavResponse(trialIncorrCond{i,j});%record down the responses made during incorrect trials
    end
end

%calculate mean performance for each location:
%1) average across all trials at each location
numCorrTrialsLoc1=zeros(size(allLetters,1),1);
numIncorrTrialsLoc1=zeros(size(allLetters,1),1);
for i=1:size(allLetters,1)
    for j=1:size(allLetters,2)
        numCorrTrialsLoc1(i)=numCorrTrialsLoc1(i)+numCorrTrialsCond(i,j);
        numIncorrTrialsLoc1(i)=numIncorrTrialsLoc1(i)+numIncorrTrialsCond(i,j);
    end
    meanPerfCond1(i,1)=numCorrTrialsLoc1(i)/(numCorrTrialsLoc1(i)+numIncorrTrialsLoc1(i));
end
%2) average across mean for each condition at each location
meanPerfCond2=mean(meanPerfCond,2);

%analyse effects of sample location:
fig1=figure;
fig2=figure;
for i=1:size(allLetters,1)
    for j=1:size(allLetters,2)
        locXcorr{i,j}=allSampleX(trialCorrCond{i,j});%x-coordinates for correct trials
        locXincorr{i,j}=allSampleX(trialIncorrCond{i,j});%x-coordinates for incorrect trials
        locYcorr{i,j}=allSampleY(trialCorrCond{i,j});%y-coordinates for correct trials
        locYincorr{i,j}=allSampleY(trialIncorrCond{i,j});%y-coordinates for incorrect trials
        resCorr{i,j}=allVisualHeightResolution(trialCorrCond{i,j});%number of subdivisions in the image height for correct trials
        resIncorr{i,j}=allVisualHeightResolution(trialIncorrCond{i,j});%number of subdivisions in the image height for incorrect trials
        figure(fig1)
        plot(locXcorr{i,j},locYcorr{i,j},'ko');hold on
        plot(locXincorr{i,j},locYincorr{i,j},'rx');
        figure(fig2)
        plot(locXcorr{i,j},locYcorr{i,j},'o','MarkerEdgeColor',[0.25*i 0.5 0.25*j]);hold on
        plot(locXincorr{i,j},locYincorr{i,j},'x','MarkerEdgeColor',[0.25*i 0.5 0.25*j]);
    end
end
%analyse effects of simulated phosphene resolution:
for i=1:size(allLetters,1)
    for j=1:size(allLetters,2)
        minResCorr(i,j)=min(resCorr{i,j});
        minResIncorr(i,j)=min(resIncorr{i,j});
        maxResCorr(i,j)=max(resCorr{i,j});
        maxResIncorr(i,j)=max(resIncorr{i,j});
    end
end
minRes=min(minResCorr(:));
maxRes=max(maxResCorr(:));
tallyCorr=zeros(1,maxRes-minRes+1);
tallyIncorr=zeros(1,maxRes-minRes+1);
tallyCorrCond=zeros(size(allLetters,1),size(allLetters,2),maxRes-minRes+1);
tallyIncorrCond=zeros(size(allLetters,1),size(allLetters,2),maxRes-minRes+1);
for k=1:maxRes-minRes+1
    resolutions=minRes:maxRes;
    for i=1:size(allLetters,1)
        for j=1:size(allLetters,2)
            tallyCorr(k)=tallyCorr(k)+sum(resCorr{i,j}==resolutions(k));
            tallyIncorr(k)=tallyIncorr(k)+sum(resIncorr{i,j}==resolutions(k));
            tallyCorrCond(i,j,k)=tallyCorrCond(i,j,k)+sum(resCorr{i,j}==resolutions(k));
            tallyIncorrCond(i,j,k)=tallyIncorrCond(i,j,k)+sum(resIncorr{i,j}==resolutions(k));
        end
    end
    perfRes(k)=tallyCorr(k)/(tallyCorr(k)+tallyIncorr(k));%calculate average performance across conditions, for each resolution. Equal weight given to each condition regardless of number of trials for each condition
end
fig3=figure;
plot(perfRes);
perfResCond=tallyCorrCond./(tallyCorrCond+tallyIncorrCond);%mean perf for each resolution, for each condition
numConds=size(allLetters,1)*size(allLetters,2);
fig4=figure;
for k=1:maxRes-minRes+1
    formattedPerfResCond=perfResCond(:,:,k);
    formattedPerfResCond=formattedPerfResCond(:);
    meanPerfRes(k)=mean(formattedPerfResCond);
    stdPerfRes(k)=std(formattedPerfResCond);
    errorbar(resolutions(k),meanPerfRes(k),stdPerfRes(k));hold on
end
fig5=figure;
for i=1:size(allLetters,1)
    for j=1:size(allLetters,2)
        perfCond=[];
        for k=1:maxRes-minRes+1
            perfCond=[perfCond perfResCond(i,j,k)];
        end
        plot(minRes:maxRes,perfCond,'MarkerEdgeColor',[i/numConds i/numConds j/numConds]);hold on
    end
end

