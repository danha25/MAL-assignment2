%% Assignment MAL: Recording Speech, Fourier Transform, Cluster Analysis, and Speech recognition
%
%%
% The assignments  below should be solved and documented as a mini-project
% that will form the basis for the examination. When solving the exercises
% it is important that
%
%%
% * you document all relevant results and analyses that you have obtained/performed during the exercises
% * try to relate your results to the theoretical background of the methods being applied.
%
%%
% The documentation should be integrated (by adding new code/markdown cells)
% in this script which contains the exercises. This script is to be
% converted into a .pdf by using the markdown option in Matlab. Please see
% the following documentation for more information about markdown options
%%
% <https://se.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html>
%%
% The mini-project must be uploaded to Studienet as a single pdf-file.
% You can create a .pdf of your script/markdown by saving the script as an
% m file (e.g. 'MyFile.m') and the running the following code:
%
%%
%   publish('MyFile.m','pdf')
%%
% You can continuously check your output by clicking the "Publish" button
% under the 'Publish' tab. This will give you a quick html version of your document.
% Also, see the above documentation above on how to include images, graphs,
% tables, etc. in your file. If you have any questions about the exercises,
% you are strongly encouraged to have talk with your fellow students first.
% If you are reading the pdf version of this document, please note that you
% must find the script version and complete the assignments in the script.
%% 1. Recording Speech
% Using the MATLAB recorder that has been provided during the course, Record 100 speech signals of each of these
% words: Dog, Cat, Bird, Horse, Cow.

Fs=8000; %Sampling frequency
nBits=16;
nChannels=1;
time=3; %Recording time in seconds

word='bird'; % word you are recording
filePath='recordings/bird/%s%d.mat';
startIndex=1;   % start index of recording number 
repeat=25; % repeat recording


recObj = audiorecorder(Fs,nBits,nChannels);

recordingIndex=startIndex;
while recordingIndex<(startIndex+repeat)
    
    % Record audio signal
    disp('Start speaking.');
    recordblocking(recObj, time);
    disp('End of Recording.');

    % Play back the recording.
    play(recObj); % Listen if the recording is ok

    % Save the recording
    saveRecord = input('Do you want to save recording? Y/N [Y]: ', 's');
    if (strcmpi(saveRecord,'y') || isempty(saveRecord))
        % Store data in double-precision array.
        myRecording = getaudiodata(recObj);
        
        %Save recording to file. Filename format: "word""recording number".mat, ex. hello1.mat
        fileName=sprintf(filePath,word,recordingIndex);
        save(fileName,'myRecording') % Remember that the variable name will still be myRecording when you load the file
        disp(strcat(fileName, ' saved sucesfully'));
        recordingIndex=recordingIndex+1;
    end
end    
    

% You have to extract features and analyze these words using two different
% methods, and compare the performance of the two methods in the end.
% 
%% 2. Method 1
% *1. Calculating the FFT*
% clear; clc;
Fs = 8000;
wordSet = ["bird","cat","cow","dog","horse"];
dataSet = [];
for r = 1:5  
    for z = 1:100
        filePath = 'recordings/%s/%s%d.mat';
        fileName = sprintf(filePath, wordSet(r), wordSet(r), z);
        load(fileName);

        % * Calculate the FFT of your signals. 
        % * Compute the two-sided and one-sided spectrum of your recordings.
        f = fft(myRecording);
        L = length(f);

        psdTwoSided = abs(f/L);
        psdOneSided = psdTwoSided(1 : L/2 + 1);
        psdOneSided(2:end - 1) = psdOneSided(2:end - 1) * 2;

        freqLabel = Fs * (0:(L/2))/L; % Frequency labels on the x-axis

        %  figure() % plotting spectrums
        %  subplot(2,1,1)
        %  plot(psdTwoSided)
        %  title('Twosided')
        %  subplot(2,1,2)
        %  plot(psdOneSided)
        %  subplot(2,1,2)
        %  plot(freqLabel,psdOneSided)
        %  title('Onesided')
        
        % * Go through the one-sided spectrum in bins of 100Hz 
        % (0Hz – 100Hz, 101Hz-200Hz etc.), and find the most prominent peak 
        % (and index) in each of these frequency bins.
        % * From the most prominent peaks that you just extracted, take the 10
        % higest peak and the corresponding frequencies. The frequencies are the
        % features that we will use in the following cluster analysis.
        
        NyquistFreq = Fs/2;

        freqBinSize = 100; % frequency steps
        freqBin = round((length(psdOneSided) / NyquistFreq) * freqBinSize); %freqBinSize in samples
        k = 1;
        for i = 2:freqBin:length(psdOneSided) - freqBin
            [mag(k) ind(k)]= max(psdOneSided(i:i + freqBin));
            freq(k) = freqLabel(i + ind(k));
            k = k + 1;
        end
        totalmag = [];
        totalInd = [];
        for q = 1:10
            [mag2(q), ind2(q)] = max(mag);
            mag(ind2(q)) = 0;
            totalmag = [totalmag mag2(q)];
            totalInd = [totalInd ind2(q)];
        end
        for t = 1:length(mag2)
            freqPeaks(t) = freq(ind2(t));
        end
        freqPeaks = sort(freqPeaks);
        dataSet = [dataSet; freqPeaks];
    end
end

% *2. Cluster Analysis*
%
% * Use either K-means or hierarchical cluster anaysis to cluster the data.
% * Pick the parameters that you use for the cluster analysis and explain your choices.
% * Show the "elbow" plot for the cluster analysis.
% * Determine the specificity and sensitivity of the clustering.
K = 10;
kmeansWSS(K, dataSet)
K = 5;
rng(1);
opts = statset('Display', 'final');
[idx, c] = kmeans(dataSet, K, 'Distance', 'cityblock', 'Replicates', 200, 'Options', opts, 'Start', 'cluster');

sensitivity1 = [];
sensitivity1 = [sensitivity1 length(find(idx(1:100, :) == 1)) / 100];
sensitivity1 = [sensitivity1 length(find(idx(100:200, :) == 2)) / 100];
sensitivity1 = [sensitivity1 length(find(idx(200:300, :) == 3)) / 100];
sensitivity1 = [sensitivity1 length(find(idx(300:400, :) == 4)) / 100];
sensitivity1 = [sensitivity1 length(find(idx(400:500, :) == 5)) / 100]

specificity1 = [];
specificity1 = [specificity1 length(find(idx([100:500], :) ~= 1))/400];
specificity1 = [specificity1 length(find(idx([1:100 200:500], :) ~= 2))/400];
specificity1 = [specificity1 length(find(idx([1:200 300:500], :) ~= 3))/400];
specificity1 = [specificity1 length(find(idx([1:300 400:500], :) ~= 4))/400];
specificity1 = [specificity1 length(find(idx([1:400], :) ~= 5))/400]

% figure;
% plot3(dataSet(idx==1,1), dataSet(idx==1,2),dataSet(idx==1,3), 'b.', 'MarkerSize',12)
% hold on
% plot3(dataSet(idx==2,1), dataSet(idx==2,2),dataSet(idx==2,3), 'r.', 'MarkerSize',12)
% plot3(dataSet(idx==3,1), dataSet(idx==3,2),dataSet(idx==3,3), 'y.', 'MarkerSize',12)
% plot3(dataSet(idx==4,1), dataSet(idx==4,2),dataSet(idx==4,3), 'g.', 'MarkerSize',12)
% plot3(dataSet(idx==5,1), dataSet(idx==5,2),dataSet(idx==5,3), 'c.', 'MarkerSize',12)
% plot3(c(:,1), c(:,2),c(:,3),'kx','MarkerSize',15, 'LineWidth', 3) 
% legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4', 'Cluster 5','Centroids','Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off

% Silhouette plot
figure
[silh3, ~] = silhouette(dataSet, idx, 'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value'
ylabel 'Cluster'

%% 3. Method 2

% *1. Calculating the MFCC*
% * Calculate the MFCCs of your signals. Extract 10 feature values for each.

% *2. Model Generation*
% 
% * Use your different arrays of features to generate a GMM of each word.
% * Determine the specificity and sensivity of the models, by comparing each
% recording against all models in a posterior analysis.

% clc;
wordSet = ["bird", "cat", "cow", "dog", "horse"];
wordMFCC = zeros(500, 10);

s = 1;
for r = 1:500
    filePath = 'recordings/%s/%s%d.mat';
    fileName = sprintf(filePath, wordSet(s), wordSet(s), (mod(r, 100)== 0) * 100 + (mod(r, 100)~= 0) * mod(r, 100)); % ternary op, if mod(r,100) is 0 we are either at 100, 200, 300, 400 or 500, we get recording 100 for the respective word, else we get mod(r, 100) = r
    load(fileName);
    wordMFCC(r, :) = mfcc(myRecording, 8000, 'WindowLength', length(myRecording), 'LogEnergy', 'Ignore', 'NumCoeffs', 10);
    if (mod(r, 100) == 0)       % r is 100, 200, 300, 400, 500
        s = s + 1;              % counter for switching words from wordSet
        mfccArr = wordMFCC(r - 99:r);
        GMMs = fitgmdist(mfccArr(:), 2);
        save(sprintf('GMMs/gmm%s.mat', wordSet(r / 100)), 'GMMs');
    end
end


%%
% clc;

posteriorRes = zeros(500, 5); % [nLogLBird, nLogLCat, nLogLCow, nLogLCow, nLogLHorse] mxn(m is no. of rec; n is likelihood which word it is)

for v = 1:5
    filename = sprintf('GMMs/gmm%s.mat', wordSet(v));
    varname = sprintf('gmm%s', wordSet(v));
    load(filename);
    for k = 1:500
        [~, nLogLWord] = posterior(GMMs, wordMFCC(k, :)');
        posteriorRes(k, v) = nLogLWord;
    end
end

minVals = [];
sensitivity2 = [];
specificity2 = [];
for k5=0:4
    minValsHelp = [];
    for k3 = 1:100
        [m i] = min(posteriorRes((100*k5)+k3, :));
        minValsHelp = [minValsHelp; m i];
    end
    minVals = [minVals; minValsHelp];
    countFirstWord = length(find(minValsHelp(:, 2) == k5+1));
    sensitivity2 = [sensitivity2 countFirstWord/100];
end

specificity2 = [specificity2 length(find(minVals(100:500, 2) ~= 1))/400];
specificity2 = [specificity2 length(find(minVals([1:100 200:500], 2) ~= 2))/400];
specificity2 = [specificity2 length(find(minVals([1:200 300:500], 2) ~= 3))/400];
specificity2 = [specificity2 length(find(minVals([1:300 400:500], 2) ~= 4))/400];
specificity2 = [specificity2 length(find(minVals([1:400], 2) ~= 5))/400]
sensitivity2

%% 4. Comparing the methods
% Use the calculated specificity and sensitivity for each method to compare
% their performance - which method performs the best? Describe why you
% think that this method is better than the other. List the pros and cons
% for each of these methods when doing speech recognition.
clc;

% Judging from the observations below, Kmeans proves to be more accurate.
sensitivity1
sensitivity2

specificity1
specificity2