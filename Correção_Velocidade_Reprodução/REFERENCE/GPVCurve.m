function [averageCurve, weightedCurve] = GPVCurve(EstNote, nFrames, weightCoeff)

EstNoteTransp = EstNote;
EstNoteTranspNorm = EstNote;
trackSet = NaN(nFrames,length(EstNote));
trackSetNorm = trackSet;
trackMags = trackSet;

% visualising track set
for j = 1:length(EstNote)
    EstNoteTransp{j} = EstNote{j}';
    trackSet(EstNoteTransp{j}(:,1),j) = EstNoteTransp{j}(:,2);
end

% routine for calculating the global curve
for i = 1:length(EstNote)
    EstNote{i}(2,:) = (EstNote{i}(2,:) - mean(EstNote{i}(2,:)))/mean(EstNote{i}(2,:));  % normalise each note-let by the mean
    EstNoteTranspNorm{i} = EstNote{i}';                                                     % transposing to a column representation
    trackSetNorm(EstNoteTranspNorm{i}(:,1),i) = EstNoteTranspNorm{i}(:,2);                          % trackset
    trackMags(EstNoteTranspNorm{i}(:,1),i) = EstNoteTranspNorm{i}(:,3);
end

% average curve
avgNum = sum(~isnan(trackSetNorm),2);
avgNum (avgNum==0) = 1;

trackSetNorm(isnan(trackSetNorm)) = 0;
averageCurve = sum(trackSetNorm,2);

averageCurve = averageCurve ./ avgNum;

% weighted curve
weightedMags = trackMags .^ weightCoeff;
weightedSum = trackSetNorm .* weightedMags;
weightedSum (isnan(weightedSum)) = 0;
weightedCurve = sum(weightedSum,2);

weightedMags(isnan(weightedMags)) = 0;
wtdNum = sum(weightedMags,2);
wtdNum (wtdNum == 0) = 1;

weightedCurve = weightedCurve ./ wtdNum;



end