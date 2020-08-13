function tonalscore = sts (v, epsilon)
% Computation of the specific tonal score function of a single frame or block according
% with "THE TONALNESS SPECTRUM: FEATURE-BASED ESTIMATION OF TONAL
% COMPONENTS, DAFx2013, Kraft and Zolzer and Lerch" 
% inputs:
%       v:          specific feature of interest (vector)
%       epsilon:    normalisation constant
% output:
%       tonalscore: specific tonal score of the feature
 
tonalscore = exp(-1*epsilon * v.^2);

end
