function output = movingstd(input,n)

% Calculates at each timepoint of 'input' the local standard deviation across 'n' timepoints.
% __________________________________
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024

sumin=filter(ones(n,1),1,input);
sumsq=filter(ones(n,1),1,input.^2);

output=sqrt((sumsq-(sumin.^2)/n)/(n-1));
output=output(n:end);

