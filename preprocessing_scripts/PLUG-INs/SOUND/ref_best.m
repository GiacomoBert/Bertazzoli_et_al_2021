function [data_bestC] = ref_best(data, bestC)
% This function re-references the data to the channel with the least noise.
%
% .........................................................................
% 20 March 2017 : Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

data_bestC = data - repmat(data(bestC,:),[size(data,1),1]);

end
