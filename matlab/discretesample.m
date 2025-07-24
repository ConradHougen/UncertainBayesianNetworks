function [X] = discretesample(p, NS)

%   discretesample returns a column of NS random samples from a discrete
%   distribution p
%   p is given as a column of numbers in [0,1] 


%l=length(p);
l = size(p,1);   %%LMK Change 10/12/2018

pcumulative = cumsum(p,1); % cumulative sum of the elements of p

if nargin<2,  %%LMK chage 10/12/2018
    NS = size(p,2);
end

RanCol=rand(NS,1); % a column of NS random numbers in [0,1]

RanColmatrix = repmat(RanCol', l, 1); % puts RanCol in a row and replicates it l times 
                                      % returning a l times NS matrix with
                                      % repeating rows
if nargin>1,                     %%Added my LMK 10/12/2018                  
    pmatrix = repmat(pcumulative,1,NS); % replicates pcumulative NS times returning 
                                    % l times NS matrix with repeating
                                    % columns
else
    pmatrix = pcumulative;
end

S=(RanColmatrix-pmatrix)<0; % each column of the difference matrix is a comparison with 
                            % one random number with the pcumulative; this
                            % determines where the random number falls in
                            % [0,1] thus giving a way to sample a value by
                            % the next two commands 
Sum=sum(S,1);
X=(l-Sum)'; % in correspondence with Sum, one of the values 0,1,...,l-1 is chosen as the sampled value



 
          

end

