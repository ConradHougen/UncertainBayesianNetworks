function m_out = prob_forward(mc,m_in,card)

% m_in is the matrix of the messages of the parents
% a row for each setting of the parents

% mc is the matrix containing one pd in each column (one for each setting
% of the parents of a node X

% card is the cardinality of each parent variable

%N = length(m_in); % N is the number of parents of the node X 
%Ns = 2^N; % the number of conditional distributions at X

%%LMK - I do understand and commented it out
%pm_in=prod(m_in); %product along the rows gives a column of products,
%one for each setting of the parents


%mo = ones(1,Ns);


N = length(card);
if N==0,
    m_out = mc';
else
    Ns = prod(card);

    v = zeros(1,N);
    mo = ones(Ns,1);
    for i=0:Ns-1,
        for j=1:N,
            mo(i+1) = mo(i+1)*m_in(j,v(j)+1);
        end
        v(N) = v(N)+1;
        for j=N:-1:2,
            if v(j)==card(j),
                v(j) = 0;
                v(j-1) = v(j-1)+1;
            else
                break;
            end
        end
    end


    %for i=1:N,
        %ix = bitand(bitshift(0:Ns-1,i-N),1);
        %mo = mo.*((1-m_in(i)).^(1-ix)).*(m_in(i).^ix);
    %end

    %m_out = mc'*mo';

    m_out = mo'*mc';
    %m_out = mc'*pm_in; 
end






    
    