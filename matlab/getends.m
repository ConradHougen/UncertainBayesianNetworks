function loc = getends(node)

%creates the starting list of active nodes - those that have one parent or
%child

%same for binary and multinomial

N = length(node);
loc = [];
for i=1:N,
    %if (length(node(i).children)==0)||(length(node(i).parents)==0),
    if length(node(i).children)+length(node(i).parents)==1,
        loc = [loc; i];
    end
end