function obs = getobs(net)

N = length(net);
obs = zeros(N,2);
for i=1:N,
    if ~isempty(net(i).value),
        obs(i,:) = [i net(i).value+1];
    end
end
obs = obs(obs(:,1)>0,:);