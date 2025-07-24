% clear;
% for i=1:10,
%     test_bp('store_net3_10.mat','inferMonteBP',10*i); 
% end

dbp = load('inferSO4_store_net3.mat');
mr = dbp.we*[1;0;0.5];
sr = 2./dbp.we(:,3);

for i=1:10,
    dm = load(sprintf('inferMonteBP_%d_store_net3.mat',10*i));
    m(:,i) = dm.we*[1;0;0.5];
    s(:,i) = 2./dm.we(:,3);
end


m = [mr m];
s = [sr s];
v = m.*(1-m)./(s+1);

rho = corr(v);

clear m s

dbp = load('inferSO4_store_net3_10.mat');
mr = dbp.we*[1;0;0.5];
sr = 2./dbp.we(:,3);

for i=1:10,
    dm = load(sprintf('inferMonteBP_%d_store_net3_10.mat',10*i));
    m(:,i) = dm.we*[1;0;0.5];
    s(:,i) = 2./dm.we(:,3);
end


m = [mr m];
s = [sr s];
v = m.*(1-m)./(s+1);

rho2 = corr(v);