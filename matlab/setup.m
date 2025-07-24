function meta = setup(circ,meta)

N = size(meta.valuemap,1);
Ng = length(meta.varmap);
card = max(meta.valuemap,[],2);

nnodes = size(meta.valuemap,2);

meta.Li = cell(5,1);
for i=1:N,
    meta.Li{i} = [meta.varmap==i; ones(card(i),1)*meta.valuemap(i,1:Ng) == (1:card(i))'*ones(1,Ng)]';
    meta.Li{i}(:,2:card(i)+1) = (meta.Li{i}(:,1)*ones(1,card(i))).*meta.Li{i}(:,2:card(i)+1);
end

act = 1:Ng;
open = Ng+1:nnodes;

cnt = 1;
cnt2 = 0;
str = {'or','and'};


while ~isempty(open),

    addset = [];
    meta.liptf{cnt}  = [];
    M = zeros(nnodes,nnodes);
    for i=open,
        if strcmp(circ(i).fun,str{cnt2+1})&&(length(intersect(circ(i).children,act))==length(circ(i).children)),
            addset = [addset i];
            meta.liptf{cnt} = union(meta.liptf{cnt},circ(i).children)';
            M(i,circ(i).children) = 1;
        end
    end
    meta.levelf{cnt} = addset;
    meta.Lf{cnt} = M(addset,meta.liptf{cnt});
    act = union(act,addset);
    open = setdiff(open,addset);
    cnt2 = xor(cnt2,1);
    cnt = cnt+1;
end

Nlvl = cnt-1;

cnt2 = mod(Nlvl+1,2);

for l=Nlvl:-1:1,
    meta.liptb{l} = meta.levelf{l}+nnodes;
    meta.levelb{l} = [];
    %M = zeros(2*nnodes,2*nnodes);
    M = [];
    M3 = [];
   % idx = [];
    for i = meta.levelf{l},
        nc = length(circ(i).children);
        M2 = zeros(nc,2*nnodes);
        M2(:,i+nnodes) = 1;
        M4 = zeros(nnodes,nc);
        M4(circ(i).children,:) = eye(nc);
        %idx = [idx circ(i).children+nnodes];
        %M(circ(i).children+nnodes,i+nnodes) = 1;
        meta.levelb{l} = union(meta.levelb{l},circ(i).children+nnodes)';
        if cnt2,
            meta.liptb{l} = union(meta.liptb{l},circ(i).children)';
            for j=1:length(circ(i).children),
                M2(j,setdiff(circ(i).children,circ(i).children(j))) = 1;
               % M(j+nnodes,setdiff(circ(i).children,j)) = 1;
            end
        end
        M = [M;M2];
        M3 = [M3 M4];
    end
   % Lb{l} = M(levelb{l},liptb{l});
    meta.Lb{l} = M(:,meta.liptb{l});
    meta.Lb2{l} = M3(meta.levelb{l}-nnodes,:);
    if ~cnt2,
        meta.Lb{l} = meta.Lb2{l}*meta.Lb{l};
    end
    cnt2 = xor(cnt2,1);
end

meta.Nlvl = Nlvl;

start = 2*Nlvl*ones(2*nnodes,1);
start(1:Ng) = 1;
stop = zeros(2*nnodes,1);

for i=1:Nlvl,
    stop(meta.liptf{i}) = max(stop(meta.liptf{i}),i);
    stop(meta.liptb{i}) = max(stop(meta.liptb{i}),2*Nlvl+1-i);
    start(meta.levelf{i}) = min(start(meta.levelf{i}),i);
    start(meta.levelb{i}) = min(start(meta.levelb{i}),2*Nlvl+1-i);
end
stop(stop==0) = 2*Nlvl;
stop(1:Ng) = 2*Nlvl;

meta.activef = cell(Nlvl,1);
meta.activeb = cell(Nlvl,1);
for i=1:Nlvl,
    meta.activef{i} = setdiff(find((start<=i)&(stop>=i)),meta.levelf{i});
    meta.activeb{i} = setdiff(find((start<=2*Nlvl+1-i)&(stop>=2*Nlvl+1-i)),meta.levelb{i});
end

    
    