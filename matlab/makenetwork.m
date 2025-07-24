function node = makenetwork(filename)

%this function creates network (structure) 
%from a list of parent-children pairs
%same for binary and multinomial 

fp = fopen(filename,'r');

% opens the file under filename (to read)?


node(1).children = [];
node(1).parents = [];


index = fscanf(fp,'%d');
frewind(fp);

% reads its content row by row and saves it in a
% column named index

N = max(index);


         
node(N).children = [];



while ~feof(fp),
    index = fscanf(fp,'%d',2);
    node(index(1)).children = [node(index(1)).children index(2)];
    node(index(2)).parents = [node(index(2)).parents index(1)];
end

fclose(fp);

% node(1).children = [];
% node(1).parents = [];
% 
% index = fscanf(fp,'%d');
% frewind(fp);
% 
% N = max(index);
% 
% %node(N).children = [];
% 
% 
% while ~feof(fp),
%     index = fscanf(fp,'%d',2);
%     node(index(1)).children = [node(index(1)).children index(2)];
%     node(index(2)).parents = [node(index(2)).parents index(1)];
% end
% 
% % reads the index 2 by 2 elements and creates 
% % rows of parents and rows of children for each of the nodes.
% 
% 
% for i=1:N,
%     node(i).parents = sort(node(i).parents);
%     node(i).children = sort(node(i).children);
% end
% 
% % sorts in increasing order.
% % % 
%fclose(fp);
    