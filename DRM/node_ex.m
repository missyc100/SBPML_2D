function index = node_ex(node,index)
n = length(index);
for i = 1:n
    m = find(node(:,1) == index(i));
    index(i,2) = node(m,2);
    index(i,3) = node(m,3);
end