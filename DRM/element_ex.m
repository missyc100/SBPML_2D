function index = element_ex(ele,index)
n = length(index);
for i = 1:n
    m = find(ele(:,1) == index(i));
    index(i,2:5) = ele(m,2:5);
end