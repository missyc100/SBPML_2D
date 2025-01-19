function  [LaL,FbL,jj4]=solve_eig_L(FeiL,ladL,re_degree,esp2)
%UNTITLED 此处显示有关此函数的摘要
%%%%%%%%%%%%%每一频率下的特征向量矩阵、和特征值矩阵
LaL=zeros(re_degree,re_degree);
Fb_L=zeros(2*re_degree,re_degree);
jj4=0;
for jj2=1:2*re_degree%%%%%%%找出模态中虚特征值小于0对应的的特征向量
    if  ((abs(imag(ladL(jj2,jj2)))<esp2) && (real((ladL(jj2,jj2)))<-esp2 ))...
                     || (imag(ladL(jj2,jj2))>esp2)
        jj4=1+jj4;
        LaL(jj4,jj4)=(ladL(jj2,jj2));
        Fb_L(:,jj4)=FeiL(:,jj2);%%%%%%%%%%找出特征值虚部小于零的特征向量
    end
end
FbL=Fb_L(1:re_degree,:);%
end

