function  [LaR,FbR,jj4]=solve_eig_R(FeiR,ladR,re_degree,esp2)
%UNTITLED 此处显示有关此函数的摘要
%%%%%%%%%%%%%每一频率下的特征向量矩阵、和特征值矩阵
LaR=zeros(re_degree,re_degree);
Fb_R=zeros(2*re_degree,re_degree);
jj4=0;
for jj2=1:2*re_degree%%%%%%%找出模态中虚特征值小于0对应的的特征向量
    if    ((abs(imag(ladR(jj2,jj2)))<esp2) && (real((ladR(jj2,jj2)))>esp2 ))...
                     || (imag(ladR(jj2,jj2))<-esp2)
        jj4=1+jj4;
        LaR(jj4,jj4)=(ladR(jj2,jj2));
        Fb_R(:,jj4)=FeiR(:,jj2);%%%%%%%%%%找出特征值虚部小于零的特征向量
    end
end
FbR=Fb_R(1:re_degree,:);%
end

