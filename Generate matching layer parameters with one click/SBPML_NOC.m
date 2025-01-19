function [ Inf_NOC_Na, Inf_NOC_Nb,Inf_NOC_Nele,ele_first] = SBPML_NOC( NOC_pm,XofN_pm,rows,Set_b)
%% ��Ԫ����Ϣ
NNpm=size(XofN_pm,1);
NEpm=size(NOC_pm,1);
%% �����߶�Ӧ����
% % 12����
SidesInfor=[1,2;2,3;3,4;4,1];
Sides_No=zeros(NEpm*4,2); % ÿ����Ԫ��4���߽���洢��Sides_No
for jj=1:NEpm
    NOC_ele=NOC_pm(jj,1:4)';
    Snode=NOC_ele(SidesInfor(:,1),:);
    Enode=NOC_ele(SidesInfor(:,2),:);
    Sides_No((jj-1)*4+1:(jj)*4,:)=[Snode,Enode];
end
[Sides_No1,I] = sort(Sides_No,2,'ascend');
%  ÿ����Ԫ��ÿ���⡢С������ǰ��������ں�
[Sides_No2] = unique(Sides_No1(:,1:2),'rows');
%  ��Ԫ�乲���⣨�ߣ������ȥ���ظ���
%% ���������Ϣ
% ͨ���⣬�ҳ�ÿ��������ӵĽ�㣬�����6�����������
Infor_NIN=zeros(NNpm,4);
esp=1.0e-9;
for ii=1:NNpm
    [ro1,co]=find(abs(Sides_No2(:,1)-ii)<esp);
    [ro2,co]=find(abs(Sides_No2(:,2)-ii)<esp);
    % �ҳ��ý���������������ϵ��ס�β���
    XX=[Sides_No2(ro1,2); Sides_No2(ro2,1)].';
    co1=size(XX,2);
    Infor_NIN(ii,1:co1)=XX;
end
% Infor_NIN Ϊÿ����������ӵĽ�㣬�������ýڵ㱾��
%% ����Ϣ

% ����ÿһ���㣬Ȼ���ҳ���ͨ�Ե�Ԫ����ÿ����Ԫ���ӵ����ڲ㼰�������
%% �����ҵڶ���
Infor_boun=Infor_NIN(Set_b,:);             % �߽��������ӵĽ��
In_plane_boun=zeros(size(Set_b,1),3);   % �ý����ƽ����������ĸ����������  ��һ��Ϊ��㱾��
Out_plane_boun=zeros(size(Set_b,1),3);% �ý�����������������������
In_plane_boun(:,1)=Set_b;
Out_plane_boun(:,1)=Set_b;

for ii=1:size(Set_b,1)
    Layer1=setdiff(Infor_boun(ii,:),Set_b);    % ÿ���߽������ӵķǱ�����
    Layer1(find(Layer1==0))=[];
    co1=size(Layer1,2);
    Out_plane_boun(ii,2:2+co1-1)=Layer1;
    %     Layer0=intersect(Infor_boun(ii,:),Set_b);  % ÿ���߽������ӵı�����
    Layer0=setdiff(Infor_boun(ii,:),Layer1);  % ÿ���߽������ӵı�����
    Layer0(find(Layer0==0))=[];
    co2=size(Layer0,2);
    In_plane_boun(ii,2:2+co2-1)=Layer0;
end
In_boundary=In_plane_boun;
Set_2=Out_plane_boun(:,2);
%% �������ڲ�����
Set=Set_2;
Setboundary=Set_b;
In_layer_setdiff_inboundary=cell(1,rows);          % ��ÿ�������ӵ����ڽ��
radial_setdiff_setb2=zeros(size(Set_b,1),rows);  % ��ÿ�������ӵ��ڲ�����

for jj=1:rows-2
    In_plane_boun=zeros(size(Set,1),3);
    Out_plane_boun=zeros(size(Set,1),3);
    In_plane_boun(:,1)=Set;
    Out_plane_boun(:,1)=Set;
    Line_Layer_node=Infor_NIN(Set,:);
    for ii=1:size(Set,1)
        Layer1=setdiff(Line_Layer_node(ii,:),Set);
        Layer1(find(Layer1==0))=[];
        Layer2=setdiff(Layer1,[Setboundary]);  % �޳�������ƽ�������
        co1=size(Layer2,2);
        Out_plane_boun(ii,2:2+co1-1)=Layer2;
        
        Layer0=setdiff(Line_Layer_node(ii,:),Layer1);   % �޳�������ƽ����������ڲ����ⲿ
        Layer0(find(Layer0==0))=[];
        co2=size(Layer0,2);
        In_plane_boun(ii,2:2+co2-1)=Layer0;
    end
    In_layer_setdiff_inboundary{:,rows-jj}=In_plane_boun;
    radial_setdiff_setb2(:,rows-jj-1)=Out_plane_boun(:,2);
    Setboundary=Set;
    Set=Out_plane_boun(:,2);
end
%% �����ڲ�
Line_Layer_node=Infor_NIN(Set,:);
In_plane_boun=zeros(size(Set,1),3);
In_plane_boun(:,1)=Set;
for ii=1:size(Set,1)
    Layer0=setdiff(Line_Layer_node(ii,:),Setboundary);   % �ų����ڲ��
    Layer0(find(Layer0==0))=[];
    co2=size(Layer0,2);
    In_plane_boun(ii,2:2+co2-1)=Layer0;
end
In_layer_setdiff_inboundary{:,1}=In_plane_boun;

radial_setdiff_setb2(:,rows-1:rows)=[Set_2,Set_b];
In_layer_setdiff_inboundary{1,rows}=In_boundary;

radial=radial_setdiff_setb2;
In_layer=In_layer_setdiff_inboundary;
Num_first=size(radial,2);
% first_co0 = [3.145037;3.943812;4.945453;6.201686;...
%     7.776794;9.752131;12.229033;15.335184;
%     19.230184;24.114568;30.239497];
first_co=[0:1/(Num_first-1):1];
% for i = 1:length(first_co0)
%     first_co(i+1) = sum(first_co0(1:i))/sum(first_co0);
% end
Inf_NOC_Nb=zeros(NEpm,2);
Inf_NOC_Na=zeros(NEpm,2);
Inf_NOC_Nele=zeros(NEpm,4);
for ele=1:NEpm
    NOC_ele=NOC_pm(ele,1:4)';
    for ii=1:3
        [a0,b0,c0]=find(NOC_ele(ii,:)==radial);
        N_layer(ii)=b0;
    end
    ele_layer=max(N_layer);  % �õ�Ԫ��㴦�������㣬���� ele_layer Ϊ����ϴ�ֵ
    ele_first(ele,1:2)=first_co(:,ele_layer-1:ele_layer);
    In_N_layer=In_layer{1,ele_layer};   % �õ�Ԫ����������ӵĽ��
    NOC_Nb=intersect(NOC_ele,In_N_layer(:,1));
    Coo=XofN_pm(NOC_Nb.',:);  % ����������ӵĽ����ܲ�������ͨ��
% %     Coo=[Coo;Coo(1,:)];
%         Coo=[Coo];
%     XX=diff(Coo);
% %     
% %     a=zeros(4,3);
% %     for ii=1:4
% %         if ii<4
% %             a(ii,:)=cross(XX(ii,:),XX(ii+1,:));
% %         else
% %             a(ii,:)=cross(XX(ii,:),XX(1,:));
% %         end
% %     end
% %     %     a1= a(find(a~=0));
% %     [roo,coo]=find(abs(a)==max(max(abs(a))));
% %     a1= a(:,coo(1));
% %     if all(a1(:)<0)==0 && all(a1(:)>0)==0
% %         b=[2,3,4,1]';
% %         cc=[b,NOC_Nb(b,1),a1];
% %         [roa,b0,c]=find(cc(:,3)<0);
% %         value=cc(roa(1),2);
% %         cc(roa(1),2)= cc(roa(2),2);
% %         cc(roa(2),2)=value;
% %     else 
% %         b=[2,3,4,1]';
% %         cc=[b,NOC_Nb(b,1),a1];
% %     end
%     %     Inf_NOC_Nb(ele,:)
%     %     find(radial()=cc(:,2);
    [tf, index] = ismember(NOC_Nb(:,1), radial(:,ele_layer));
    
    Inf_NOC_Nele(ele,:)=[radial(index,ele_layer-1).',radial(index,ele_layer).'];
    Inf_NOC_Nb(ele,:)=radial(index,rows);
    Inf_NOC_Na(ele,:)=radial(index,1);
end
end

