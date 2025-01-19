function [AMs,ACs,AKs] = integrate_S(XofNS,NOCS,cx,lamds,Gs,rs)
NNs=length(XofNS);
NEs=NNs-1;
NOC_B=NOCS(:,2:3);
    %%%%%%形成内域总刚阵和总质量阵
    AMs= zeros(2*NNs,2*NNs);
    ACs =zeros(2*NNs,2*NNs);
    AKs =zeros(2*NNs,2*NNs);
    node=zeros(1,2);
    node1=zeros(1,2);
    for n = 1:NEs
        Ms=MMs(n,XofNS,NOCS,cx,lamds,Gs,rs);
        Ms = diag(sum(Ms,2));
        Cs=CCs(n,XofNS,NOCS,cx,lamds,Gs,rs);
        Ks=KKs(n,XofNS,NOCS,cx,lamds,Gs,rs);
        for j1 = 1 : 2;
            node(j1) = NOC_B(n,j1);
            for i1=1:2;
                node1(i1)=NOC_B(n,i1);
                AMs(2*node(j1)-1:2*node(j1),2*node1(i1)-1:2*node1(i1)) =  AMs(2*node(j1)-1:2*node(j1),2*node1(i1)-1:2*node1(i1)) + Ms(2*j1-1:2*j1,2*i1-1:2*i1);
                ACs(2*node(j1)-1:2*node(j1),2*node1(i1)-1:2*node1(i1)) =  ACs(2*node(j1)-1:2*node(j1),2*node1(i1)-1:2*node1(i1)) + Cs(2*j1-1:2*j1,2*i1-1:2*i1);
                AKs(2*node(j1)-1:2*node(j1),2*node1(i1)-1:2*node1(i1)) =  AKs(2*node(j1)-1:2*node(j1),2*node1(i1)-1:2*node1(i1)) + Ks(2*j1-1:2*j1,2*i1-1:2*i1);
            end
        end
    end
end