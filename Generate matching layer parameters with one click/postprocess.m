function [ re_reslut ] = postprocess( data )
t=data(:,1);  
loading=data(:,2);
value_max=max(loading,[],1);
loading=loading./value_max;
chazhi=4;
ntime=length(data);
m=(ntime-1)/chazhi;
deltat=0.0001;
for i=1:m
    j=(i-1)*3+1:i*3+1;
    ti=t((i-1)*3+1):deltat:t(i*3+1); %待求未知因变量点
    t1=t(j);
        yt=loading(j);
        yti=lagrange(t1',yt',ti);%应用lagrange插值求得待求点
        number=length(yti);
        tt((i-1)*number+1:(i)*number)=ti;
        Earthquake((i-1)*number+1:(i)*number)=yti;     
end 
Earthq=Earthquake';
True_t=tt';
L0=1/400;esp=(L0*0.1);
t0=1;i=1;
max_ro=300;
TT=zeros(max_ro,1);ex=zeros(max_ro,1);
for I=1:length(Earthq)
    if abs(sqrt((True_t(I)-True_t(t0))^2+(Earthq(I)-Earthq(t0))^2)-L0)<esp
        i=i+1;
        TT(i,1)=True_t(I);
        ex(i,1)=Earthq(I);
        t0=I;
    end
end
re_reslut=[TT,ex.*value_max];
end

