    function [ Ah,Ch,Ln ] =Continue_fraction_Do(mm,CsH,GH)
    ele_Co1=[1,1;1,1];
    ele_Co2=[1,0;0,1];
    ele_Co3=[2,1;1,2];

    
    Co1=zeros(mm+1);
    Co2=zeros(mm+1);
    Co3=zeros(mm+1);
    if mm~=0
        Noc=zeros(mm,2);
        Noc(:,1)=1:mm;
        Noc(:,2)=2:mm+1;
        for n = 1:mm
            degree = Noc(n,:);
            Co1(degree,degree) =  Co1(degree,degree) + ele_Co1;
            Co2(degree,degree) =  Co2(degree,degree) + ele_Co2;
            Co3(degree,degree) =  Co3(degree,degree) + ele_Co3;
        end
        Co1(end,end)=2;
        Co2(end,end)=2;
    else
        Co1=1;
        Co2=1;
        Co3=3;

        
    end
    Ahs=GH*CsH/2;
    Chs=GH/CsH;
    Ah=Ahs*Co1;
    Ch=Chs*Co2;
    
    Ln=1/6*Co3;
    end

