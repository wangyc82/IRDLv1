function [V]=codonComposition(Seq)
%%%% compute the frequence of the 20 amino acid ;
%%%% the order of the 20 AA: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y

n=length(Seq);

A=zeros(n,1);

m=1;
for i=1:4
    for j=1:4
        M(m)=10*i+j;
        m=m+1;
    end
end

l=1;
for i=1:4
    for j=1:4
        for k=1:4
            C(l)=100*i+10*j+k;
            l=l+1;
        end
    end
end

Class1={'A'};
Class2={'T'};
Class3={'C'};
Class4={'G'};

pp=ismember(Seq',Class1);
pp=find(pp==1);f0(1)=length(pp);
A(pp)=1;

pp=ismember(Seq',Class2);
pp=find(pp==1);f0(2)=length(pp);
A(pp)=2;

pp=ismember(Seq',Class3);
pp=find(pp==1);f0(3)=length(pp);
A(pp)=3;

pp=ismember(Seq',Class4);
pp=find(pp==1);f0(4)=length(pp);
A(pp)=4;

f=zeros(1,16);

for i=1:16
    for j=1:n-1
        if(M(i)==A(j)*10+A(j+1))
            f(i)=f(i)+1;
        end
    end
end

f1=zeros(1,64);

for i=1:64
    for j=1:n-2
        if(C(i)==A(j)*100+A(j+1)*10+A(j+2)) 
            f1(i)=f1(i)+1;
        end
    end
end
% V=[];
% V=(f-min(f))/max(f);
 V=[f0/sum(f0),f/sum(f),f1/sum(f1)];

% g=zeros(1,20);
% g(1)=f(61)+f(62)+f(63)+f(64);
% g(2)=f(57)+f(58)+f(59)+f(60);
% g(3)=f(53)+f(54)+f(55)+f(56);
% g(4)=f(50)+f(51);
% g(5)=f(49)+f(52);
% g(6)=f(45)+f(46)+f(47)+f(48)+f(13)+f(16);
% g(7)=f(41)+f(42)+f(43)+f(44);
% g(8)=f(37)+f(38)+f(39)+f(40)+f(21)+f(24);;
% g(9)=f(34)+f(35);
% g(10)=f(33)+f(36);
% g(11)=f(30)+f(31);
% g(12)=f(32);
% g(13)=f(25)+f(26)+f(27)+f(28)+f(14)+f(15);
% g(14)=f(22)+f(23);
% g(15)=f(18)+f(19);
% g(16)=f(9)+f(10)+f(11)+f(12);
% g(17)=f(5)+f(6)+f(7);
% g(18)=f(8);
% g(19)=f(2)+f(3);
% g(20)=f(1)+f(4);
% % 
% % % V=g;
% %  V=(g-min(g))/max(g);
% V=g/sum(g);
