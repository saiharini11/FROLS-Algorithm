clear all;
clc;
%Total no of data points
N=200;
%Initializing u(k) randomly between [-1,1]
u(1)=-1+2*rand;
u(2)=-1+2*rand;
y(1)=0;
y(2)=0;
%White gausian noise with 0 mean and 0.01 variance
%(0.2)^2*randn
%Generating output data
for k=3:200
    y(k)=-0.605*y(k-1)-0.163*(y(k-2)^2)+0.588*u(k-1)-0.240*u(k-2)+0.04*randn;
    u(k)=-1+2*rand;
    %length(y)
end 
%z=[u' y']
ybar=y';

ny=2;
nu=2;
ne=0;
l=3;
n=ny+nu+ne;

%Total no of terms
M=factorial(n+l)/(factorial(n)*factorial(l));
%Dictionary D
pm=zeros(198,M);
p0=1;
%l=1
for k=3:200
    p(k,:)=[1 y(k-1) u(k-1) y(k-2) u(k-2) y(k-1)^2 y(k-1)*u(k-1) y(k-1)*y(k-2) y(k-1)*u(k-2) u(k-1)^2 u(k-1)*y(k-2) u(k-1)*u(k-2) y(k-2)^2 y(k-2)*u(k-2) u(k-2)^2 (y(k-1)^2)*u(k-1) (y(k-1)^2)*y(k-2) (y(k-1)^2)*u(k-2) (u(k-1)^2)*y(k-1) (u(k-1)^2)*y(k-2) (u(k-1)^2)*u(k-2) (y(k-2)^2)*y(k-1) (y(k-2)^2)*u(k-1) (y(k-2)^2)*u(k-2) (u(k-2)^2)*y(k-1) (u(k-2)^2)*u(k-1) (u(k-2)^2)*y(k-2) y(k-1)^3 u(k-1)^3 y(k-2)^3 u(k-2)^3 y(k-1)*u(k-1)*y(k-2) y(k-1)*u(k-1)*u(k-2) u(k-1)*u(k-2)*y(k-2) y(k-1)*y(k-2)*u(k-2)];
end
%step 1
s=1;
q1=p;
size(q1);
size(ybar);
sigma=ybar'*ybar;
for m=1:M
    g(m)=(ybar'*q1(:,m))/(q1(:,m)'*q1(:,m));
    err1(m)=((g(m))^2*q1(m)'*q1(m))/sigma;  
end
 [val,idx]=max(err1);
 l1=idx;
 alpha1=p(:,l1);
 q1=p(:,l1);
 a11=1;
 a(1,1)=1;
 g1=g(l1);
 g=[g1];
 err(1)=err1(l1);
 alpha=[alpha1];
 q=[q1];
 l=[l1];
 size(q);
 %step s
 esr=1;
 for x=1:4
 %while esr>0.005
      s=s+1;
      for m=1:M
          count=0;
          for i=1:length(l)
            if(m==l(i))
                l(i);
                count=count+1;
                errs(m)=0;
            end
          end
          if(count==0)
              res=0;
            for r=1:s-1
                temp=((p(:,m)'*q(:,r))/(q(:,r)'*q(:,r)))*q(:,r);
                res=res+temp; 
            end
            qm(:,m)=p(:,m)-res;
            g(m)=(ybar'*qm(:,m))/(qm(:,m)'*qm(:,m));
            errs(m)=((g(m)^2)*qm(:,m)'*qm(:,m))/sigma;
          end
      end
      [val,idx]=max(errs);
      ls=idx;
      l=[l ls];
      alphas=p(:,ls);
      alpha=[alpha alphas];
      q=[q qm(:,ls)];
      for r=1:s-1
          a(r,s)=(q(:,r)'*p(:,ls))/(q(:,r)'*q(:,r));
      end
      a(s,s)=1;
      err(s)=errs(ls);
      sum=0;
      for k=1:s
          sum=sum+err(k);
      end
      esr=1-sum;
 end
  x=[];
  for i=1:length(l)
      x=[x p(:,l(i))];
  end
  gbar=alpha; 
  %x*beta=ybar
  %beta is found using least squares method
  %Sparse model -  Only necessary terms
  disp("Sparse model terms")
  l
  beta=inv(x'*x)*x'*ybar
  size(beta);
   
      
          
              
                

            
            

