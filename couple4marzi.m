clc
clear
close all
%% constant
Pi=0;
Po=0;
Ti=40;
To=10;
r_inner=1;
r_outer=2;
v=0.3;
E=2*10^9;
%D=(E/(1-(v^2)))*[1 v;v 1];
D=(E/((1+v)*(1-2*v)))*[1-v v;v 1+v];
alpha=13*10^-12;
w1=((1+v)*(1-2*v))/(E*(1-v));
w2=((1+v)/(1-v));
w3=v/(1-v);
n=5;
%% Assemble
nodes = linspace(r_inner, r_outer, 2*n + 1);
 elements = zeros(n, 6);
 ele = zeros(n,3);
 I=1:4:(4*n-2);
 J=6:4:(4*n+2);
for i = 1:n
    ele(i, :) = [2*i-1, 2*i, 2*i+1];
 elements(i, :) = [I(i):J(i)]; %(1_LedtNode, 2_MiddleNode, 3_RightNode)
end
% Assemble global matrix for temperature problem
K_total1 = zeros(4*n+2);
F_total1=zeros(4*n+2,1);
%% temperuter equation
for e = 1:n
node_ids = ele(e, :);
node_idc = elements(e , :);
r1 = nodes(node_ids);
r_i = r1(1);
r_j = r1(3);
l = r_j - r_i; 
r_avg = (r_i + r_j)/2;
r= sym('r');
L1 = (r_j - r)/l;
L2 = (r - r_i)/l;
N1 = L1*(2*L1 -1);
N2 = 4*L1*L2;
N3 = L2*(2*L2 -1);
N = [N1 0 N2 0 N3 0;0 N1 0 N2 0 N3];
N(1,:)=0;
dn=diff(N,r);
a1=dn.'*dn;
a2=N.'*dn;
k1=int(a1,r,r_i,r_j);
k2=(1/r_avg)*(int(a2,r,r_i,r_j));
if r_i==r_inner
    k=k1;
elseif r_j==r_outer
    k=k1;
else
    k=k1-k2;
end
k=double(k); % Assemble into global matrix
K_total1(node_idc, node_idc) = K_total1(node_idc, node_idc) + k;
end

%% displacement equation
K_total2 = zeros(4*n+2);
F_total2=zeros(4*n+2,1);
for e = 1:n
node_ids = ele(e, :);
node_idc = elements(e , :);
r1 = nodes(node_ids);
r_i = r1(1);
r_j = r1(3);
l = r_j - r_i; 
r_avg = (r_i + r_j)/2;
r= sym('r');
L1 = (r_j - r)/l;
L2 = (r - r_i)/l;
N1 = L1*(2*L1 -1);
N2 = 4*L1*L2;
N3 = L2*(2*L2 -1);
N=[N1 0 N2 0 N3 0;0 N1 0 N2 0 N3];
NU=N;
NU(2,:)=0;
NT=[0 N1 0 N2 0 N3;0 0 0 0 0 0];
dNdr = diff(NU);
dNdt=diff(NT);
a1=dNdr.'*dNdr;
a2=NU.'*dNdr;
a3=NU.'*NU;
a1t=dNdr.'*dNdt;
a3t=NU.'*NT;
d1=int(a1,r,r_i,r_j);
d2=(1/r_avg)*int(a2,r,r_i,r_j);
d3=(1/(r_avg^2))*int(a3,r,r_i,r_j);
 

dt=alpha*w2*r_avg * int(a1t,r,r_i,r_j);

if r_i==r_inner
    r=r_i;
    Ni=subs(NU);
    a3t=Ni.'*NT;
    a2i=Ni.'*Ni;
    r=sym('r');
    dt3=alpha*w2 * int(a3t,r,r_i,r_j);
    d4=((1/r_avg)*w3)*int(a2i,r,r_i,r_j);
    ke=d1+d3-d2+d4-dt+dt3;
elseif r_j==r_outer
    r=r_j;
    Nj=subs(NU);
    a3t=Nj.'*NT;
    a2j=Nj.'*Nj;
    r=sym('r');
    dt3=alpha*w2 * int(a3t,r,r_i,r_j);
    d4=((1/r_avg)*w3)*int(a2j,r,r_i,r_j);
    ke=d1+d3-d2+d4-dt+dt3;
else
ke=d1-d2+d3-dt;
end
ke=double(ke);
K_total2(node_idc,node_idc)=K_total2(node_idc,node_idc)+ke;
end
F_total2(1)=w1*Pi;
F_total2(4*n+1)=w1*Po;
%% merge two equation
K_tot=K_total2+K_total1;
K=K_tot;
F=F_total2;
K(2,:)=0;
K(:,2)=0;
K(2,2)=1;
K(4*n+2,:)=0;
K(:,4*n+2)=0;
K(4*n+2,4*n+2)=1;
F(2)=Ti;
F(4*n+2)=To;
for w=setdiff(1:4*n+1,2)
    F(w)=F(w)-(K_tot(2,w)*F(2))-(K_tot(4*n+2,w)*F(4*n+2));
end
a=inv(K)*F;
%disp(a)
t=zeros(2*n+1,1);
for i=1:2*n+1
    t(i)=a(2*i);
end
figure;
plot(nodes,t,'r-o')
u=zeros(2*n+1,1);
for j=1:2*n+1
    u(j)=a(2*j-1);
end
%% stress analysis
ee=zeros(2,1);
s_rr=zeros(n,1);
s_oo=zeros(n,1);
for i = 1:n
node_ids = ele(i, :);
r1 = nodes(node_ids);
r_i = r1(1);
r_j = r1(3);
l = r_j - r_i;
r_avg=(r_j+r_i)/2;
r=sym('r');
L1 = (r_j - r)/l;
L2 = (r - r_i)/l;
N1 = L1*(2*L1 -1);
N2 = 4*L1*L2;
N3 = L2*(2*L2 -1);
N = [N1, N2, N3];
dn=diff(N);
nn=N/r;
B=[dn;nn];
r=r_i;
b=subs(B);
ee=b*u(2*i-1:2*i+1);
sigma=D*ee;
%sigma=sigma*10^-9;
sigma_rr=sigma(1);
sigma_oo=sigma(2);
s_rr(i)=s_rr(i)+sigma_rr;
s_oo(i)=s_oo(i)+sigma_oo;
end
x=linspace(r_inner,r_outer,n);
%disp(t)
%disp(s_rr)
%disp(s_oo)
%s_rr/10^8
disp(F)
figure;

plot(x,s_rr,'--b');
hold on
plot(x,s_oo,'--r');
grid on
legend('sigma_r','sigma_t')
