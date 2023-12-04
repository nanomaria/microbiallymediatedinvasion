
close all
close all

global qhp qcm qcp qhm mup mum d N alpha beta cm


mup = 0.1;
mum = 0.1;
alpha = 0.4;
beta = 0.4;
qcm = 1;
qhm = 1;
qcp = 1;
qhp = 5;
d = 2;
cm = 0;

Q = qhp*alpha/d - qcp*beta;
Q1 = qcm*beta-qhm*alpha/d;

pmax = 7;
mmax = 7;


%syms pv mv
figure(2);
for i=1:10
    N=i;
n1 = fimplicit(@(pv,mv)  pv.*mv*N*Q./(pv/d+N*mv) - mup.*pv^2, [0,pmax],'Color', [0.6 0.3 0],'Linestyle','--','Linewidth',2-0.15*i);
hold on
n2 = fimplicit(@(pv,mv) mv.*pv*Q1./(pv/d+N*mv) - mum.*mv^2, [0,pmax],'Color', [0.1 0.1 0.7],'Linestyle','--','Linewidth',2-0.15*i);
hold on
line([Q/mup, Q/mup],[0,mmax],'Color',[0 0 0], 'Linestyle',':','Linewidth',1.2)



end
set(gca,'fontsize',14)
xlabel('p')
ylabel('m')




N = 10;

alpha = repelem(0.4,N);
beta = repelem(0.4,N);


d = 2;
p0 = 0.1;
m0 = repelem(0,N);

Tfin = 100;

p_star = zeros(1,N);
m_star = zeros(1,N);

for i =1:N

m0(i) = 0.1;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y] = ode45(@basic_eq_v, 0:.1:Tfin, [p0, m0], options);
p_star(i) = Y(end,1);
m_star(i) = Y(end,2);

%figure(2)
%plot(T,Y(:,1),'g')
%hold on

%figure(1)
%plot(i,p_star(i),'*','Linewidth',1.5)
%hold on
%set(gca,'fontsize',14)

end




plot(p_star,m_star,'s','MarkerSize',8,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])


%m = 0:.05:mmax;
%p = 0:.05:pmax;

%n1_e = mup*p.^2./(d*Q-d*mup*p);
%n2_e = -p./(2*d) + 0.5*sqrt((p/d).^2 + 4*p*Q1/mum);


figure(2);
%plot(p,n1_e, 'k--')
%hold on
%plot(p,n2_e,'k--')
plot([0,0],[0,mmax], 'Color', [0.6 0.3 0],'Linestyle','--','Linewidth',2)
line([0,pmax],[0,0], 'Color', [0.1 0.1 0.7],'Linestyle','--','Linewidth',2)
axis([0 pmax 0 mmax])



%plot()