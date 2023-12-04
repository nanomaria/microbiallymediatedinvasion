close all
clear all

global qhp qcm qcp qhm mup mum d alpha beta rp



% model parameters
mup = 0.1;
mum = 0.1;
alpha = 0.4;
beta = 0.4;
rp = 0.1;
qcm = 1;
qhm = 1;
qcp = 1;
qhp = 5;
d = 2;

pmax = 7;
mmax = 7;


Q = qhp*alpha/d - qcp*beta;
Q1 = qcm*beta-qhm*alpha/d;

f = @basic_eq_rp;
t1 = 7;
t2 = 7;
y1 = linspace(0,t1,16); % from zero to t1, 30 arrows
y2 = linspace(0,t2,16);
[x,y] = meshgrid(y1, y2);
size(x);
size(y);

u = zeros(size(x));
v = zeros(size(y));

t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

syms un vn

un = u./sqrt(u.^2 + v.^2);
vn = v./sqrt(v.^2 + v.^2);


figure(1);
subplot(1,2,1)
quiver(x, y, un, vn, 0.7,'Color',[0.7 0.7 0.7],'linewidth',1 )
syms mv pv;
hold on
n1 = fimplicit(@(pv,mv)  rp.*pv+pv.*mv.*Q./(pv/d+mv) - mup.*pv^2, [0,pmax],'Color', [0.6 0.3 0],'Linestyle','--','Linewidth',2);
hold on
n2 = fimplicit(@(pv,mv) mv.*pv*Q1./(pv/d+mv) - mum.*mv^2, [0,pmax],'Color', [0.1 0.1 0.7],'Linestyle','--','Linewidth',2);
hold on
axis([ 0 pmax 0 mmax]);
hold on
y10 = 0.1;
x10 = 0.2;
[ts,ys] = ode45(@basic_eq_rp,0:0.1:1000,[x10;y10]);
plot(ys(:,1),ys(:,2),'Color',[1 0.8 0],'linewidth',1.5)
plot(x10,y10,'Color',[1 0.8 0],'MarkerSize',20)
plot(ys(end,1),ys(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
hold on
y10_rp = 0;
x10_rp = 0.1;
[ts,ys_rp] = ode45(@basic_eq_rp,0:0.1:1000,[x10_rp;y10_rp]);
plot(ys_rp(end,1),ys_rp(end,2),'s','MarkerSize',15,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
plot([0,0],[0,mmax], 'Color', [0.6 0.3 0],'Linestyle','--','Linewidth',2)
hold on
line([0,pmax],[0,0], 'Color', [0.1 0.1 0.7],'Linestyle','--','Linewidth',2)
hold on
line([Q/mup, Q/mup],[0,mmax],'Color',[0 0 0], 'Linestyle',':','Linewidth',1.2)
set(gca,'fontsize',14)
xlabel('p')
ylabel('m')
%Squared


t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end




[T,Y] = ode45(@basic_eq_rp,0:0.1:1000,[x10;y10]);

figure(1)
subplot(1,2,2)
plot(T,Y(:,1),'Color', [0.6 0.3 0],'Linestyle','-','Linewidth',1.5)
hold on
plot(T,Y(:,2),'Color', [0.1 0.1 0.7],'Linestyle','-','Linewidth',1.5)
axis([0 100 0 max(Y(:,1))*1.1])
set(gca,'fontsize',14)
legend('p','m')
xlabel('Time')
ylabel('p,m')

