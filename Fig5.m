% LINKED INVASION PROJECT
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------
clear all
close all

global rp_n rp_i  mup_n mup_i mum_n mum_i d ...
       qhp_n qhp_i qcp_n qcp_i qhm_n qhm_i qcm_n qcm_i ...
       alpha_nn alpha_ni alpha_in alpha_ii beta_nn beta_ni beta_in beta_ii ...
       cp_ni cp_in cm_ni cm_in


% Standard parameters
rp_n = 0.02;
rp_i = 0.02;
mup_n = 0.1;
mup_i = 0.1;
mum_n = 0.1;
mum_i = 0.1;
d = 2;

% Conversion efficiencies
qhp_n = 5;
qhp_i = 5;
qcp_n = 1;
qcp_i = 1;
qcm_n = 1;
qcm_i = 1;
qhm_n = 1;
qhm_i = 1;


% Competition parameters

cp_ni = 0.0;
cp_in = 0.0;
cm_ni = 0.0;
cm_in = 0.0;

% Resource exchange parameters

% native to native
alpha_nn = 0.4;
beta_nn = 0.4;
% invasive to invasive
alpha_ii = 0.4;
beta_ii = 0.4;
% native to invasive
alpha_ni = 0;
beta_in = 0;
% invasive to invasive
alpha_in = 0;
beta_ni = 0;







% initial conditions of pn and mn (based on steady state value, in the
% absence of pi and mi

% initial conditions

dt = 0.1;
Tfin = 100;

pn0 = 0.1;
mn0 = 0.1;
pi0 = 0;
mi0 = 0;

pn = zeros(1,Tfin/dt);
mn = zeros(1,Tfin/dt);
pi = zeros(1,Tfin/dt);
mi = zeros(1,Tfin/dt);

pn(1) = pn0;
mn(1) = mn0;
pi(1) = pi0;
mi(1) = mi0;


Qpnn = qhp_n*alpha_nn/d - qcp_n*beta_nn;
Qpii = qhp_i*alpha_ii/d - qcp_i*beta_ii; 
Qpin = qhp_n*alpha_in/d - qcp_n*beta_ni; 
Qpni = qhp_i*alpha_ni/d - qcp_i*beta_in; 

Qmnn = qcm_n*beta_nn-qhm_n*alpha_nn/d;
Qmii = qcm_i*beta_ii-qhm_i*alpha_ii/d;
Qmin = qcm_n*beta_in-qhm_n*alpha_ni/d;
Qmni = qcm_i*beta_ni-qhm_i*alpha_in/d;


for i = 1:(Tfin/dt-1)



pn(i+1) = pn(i) + (rp_n*pn(i)+(mn(i)*Qpnn+mi(i)*Qpin)*pn(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cp_in*pn(i)*pi(i)-mup_n*pn(i)^2)*dt;
mn(i+1) = mn(i) + ((pn(i)*Qmnn + pi(i)*Qmin)*mn(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cm_in*mn(i)*mi(i)-mum_n*mn(i)^2)*dt ;
pi(i+1) = pi(i) + (rp_i*pi(i)+(mi(i)*Qpii+mn(i)*Qpni)*pi(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cp_ni*pn(i)*pi(i)-mup_i*pi(i)^2)*dt;
mi(i+1) = mi(i) + ((pi(i)*Qmii + pn(i)*Qmni)*mi(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cm_ni*mn(i)*mi(i)-mum_i*mi(i)^2)*dt ;


end


T = dt:dt:Tfin;

pn0 = pn(end);
mn0 = mn(end);

%[T,Y] = ode45(@eq_2p2f, 0:.1:Tfin, [0; 0; pi0; mi0], options);

pi0 = 0.1;%Y(end,3);
mi0 = 0.1;%Y(end,4);



% Examine choice:
% 1: RISK: Symbiont spillover and co-invasion 
% 2: RISK: Symbiont spillback and microbial invasion

choice = 2;


if choice == 1
dt = 0.1;
Tfin = 150;



% 1. RISK: Symbiont spillback and co-invasion
% invasive to native
alpha_in = 0;
beta_ni = 0;
% native to invasive
alpha_ni = 0.4; %0.4
beta_in = 0.3;  %0.3
% invasive to invasive
alpha_ii = 0.3;
beta_ii = 0.3;
% invasive microbes are stronger competitors
cm_in = 0.2;






pn = zeros(1,Tfin/dt);
mn = zeros(1,Tfin/dt);
pi = zeros(1,Tfin/dt);
mi = zeros(1,Tfin/dt);

pn(1) = pn0;
mn(1) = mn0;
pi(1) = pi0;
mi(1) = mi0;


Qpnn = qhp_n*alpha_nn/d - qcp_n*beta_nn;
Qpii = qhp_i*alpha_ii/d - qcp_i*beta_ii; 
Qpin = qhp_n*alpha_in/d - qcp_n*beta_ni; 
Qpni = qhp_i*alpha_ni/d - qcp_i*beta_in; 

Qmnn = qcm_n*beta_nn-qhm_n*alpha_nn/d;
Qmii = qcm_i*beta_ii-qhm_i*alpha_ii/d;
Qmin = qcm_n*beta_in-qhm_n*alpha_ni/d;
Qmni = qcm_i*beta_ni-qhm_i*alpha_in/d;


for i = 1:(Tfin/dt-1)



pn(i+1) = pn(i) + (rp_n*pn(i)+(mn(i)*Qpnn+mi(i)*Qpin)*pn(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cp_in*pn(i)*pi(i)-mup_n*pn(i)^2)*dt;
mn(i+1) = mn(i) + ((pn(i)*Qmnn + pi(i)*Qmin)*mn(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cm_in*mn(i)*mi(i)-mum_n*mn(i)^2)*dt ;
pi(i+1) = pi(i) + (rp_i*pi(i)+(mi(i)*Qpii+mn(i)*Qpni)*pi(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cp_ni*pn(i)*pi(i)-mup_i*pi(i)^2)*dt;
mi(i+1) = mi(i) + ((pi(i)*Qmii + pn(i)*Qmni)*mi(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cm_ni*mn(i)*mi(i)-mum_i*mi(i)^2)*dt ;


end


T = dt:dt:Tfin;




figure(1)
%subplot(2,1,1)
h1 = plot(T,pn,'Color',[0.294  0.651 0.184],'linewidth',1.5);
hold on
h2 = plot(T,mn,'Color',[0.294  0.651 0.184],'Linestyle','--','linewidth',1.5); 
hold on
h3 = plot(T,pi,'k-','linewidth',1.5) ;
hold on 
h4 = plot(T,mi,'k--','linewidth',1.5) ;
xlabel('Time [arbitrary unit]')
ylabel('Biomass [arbitrary unit]')
%title('1 plant, X fungi')
set(gca,'fontsize',12)
hold on
line([0, Tfin] ,[ pn0 ,pn0],'Color',[0.294  0.651 0.184],'Linestyle',':')
hold on
line([0, Tfin] ,[ mn0 ,mn0],'Color',[0.294  0.651 0.184],'Linestyle',':')
legend([h1,h2,h3,h4],{'p_n', 'm_n','p_i','m_i'}, 'Location','northeast', 'FontSize',12,'Orientation','horizontal')
legend boxoff



end

if choice == 2

Tfin = 80;


% 2. RISK: Symbiont spillover and microbial invasion
% invasive to native
alpha_in = 0.3;
beta_ni = 0.4;
% native to invasive
alpha_ni = 0;
beta_in = 0;
% invasive to invasive
alpha_ii = 0.3;
beta_ii = 0.3;

cp_ni = 0.2;
cp_in = 0.2;
cm_in = 0.2;





dt = 0.1;

pn = zeros(1,Tfin/dt);
mn = zeros(1,Tfin/dt);
pi = zeros(1,Tfin/dt);
mi = zeros(1,Tfin/dt);

pn(1) = pn0;
mn(1) = mn0;
pi(1) = pi0;
mi(1) = mi0;


Qpnn = qhp_n*alpha_nn/d - qcp_n*beta_nn;
Qpii = qhp_i*alpha_ii/d - qcp_i*beta_ii; 
Qpin = qhp_n*alpha_in/d - qcp_n*beta_ni; 
Qpni = qhp_i*alpha_ni/d - qcp_i*beta_in; 

Qmnn = qcm_n*beta_nn-qhm_n*alpha_nn/d;
Qmii = qcm_i*beta_ii-qhm_i*alpha_ii/d;
Qmin = qcm_n*beta_in-qhm_n*alpha_ni/d;
Qmni = qcm_i*beta_ni-qhm_i*alpha_in/d;


for i = 1:(Tfin/dt-1)



pn(i+1) = pn(i) + (rp_n*pn(i)+(mn(i)*Qpnn+mi(i)*Qpin)*pn(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cp_in*pn(i)*pi(i)-mup_n*pn(i)^2)*dt;
mn(i+1) = mn(i) + ((pn(i)*Qmnn + pi(i)*Qmin)*mn(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cm_in*mn(i)*mi(i)-mum_n*mn(i)^2)*dt ;
pi(i+1) = pi(i) + (rp_i*pi(i)+(mi(i)*Qpii+mn(i)*Qpni)*pi(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cp_ni*pn(i)*pi(i)-mup_i*pi(i)^2)*dt;
mi(i+1) = mi(i) + ((pi(i)*Qmii + pn(i)*Qmni)*mi(i)/(mn(i)+mi(i)+pn(i)/d+pi(i)/d) - cm_ni*mn(i)*mi(i)-mum_i*mi(i)^2)*dt ;


end


T = dt:dt:Tfin;




figure(1)
%subplot(2,1,2)
plot(T,pn,'Color',[0.294  0.651 0.184],'linewidth',1.5)
hold on
plot(T,mn,'Color',[0.294  0.651 0.184],'Linestyle','--','linewidth',1.5) 
hold on
plot(T,pi,'k-','linewidth',1.5)
hold on
plot(T,mi,'k--','linewidth',1.5)
%legend({'p_n', 'm_n','p_i','m_i'}, 'Location','Southeast', 'FontSize',14)
xlabel('Time [arbitrary unit]')
ylabel('Biomass [arbitrary unit]')
%title('1 plant, X fungi')
set(gca,'fontsize',12)
hold on

line([0, Tfin] ,[ pn0 ,pn0],'Color',[0.294  0.651 0.184],'Linestyle',':')
hold on
line([0, Tfin] ,[ mn0 ,mn0],'Color',[0.294  0.651 0.184],'Linestyle',':')


end
