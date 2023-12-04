% A UNIFYING THEORETICAL FRAMEWORK FOR MICROBIAL-MEDIATED INVASION
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




choice = 9;
% choice can be as shown below (compare with Fig. 1 and Table 1).
% 1: (-,+),
% 2.1: (0,+) - symbiont competition
% 2.2: (0,+) - plant competition
% 3: (+,+),
% 4: (-,0)
% 5: (0,0),
% 6.1: (+,0) - symbiont competition
% 6.2: (+,0) - plant competition
% 7: (-,-),
% 8: (0,-)
% 9: (+,-)




% Standard parameters
rp_n = 0;
rp_i = 0;
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

cp_ni = 0.4;
cp_in = 0.4;
cm_ni = 0.4;
cm_in = 0.4;

% Resource exchange parameters

% native to native
alpha_nn = 0.4;
alpha_ii = 0.4;
% invasive to invasive
beta_nn = 0.4;
beta_ii = 0.4;


% initial conditions

pn0 = .1;
mn0 = .1;
pi0 = .1;
mi0 = .1;

% Length of simulation

Tfin = 200;





if choice == 1 % scenario (-,+)
% invasive to native
alpha_in = 0;
beta_ni = 0.4;
% native to invasive
alpha_ni = 0.4;
beta_in = 0.4;
end


if choice == 2.1 % scenario (0,+) - symbiont competition
% invasive to native
alpha_in = 0;
beta_ni = 0;
% native to invasive
alpha_ni = 0.4;
beta_in = 0.4;
% competition
cp_in = 0;
cp_ni = 0;
end


if choice == 2.2 % scenario (0,+) - plant competition
% invasive to native
alpha_in = 0;
beta_ni = 0;
% native to invasive
alpha_ni = 0.4;
beta_in = 0.4;
% competition
cm_in = 0;
cm_ni = 0;
end


if choice == 3   % scenario (+,+) 
% invasive fungus to native plant
alpha_in = 0.4;
beta_ni = 0.4;
% native fungus to invasive plant
alpha_ni = 0.4;
beta_in = 0.4;
end




if choice == 4 % scenario (-,0)
% invasive fungus to native plant
alpha_in = 0;
beta_ni = 0.4;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0;
end




if choice == 5 % scenario (0,0)
% invasive fungus to native plant
alpha_in = 0;
beta_ni = 0;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0;
end



if choice == 6.1  % scenario (+,0) - symbiont competition
% invasive fungus to native plant
alpha_in = 0.4;
beta_ni = 0.4;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0;
%competition
cp_ni = 0;
cp_in = 0;
end



if choice == 6.2  % scenario (+,0) - plant competition
% invasive fungus to native plant
alpha_in = 0.4;
beta_ni = 0.4;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0;
%competition
cm_ni = 0;
cm_in = 0;
end



if choice == 7 % scenario (-,-)
% invasive fungus to native plant
alpha_in = 0;
beta_ni = 0.4;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0.4;
end



if choice == 8  % scenario (0,-)
% invasive fungus to native plant
alpha_in = 0;
beta_ni = 0;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0.4;
end



if choice == 9  % scenario (+,-)
% invasive fungus to native plant
alpha_in = 0.4;
beta_ni = 0.4;
% native fungus to invasive plant
alpha_ni = 0;
beta_in = 0.4;
end


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
plot(T,pn,'g','linewidth',1.5)
hold on
plot(T,mn,'g--','linewidth',1.5) 
hold on
plot(T,pi,'k-','linewidth',1.5)
hold on
plot(T,mi,'k--','linewidth',1.5)
legend({'p_n', 'm_n','p_i','m_i'}, 'Location','Northeast', 'FontSize',14)
xlabel('Time [arbitrary unit]')
ylabel('Biomass [arbitrary unit]')
%title('1 plant, X fungi')
set(gca,'fontsize',14)





