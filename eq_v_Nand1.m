% LINKED INVASION PROJECT
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = eq_v_Nand1(t,y)

global qhp qcm qcp qhm mup mum d N alpha beta cm

p = y(1);
m = y(2:(N+1));
dy = zeros(N+1,1);

Q = transpose(qhp*alpha/d - qcp*beta); 
Q1 = transpose(qcm*beta-qhm*alpha/d);



dy(1) = (sum(m.*Q))*p/(sum(m)+p/d) - mup*p^2;
dy(2) = p*m(1)*Q1(1)/(sum(m)+p/d) - cm*m(1)*sum(m(2:N))- mum*m(1)^2 ;
for i = 2:(N)
dy(i+1) = p*m(i)*Q1(i)/(sum(m)+p/d) - cm*m(i)*m(1)- mum*m(i)^2 ;
end

end