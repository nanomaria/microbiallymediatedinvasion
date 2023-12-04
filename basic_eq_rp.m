% LINKED INVASION PROJECT
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = basic_eq_rp(t,y)

global qhp qcm qcp qhm mup mum d alpha beta rp

p = y(1);
m = y(2);
dy = zeros(2,1);

dy(1) = rp*p+ qhp*alpha*m*p/(m*d+p) - qcp*beta*m*p/(p/d+m) - mup*p^2;
dy(2) = qcm*beta*m*p/(p/d+m) - qhm*alpha*m*p/(m*d+p) - mum*m^2 ;

end