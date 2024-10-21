function F = DoubleCrank(X,theta_2,theta_2_dot,theta_2_wdot,alfa,l2,l5,AB,AC,h1,h2)
%
F(1)=l2*cos(theta_2)+AB*cos(X(1))-X(3);
F(2)=l2*sin(theta_2)+AB*sin(X(1))+h1;
F(3)=l2*cos(theta_2)+AC*cos(X(1)+alfa)+l5*cos(X(2))-h2;
F(4)=l2*sin(theta_2)+AC*sin(X(1)+alfa)+l5*sin(X(2))-X(4);
%
F(5)=l2*theta_2_dot*cos(theta_2)+AB*X(5)*cos(X(1));
F(6)=-l2*theta_2_dot*sin(theta_2)-AB*X(5)*sin(X(1))-X(7);
F(7)=l2*theta_2_dot*cos(theta_2)+AC*X(5)*cos(X(1)+alfa)+l5*X(6)*cos(X(2))-X(8);
F(8)=-l2*theta_2_dot*sin(theta_2)-AC*X(5)*sin(X(1)+alfa)-l5*X(6)*sin(X(2));
%
F(9)=l2*(-theta_2_wdot*sin(theta_2)-theta_2_dot^2*cos(theta_2))+AB*(-X(9)*sin(X(1))-X(5)^2*cos(X(1)))-X(11);
F(10)=l2*(theta_2_wdot*cos(theta_2)-theta_2_dot^2*sin(theta_2))+AB*(X(9)*cos(X(1))-X(5)^2*sin(X(1)));
F(11)=l2*(-theta_2_wdot*sin(theta_2)-theta_2_dot^2*cos(theta_2))+AC*(-X(9)*sin(X(1)+alfa)-X(5)^2*cos(X(1)+alfa))+l5*(-X(10)*sin(X(2))-X(6)^2*cos(X(2)));
F(12)=l2*(theta_2_wdot*cos(theta_2)-theta_2_dot^2*sin(theta_2))+AC*(X(9)*cos(X(1)+alfa)-X(5)^2*sin(X(1)+alfa))+l5*(X(10)*cos(X(2))-X(6)^2*sin(X(2)))-X(12);
end

