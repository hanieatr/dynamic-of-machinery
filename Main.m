clear
clc
%% Parameters (mm, rad)
h1=0;
h2=50;
l2=50;
l5=250;
alfa=30*pi/180;
AB=250;
AC=100;
W_2=20;
acc_2=0;
theta_2_0=-16.64014122*pi/180;
Ending_time=1;
step=0.001;
Ro=150;
theta_Ro=pi/12;
%%
X_0=[0,2,297.5,276.9,-3.92,1.86,275.47,490,-17.76,-82.63,-23047,8239];
II=Ending_time/step;
for i=0:II
   time(i+1)=i*step;
   theta_2(i+1)=0.5*acc_2*time(i+1)^2+W_2*time(i+1)+theta_2_0;
   theta_2_dot(i+1)=acc_2*time(i+1)+W_2;
   theta_2_wdot(i+1)=acc_2;
   ans1=fsolve(@DoubleCrank,X_0,[],theta_2(i+1),theta_2_dot(i+1),theta_2_wdot(i+1),alfa,l2,l5,AB,AC,h1,h2)
   answer(i+1,:)=ans1;
   X_0=ans1;
% vel and acc
% Link 2
V_G_2(:,i+1)=l2/2*theta_2_dot(i+1)*[-sin(theta_2(i+1))
                                    cos(theta_2(i+1))];
A_G_2(:,i+1)=[-(l2/2*theta_2_wdot(i+1))*sin(theta_2(i+1))-(l2/2*theta_2_dot(i+1)^2)*cos(theta_2(i+1))
               (l2/2*theta_2_wdot(i+1))*cos(theta_2(i+1))-(l2/2*theta_2_dot(i+1)^2)*sin(theta_2(i+1))];
V_A(:,i+1)=l2*theta_2_dot(i+1)*[-sin(theta_2(i+1))
                                cos(theta_2(i+1))];
A_A(:,i+1)=[-(l2*theta_2_wdot(i+1))*sin(theta_2(i+1))-(l2*theta_2_dot(i+1)^2)*cos(theta_2(i+1))
             (l2*theta_2_wdot(i+1))*cos(theta_2(i+1))-(l2*theta_2_dot(i+1)^2)*sin(theta_2(i+1))];
% Link3
V_G_3(:,i+1)=V_A(:,i+1)+Ro*ans1(5)*[-sin(ans1(1)+theta_Ro);cos(ans1(1)+theta_Ro)];
A_G_3(:,i+1)=A_A(:,i+1)+[-Ro*ans1(9)*sin(ans1(1)+theta_Ro)-Ro*ans1(5)^2*cos(ans1(1)+theta_Ro);Ro*ans1(9)*cos(ans1(1)+theta_Ro)-Ro*ans1(5)^2*sin(ans1(1)+theta_Ro)];
V_B(:,i+1)=V_A(:,i+1)+AB*ans1(5)*[-sin(ans1(1));cos(ans1(1))];
A_B(:,i+1)=A_A(:,i+1)+[-AB*ans1(9)*sin(ans1(1))-AB*ans1(5)^2*cos(ans1(1));AB*ans1(9)*cos(ans1(1))-AB*ans1(5)^2*sin(ans1(1))];
V_C(:,i+1)=V_A(:,i+1)+AC*ans1(5)*[-sin(ans1(1)+alfa);cos(ans1(1)+alfa)];
A_C(:,i+1)=A_A(:,i+1)+[-AC*ans1(9)*sin(ans1(1)+alfa)-AC*ans1(5)^2*cos(ans1(1)+alfa);AC*ans1(9)*cos(ans1(1)+alfa)-AC*ans1(5)^2*sin(ans1(1)+alfa)];
% Link5
V_G_5(:,i+1)=V_C(:,i+1)+l5/2*ans1(6)*[-sin(ans1(2));cos(ans1(2))];
A_G_5(:,i+1)=A_C(:,i+1)+[-l5/2*ans1(10)*sin(ans1(2))-l5/2*ans1(6)^2*cos(ans1(2));l5/2*ans1(10)*cos(ans1(2))-l5/2*ans1(6)^2*sin(ans1(2))];
end
% % Plots
% %link6
% plot(time,answer(:,8));
% xlabel('t(s)')
% ylabel('V_6(mm/s)')
% title('V_6 - time')
% figure
% plot(time,answer(:,12));
% xlabel('t(s)')
% ylabel('A_6(mm/s^2)')
% title('A_6 - time')
% %link4
% figure
% plot(time,answer(:,7));
% xlabel('t(s)')
% ylabel('V_4(mm/s)')
% title('V_4 - time')
% figure
% plot(time,answer(:,11));
% xlabel('t(s)')
% ylabel('A_4(mm/s^2)')
% title('A_4 - time')
% %link5
% figure
% plot(time,(V_G_5(1,:).^2+V_G_5(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('V_5(mm/s)')
% title('absolute value of V_5 - time')
% figure
% plot(time,(A_G_5(1,:).^2+A_G_5(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('A_5(mm/s^2)')
% title('absolute value of A_5 - time')
% %link3
% figure
% plot(time,(V_G_3(1,:).^2+V_G_3(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('V_3(mm/s)')
% title('absolute value of V_3 - time')
% figure
% plot(time,(A_G_3(1,:).^2+A_G_3(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('A_3(mm/s^2)')
% title('absolute value of A_3 - time')
% %Link2
% figure
% plot(time,(V_G_2(1,:).^2+V_G_2(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('V_2(mm/s)')
% title('absolute value of V_2 - time')
% figure
% plot(time,(A_G_2(1,:).^2+A_G_2(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('A_2(mm/s^2)')
% title('absolute value of A_2 - time')
% %joint A
% figure
% plot(time,(V_A(1,:).^2+V_A(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('V_A(mm/s)')
% title('absolute value of V_A - time')
% figure
% plot(time,(A_A(1,:).^2+A_A(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('A_A(mm/s^2)')
% title('absolute value of A_A - time')
% %joint B
% figure
% plot(time,(V_B(1,:).^2+V_B(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('V_B(mm/s)')
% title('absolute value of V_B - time')
% figure
% plot(time,(A_B(1,:).^2+A_B(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('A_B(mm/s^2)')
% title('absolute value of A_B - time')
% %joint C
% figure
% plot(time,(V_C(1,:).^2+V_C(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('V_C(mm/s)')
% title('absolute value of V_C - time')
% figure
% plot(time,(A_C(1,:).^2+A_C(2,:).^2).^0.5);
% xlabel('t(s)')
% ylabel('A_C(mm/s^2)')
% title('absolute value of A_C - time')
% %theta_3
% plot(time,answer(:,5));
% xlabel('t(s)')
% ylabel('\theta_3 dot(rad/s)')
% title('\theta_3 dot - time')
% figure
% plot(time,answer(:,9));
% xlabel('t(s)')
% ylabel('\theta_3 wdot(rad/s^2)')
% title('\theta_3 wdot - time')
% %theta5
% figure
% plot(time,answer(:,6));
% xlabel('t(s)')
% ylabel('\theta_5 dot(rad/s)')
% title('\theta_5 dot - time')
% figure
% plot(time,answer(:,10));
% xlabel('t(s)')
% ylabel('\theta_5 wdot(rad/s^2)')
% title('\theta_5 wdot - time')