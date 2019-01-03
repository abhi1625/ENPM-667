function s_dot=nonlinear(s,t,params,F) 
x = s(1);
x_dot = s(2);
t1 = s(3);
t1_dot = s(4);
t2 = s(5);
t2_dot = s(6);
s_dot=zeros(6,1);
s_dot(1) = x_dot;
s_dot(2) = (F-params.m1*(params.g*sin(t1)*cos(t1)+params.l1*sin(t1)*t1_dot^2)-params.m2*(params.g*sin(t2)*cos(t2)+params.l2*sin(t2)*t2_dot^2))/(params.M+params.m1*sin(t1)^2+params.m2*sin(t2)^2);
s_dot(3) = t1_dot;
s_dot(4) = (cos(t1)/params.l1)*((F-params.m1*(params.g*sin(t1)*cos(t1)+params.l1*sin(t1)*t1_dot^2)-params.m2*(params.g*sin(t2)*cos(t2)+params.l2*sin(t2)*t2_dot^2))/(params.M+params.m1*sin(t1)^2+params.m2*sin(t2)^2))-(params.g*sin(t1)/params.l1);
s_dot(5) = t2_dot;
s_dot(6) = (cos(t2)/params.l2)*((F-params.m1*(params.g*sin(t1)*cos(t1)+params.l1*sin(t1)*t1_dot^2)-params.m2*(params.g*sin(t2)*cos(t2)+params.l2*sin(t2)*t2_dot^2))/(params.M+params.m1*sin(t1)^2+params.m2*sin(t2)^2))-(params.g*sin(t2)/params.l2);
end