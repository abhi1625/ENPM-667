syms M m1 m2 g l1 l2 s t
params.M  = 1000;       % M      mass of cart 
params.m1 = 100;        % m1     mass of bob 1 
params.m2 = 100;        % m2     mass of bob 2 
params.l1 = 20;         % l1     length of link of first pendulum 
params.l2 = 10;         % l2     length of link of second pendulum 
params.g  = 10;         % g      acceleration due to gravitaion 
tspan = 0:0.1:100;      % t      time for simulation of each controller

% Linearization of the non-linear model
[s_dot,A,B] = lin_model(params);

%Stability
E1 = eig(A)
disp('The eigen values of A lie on the imaginary axis,')
disp('therefore Lyapunovs indirect method is inconclusive for the system')

% Controllability Condition
r = rank([B,A*B,(A^2)*B,(A^3)*B,(A^4)*B,(A^5)*B]);
if r == 6
   disp('The System is controllable, as rank(ctrb(A,B)) == 6')
end

%% LQR Controller

% Weight Matrices
R=0.0001;
Q = diag([100,1000,100000,100000,100000,100000]);

% Defining initial conditions

%s0 = [0; 0; 10*pi/180; 0; 10*pi/180; 0];  %For small change in theta
s0 = [1; 0; 0; 0; 0; 0];                   %For change in X

% Using lqr fn to solve the Algebraic Ricatti's eqn
[K,P,e] = lqr(double(A),double(B),Q,R);

% Checking eigen values of the closed loop system 
E2 = eig(A-B*K) 
disp('The eigenvalues  of the closed loop system lie in the left half plane.')
disp('Thus the closed loop system is stable.')

% Simulatiion of the LQR Controller


[t,y1] = ode45(@(t,y)cart_EOM(y,t,params,A,B,K),tspan,s0);

figure;
plot(t,y1(:,1),'k')
ylabel('x position of the cart')
xlabel('time in s')

%% Observer based Control

C1 = [1 0 0 0 0 0];             %For output (x(t))
C2 = [0 0 1 0 0 0;...           %For output (t1(t),t2(t)) 
      0 0 1 0 1 0];
C3 = [1 0 0 0 0 0;...           %For output (x(t),t2(t))    
      0 0 0 0 1 0];
C4 = [1 0 0 0 0 0;...           %For output (x(t),t1(t),t2(t))
      0 0 1 0 0 0;...
      0 0 0 0 0 1];
  
%Checking Observability  
ob1 = [C1;C1*A;C1*A^(2);C1*A^(3);C1*A^(4);C1*A^(5)];
ob2 = [C2;C2*A;C2*A^(2);C2*A^(3);C2*A^(4);C2*A^(5)];
ob3 = [C3;C3*A;C3*A^(2);C3*A^(3);C3*A^(4);C3*A^(5)];
ob4 = [C4;C4*A;C4*A^(2);C4*A^(3);C4*A^(4);C4*A^(5)];

observability1 = rank(ob1)
disp('The system is observable for output x(t)')
observability2 = rank(ob2)
disp('Not observable for output (t1(t),t2(t))')
observability3 = rank(ob3)
disp('The system is observable for output (x(t),t2(t))')
observability4 = rank(ob4)
disp('The system is observable for output (x(t),t1(t),t2(t))')

%Introducing Noise and Disturbances in the system
Bd = 0.001.*eye(6);             %input disturbance covarianve
Bn1 = 0.01;                      %output measurement noise
Bn3 = 0.01*eye(2);
Bn4 = 0.01*eye(3);

% Obtaining Luenberger Observers for output vectors of each observable system
[L1,P,E] = lqe(A,Bd,C1,Bd,Bn1);
[L3,P,E] = lqe(A,Bd,C3,Bd,Bn3);
[L4,P,E] = lqe(A,Bd,C4,Bd,Bn4);

% Creating Augmented Matrices for Simulation
uD = randn(6,size(tspan,2));      %input for disturbance
uN = randn(size(tspan));          %input for noise

u = 0*tspan;
u(100:length(tspan)) = 1;         % Step input at t = 10

u1 = [u; Bd*Bd*uD; uN];

Be = [B,Bd,zeros(size(B))];       %Augmented B matrix


%% State Estimator for C1 = [1,0,0,0,0,0]

sysLO1 = ss(A-L1*C1,[B L1],C1,zeros(1,2));     %State Estimator system

% Obtaining Y values for a system simulated with noise and disturbance. 
De1 = [0,0,0,0,0,0,0,Bn1];                     %Augmented D matrix

sys1 = ss(A,Be,C1,De1);
[y1,t] = lsim(sys1,u1,tspan);

%Simulating the States of the output variables 
[x1,t] = lsim(sysLO1,[u; y1'],tspan);

figure();
hold on
plot(t,y1(:,1),'k','Linewidth',2)
plot(t,x1(:,1),'r--','Linewidth',2)
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C1: output vector x(t)')
hold off

%% State Estimator for C3: For output (x(t),t2(t))

sysLO3 = ss(A-L3*C3,[B L3],C3,zeros(2,3));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De3 = [zeros(size(C3)),Bn3];                     %Augmented D matrix

sys3 = ss(A,Be,C3,De3);
[y3,t] = lsim(sys3,u1,tspan);

%Simulating the States of the output variables 
[x3,t] = lsim(sysLO3,[u; y3'],tspan);

figure();
hold on
plot(t,y3(:,1),'g','Linewidth',2)
plot(t,y3(:,2),'b','Linewidth',2)
plot(t,x3(:,1),'k--','Linewidth',1)
plot(t,x3(:,2),'r--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Estimated Response for C3: output vector (x(t),t2(t))')
hold off


%% State Estimator for C4: For output (x(t),t1(t),t2(t))

sysLO4 = ss(A-L4*C4,[B L4],C4,zeros(3,4));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De4 = [zeros(3,5),Bn4];                     %Augmented D matrix

sys4 = ss(A,Be,C4,De4);
[y4,t] = lsim(sys4,u1,tspan);

%Simulating the States of the output variables 
[x4,t] = lsim(sysLO4,[u; y4'],tspan);

figure();
hold on
plot(t,y4(:,1),'g','Linewidth',2)
plot(t,y4(:,2),'b','Linewidth',2)
plot(t,y4(:,3),'c','Linewidth',2)
plot(t,x4(:,1),'k--','Linewidth',1)
plot(t,x4(:,2),'r--','Linewidth',1)
plot(t,x4(:,3),'m--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off

%% LQG Controller for Output Vector C1 = [1,0,0,0,0,0]

Ac = A-L1*C1;
Bc = [B L1];
Cc = eye(6);
Dc = 0*[B L1];

opt = simset('solver','ode45','SrcWorkspace','Current');
sim('nonlinearlqg',tspan,opt);


%% Simulation Results

figure();
hold on
plot(tout,states(:,1),'r')
plot(tout,states(:,3),'k')
plot(tout,states(:,5),'b')
title('LQG Nonlinear System')
legend('x-position','theta1','theta2')
hold off

figure();
for k=1:1:size(states,1)
   drawsys(states(k,:),params);
end




