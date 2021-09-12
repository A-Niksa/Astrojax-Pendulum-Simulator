clear all
close all
clc

%% Defining Symbolic Parameters
syms t q1 q2 q3 q4 q5 % q1:=theta1 q2:=theta2 q3:=phi1 q4:=phi2 q5:=delta
syms dq1 dq2 dq3 dq4 dq5 % dq1=dtheta1/dt dq2=dtheta2/dt dq3=dphi1/dt
                         % dq4=dphi2/dt   dq5=ddelta/dt
syms m1 m2 g R l0 x0 y0 omega psi_0 r % mass1, mass2, gravitational acceleration, pulley radius

%% Defining Balls' Positions in Cartesian Coordinate System
x1 = q5*cos(q3)*sin(q1)  +R; 
y1 = -q5*cos(q3)*cos(q1);
z1 = -q5*sin(q3)*sin(q1);       
x2 = x1 + ((l0 + sqrt((x0-r*cos(psi_0))^2+(y0-r*sin(psi_0))^2) - sqrt((x0-r*cos(omega*t+psi_0))^2+(y0-r*sin(omega*t+psi_0))^2))-q5)*cos(q4)*sin(q2);
y2 = y1 - ((l0 + sqrt((x0-r*cos(psi_0))^2+(y0-r*sin(psi_0))^2) - sqrt((x0-r*cos(omega*t+psi_0))^2+(y0-r*sin(omega*t+psi_0))^2))-q5)*cos(q4)*cos(q2);
z2 = z1 - ((l0 + sqrt((x0-r*cos(psi_0))^2+(y0-r*sin(psi_0))^2) - sqrt((x0-r*cos(omega*t+psi_0))^2+(y0-r*sin(omega*t+psi_0))^2))-q5)*sin(q4)*sin(q2);

%% Deducing Velocities of the Balls
vx1 = diff(x1,q1)*dq1 + diff(x1,q3)*dq3 + diff(x1,q5)*dq5;
vy1 = diff(y1,q1)*dq1 + diff(y1,q3)*dq3 + diff(y1,q5)*dq5;
vz1 = diff(z1,q1)*dq1 + diff(z1,q3)*dq3 + diff(z1,q5)*dq5;
vx2 = vx1 + diff(x2,q2)*dq2 + diff(x2,q4)*dq4 + diff(x2,q5)*dq5;
vy2 = vy1 + diff(y2,q2)*dq2 + diff(y2,q4)*dq4 + diff(y2,q5)*dq5;
vz2 = vz1 + diff(z2,q2)*dq2 + diff(z2,q4)*dq4 + diff(z2,q5)*dq5;

V1 = [vx1;vy1;vz1]; % velocity matrix of ball 1
V2 = [vx2;vy2;vz2]; % velocity matrix of ball 2

%% Defining the Lagrangian
T = m1/2*(V1.')*V1 + m2/2*(V2.')*V2;
V = m1*g*y1 + m2*g*y2;
L = T-V;

%% Generalized Coordinates
q = [q1;q2;q3;q4;q5];
dq = [dq1;dq2;dq3;dq4;dq5];

%% Mass Matrix Deduction
dL_ddq = jacobian(L,dq).';
M = simplify(jacobian(dL_ddq,dq)); % mass matrix
N = simplify(jacobian(dL_ddq,q));
B = N*dq + diff(dL_ddq,t) - jacobian(L,q).';