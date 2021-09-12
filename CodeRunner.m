clear all
close all
clc
warning off

%% Assinging Numerical Values to Global Variables
global m1 m2 r_b1 r_b2 g R P omega r l0 x0 y0 psi_0
m1 = 6e-3; % kg | mass of the 1st ball
m2 = 6e-3; % kg | mass of the 2nd ball
r_b1 = 0.02; % m | radius of the 1st ball
r_b2 = 0.02; % m | radius of the 2nd ball
g = 9.81; % m/s^2 | approximate gravitational acceleration
R = 0.1; % m | radius of the pulley
P = 15; % N | external force
omega = 3; % Hz | angular velocity of the disk
r = 0.06; % m | radius of the rotating disk
x0 = 0.1614;
y0 = 0.2205;
l0 = 0.6;

%% Defining Initial Conditions
t1_0 = -0.1864325245643733464492437128467; % initial value of theta 1
t2_0 = -2.9892344895502714242518641736621; % initial value of theta 2
p1_0 = -0.15538628881836451988404101002139; % initial value of phi 1
p2_0 = 4.5359306224477402598133831828648; % initial value of phi 2
psi_0 = 2.8394156500266773931425638607477; % initial value of psi
d_0 = 0.33509187306138983735674760107519; % m | initial value of delta
t1_dot_0 = 0.12929891048518530353117941421072; % Hz | initial value of theta dot 1
t2_dot_0 = 0.22372192526950957642560332789711; % Hz | initial value of theta dot 2
p1_dot_0 = 0.51389935006164690867018098359447; % Hz | initial value of phi dot 1
p2_dot_0 = 1.1798568436174846957561345154955; % Hz | initial value of phi dot 2
psi_dot_0 = 0.4196253110965842400575314090554; % Hz | initial value of l dot
d_dot_0 = 0.41078279700730993890959182246041; % m/s | initial value of delta dot

init_conditions = [t1_0;t2_0;p1_0;p2_0;d_0;... 
    t1_dot_0;t2_dot_0;p1_dot_0;p2_dot_0;d_dot_0];

%% Solving the Corresponding Differential Equations
n = 40;
stepsize = 0.05;
t_end = 10;
tspan = [0:stepsize:t_end];
option = odeset('MaxStep',stepsize);
[t,z] = ode23s(@NumericalFunction,tspan,init_conditions,option);


%% Assigning Numerical Arrays to Coordinates
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5);
dq1 = z(:,6); dq2 = z(:,7); dq3 = z(:,8); dq4 = z(:,9); dq5 = z(:,10);

%% Converting Spherical Coordinates to Cartesian Coordinates
x1 = q5.*cos(q3).*sin(q1) + R; 
y1 = -q5.*cos(q3).*cos(q1);
z1 = -q5.*sin(q3).*sin(q1); 

x2 = x1 + ((l0 + sqrt((x0-r.*cos(psi_0)).^2+(y0-r.*sin(psi_0)).^2) - sqrt((x0-r.*cos(omega.*t+psi_0)).^2+(y0-r.*sin(omega.*t+psi_0)).^2))-q5).*cos(q4).*sin(q2);
y2 = y1 - ((l0 + sqrt((x0-r.*cos(psi_0)).^2+(y0-r.*sin(psi_0)).^2) - sqrt((x0-r.*cos(omega.*t+psi_0)).^2+(y0-r.*sin(omega.*t+psi_0)).^2))-q5).*cos(q4).*cos(q2);
z2 = z1 - ((l0 + sqrt((x0-r.*cos(psi_0)).^2+(y0-r.*sin(psi_0)).^2) - sqrt((x0-r.*cos(omega.*t+psi_0)).^2+(y0-r.*sin(omega.*t+psi_0)).^2))-q5).*sin(q4).*sin(q2);

%% Performing FFT
FFT_x1 = abs(fft(full(x1))/length(x1)); f_x1 = (240/61).*(0:(length(t)/2))/length(t);
FFT_x1 = FFT_x1(1:length(x1)/2+1);
FFT_y1 = abs(fft(full(y1))/length(y1)); f_y1 = (240/61).*(0:(length(t)/2))/length(t);
FFT_y1 = FFT_y1(1:length(y1)/2+1);
FFT_z1 = abs(fft(full(z1))/length(z1)); f_z1 = (240/61).*(0:(length(t)/2))/length(t);
FFT_z1 = FFT_z1(1:length(z1)/2+1);

FFT_x2 = abs(fft(full(x2))/length(x2)); f_x2 = (240/61).*(0:(length(t)/2))/length(t);
FFT_x2 = FFT_x2(1:length(x2)/2+1);
FFT_y2 = abs(fft(full(y2))/length(y2)); f_y2 = (240/61).*(0:(length(t)/2))/length(t);
FFT_y2 = FFT_y2(1:length(y2)/2+1);
FFT_z2 = abs(fft(full(z2))/length(z2)); f_z2 = (240/61).*(0:(length(t)/2))/length(t);
FFT_z2 = FFT_z2(1:length(z1)/2+1);

%% Deducing Velocity
vx1 = diff(x1)./diff(t);
vy1 = diff(y1)./diff(t);
vz1 = diff(z1)./diff(t);
vx2 = diff(x2)./diff(t);
vy2 = diff(y2)./diff(t);
vz2 = diff(z2)./diff(t);

%% Deducing Poincare Section
fixed_point = 0; % z = 0
poincare_x1 = zeros(1,length(t)); poincare_y1 = zeros(1,length(t)); 
poincare_x2 = zeros(1,length(t)); poincare_y2 = zeros(1,length(t));
poincare_vx1 = zeros(1,length(t)); poincare_vy1 = zeros(1,length(t));
poincare_vx2 = zeros(1,length(t)); poincare_vy2 = zeros(1,length(t));

for i = 1:length(t)-1
    if(z1(i)<0)
        if(z1(i+1)>0)
            poincare_x1(i) = interp1([z1(i) z1(i+1)],[x1(i) x1(i+1)],fixed_point,'spline');
            poincare_y1(i) = interp1([z1(i) z1(i+1)],[y1(i) y1(i+1)],fixed_point,'spline');
            poincare_vx1(i) = (x1(i+1)-x1(i))/(t(i+1)-t(i));
            poincare_vy1(i) = (x1(i+1)-y1(i))/(t(i+1)-t(i));    
        end
    end
end

for i = 1:length(t)-1
    if(z2(i)<0)
        if(z2(i+1)>0)
            poincare_x2(i) = interp1([z2(i) z2(i+1)],[x2(i) x2(i+1)],fixed_point,'spline');
            poincare_y2(i) = interp1([z2(i) z2(i+1)],[y2(i) y2(i+1)],fixed_point,'spline');
            poincare_vx2(i) = (x2(i+1)-x2(i))/(t(i+1)-t(i));
            poincare_vy2(i) = (y2(i+1)-y2(i))/(t(i+1)-t(i));
        end
    end
end

%% Coordinate with Respect to Time Diagram e.g. x-t
%x
fig1a = figure('Units','Normalized','Position',[0.25 0.145 0.55 0.45],'Visible','off'); %left bottom width height
subplot(2,1,1)
plot(t,x1,'Color',[0.2 .63 0.6],'LineWidth',1.1)
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','Fontname','Helvetica','FontSize',10)
grid on

subplot(2,1,2)
plot(t,x2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig1a.png','-dpng','-r500',fig1a);

%y
fig1b = figure('Units','Normalized','Position',[0.25 0.25 0.55 0.45],'Visible','off'); %left bottom width height
subplot(2,1,1)
plot(t,y1,'Color',[0.2 .63 0.6],'LineWidth',1.1)
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('y (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','Fontname','Helvetica','FontSize',10)
grid on

subplot(2,1,2)
plot(t,y2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('y (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig1b.png','-dpng','-r500',fig1b);

%z
fig1c = figure('Units','Normalized','Position',[0.25 0.25 0.55 0.45],'Visible','off'); %left bottom width height
subplot(2,1,1)
plot(t,z1,'Color',[0.2 .63 0.6],'LineWidth',1.1)
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('z (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','Fontname','Helvetica','FontSize',10)
grid on

subplot(2,1,2)
plot(t,z2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('z (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{2}','Fontname','Helvetica','FontSize',9)
grid on

print('Fig1c.png','-dpng','-r500',fig1c);

%% FFT Results
%x
fig2a = figure('Units','Normalized','Position',[0.25 0.25 0.55 0.45],'Visible','off'); %left bottom width height
subplot(2,1,1)
plot(f_x1,FFT_x1,'Color',[0.2 .63 0.6],'LineWidth',1.1)
xlabel('f (Hz)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('|x| (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','Fontname','Helvetica','FontSize',10)
grid on

subplot(2,1,2)
plot(f_x2,FFT_x2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('f (Hz)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('|x| (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig2a.png','-dpng','-r500',fig2a);

%y
fig2b = figure('Units','Normalized','Position',[0.25 0.25 0.55 0.45],'Visible','off'); %left bottom width height
subplot(2,1,1)
plot(f_y1,FFT_y1,'Color',[0.2 .63 0.6],'LineWidth',1.1)
xlabel('f (Hz)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('|y| (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','Fontname','Helvetica','FontSize',10)
grid on

subplot(2,1,2)
plot(f_y2,FFT_y2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('f (Hz)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('|y| (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig2b.png','-dpng','-r500',fig2b);

%z
fig2c = figure('Units','Normalized','Position',[0.25 0.25 0.55 0.45],'Visible','off'); %left bottom width height
subplot(2,1,1)
plot(f_z1,FFT_z1,'Color',[0.2 .63 0.6],'LineWidth',1.1)
xlabel('f (Hz)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('|z| (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','Fontname','Helvetica','FontSize',10)
grid on

subplot(2,1,2)
plot(f_z2,FFT_z2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('f (Hz)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('|z| (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig2c.png','-dpng','-r500',fig2c);

%% Phase Diagrams
%x
fig3a = figure('Units','Normalized','Position',[0.25 0.25 0.3 0.325],'Visible','off'); %left bottom width height
plot(x1(2:end),vx1,'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 .63 0.6],'MarkerFaceColor',[0.2 .63 0.6],'MarkerSize',1.15)
hold on
plot(x2(2:end),vx2,'LineStyle','none','Marker','square','MarkerEdgeColor',[0 0 .55],'MarkerFaceColor',[0 0 .55],'MarkerSize',1.15)
xlabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('v_{x} (^{m}/_{s})','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig3a.png','-dpng','-r500',fig3a);

%y
fig3b = figure('Units','Normalized','Position',[0.25 0.25 0.3 0.325],'Visible','off'); %left bottom width height
plot(y1(2:end),vy1,'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 .63 0.6],'MarkerFaceColor',[0.2 .63 0.6],'MarkerSize',1.15)
hold on
plot(y2(2:end),vy2,'LineStyle','none','Marker','square','MarkerEdgeColor',[0 0 .55],'MarkerFaceColor',[0 0 .55],'MarkerSize',1.15)
xlabel('y (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('v_{y} (^{m}/_{s})','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig3b.png','-dpng','-r500',fig3b);

%z
fig3c = figure('Units','Normalized','Position',[0.25 0.25 0.3 0.325],'Visible','off'); %left bottom width height
plot(z1(2:end),vz1,'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 .63 0.6],'MarkerFaceColor',[0.2 .63 0.6],'MarkerSize',1.15)
hold on
plot(z2(2:end),vz2,'LineStyle','none','Marker','square','MarkerEdgeColor',[0 0 .55],'MarkerFaceColor',[0 0 .55],'MarkerSize',1.15)
xlabel('z (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('v_{z} (^{m}/_{s})','FontName','Helvetica','FontSize',11,'FontWeight','bold')
legend('m_{1}','m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig3c.png','-dpng','-r500',fig3c);

%% Poincare Section
fig4a = figure('Units','Normalized','Position',[0.25 0.25 0.3 0.325],'Visible','off'); %left bottom width height
plot(poincare_x1,poincare_vx1,'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 .63 0.6],'MarkerFaceColor',[0.2 .63 0.6],'MarkerSize',2)
hold on
plot(poincare_x2,poincare_vx2,'LineStyle','none','Marker','square','MarkerEdgeColor',[0 0 .55],'MarkerFaceColor',[0 0 .55],'MarkerSize',2)
xlabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('v_{x} (^{m}/_{s})','FontName','Helvetica','FontSize',11,'FontWeight','bold')
title('plane @ z = 0','FontName','Helvetica','FontSize',11.5,'FontWeight','bold')
legend('m_{1}','m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig4a.png','-dpng','-r500',fig4a);

fig4b = figure('Units','Normalized','Position',[0.25 0.25 0.3 0.325],'Visible','off'); %left bottom width height
plot(poincare_y1,poincare_vy1,'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 .63 0.6],'MarkerFaceColor',[0.2 .63 0.6],'MarkerSize',2)
hold on
plot(poincare_y2,poincare_vy2,'LineStyle','none','Marker','square','MarkerEdgeColor',[0 0 .55],'MarkerFaceColor',[0 0 .55],'MarkerSize',2)
xlabel('y (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('v_{y} (^{m}/_{s})','FontName','Helvetica','FontSize',11,'FontWeight','bold')
title('plane @ z = 0','FontName','Helvetica','FontSize',11.5,'FontWeight','bold')
legend('m_{1}','m_{2}','Fontname','Helvetica','FontSize',10)
grid on

print('Fig4b.png','-dpng','-r500',fig4b);


%% 3D Trajectory
fig5a = figure('Units','Normalized','Position',[0.25 0.25 0.25 0.325],'Visible','off'); %left bottom width height
plot3(x1,z1,y1,'Color',[0.2 .63 0.6],'LineWidth',.8)
hold on
plot3(x2,z2,y2,'Color',[0 0 .55],'LineWidth',1.1)
xlabel('x (m)','FontName','Helvetica','FontSize',10,'FontWeight','bold')
ylabel('z (m)','FontName','Helvetica','FontSize',10,'FontWeight','bold')
zlabel('y (m)','FontName','Helvetica','FontSize',10,'FontWeight','bold')
legend('m_{1}','m_{2}','Fontname','Helvetica','FontSize',9)
grid on

print('Fig5a.png','-dpng','-r500',fig5a);

%% 2D Animation (upwards view + sideways view)
figure('units','normalized','outerposition',[0 0 1 1],'Visible','off')
subplot(3,6,1:2)
plot(t,x1,'Color',[0.2 .63 0.6],'LineWidth',.95)
hold on
plot(t,x2,'Color',[0 0 .55],'LineWidth',.95)
pointer_x_2D(1) = line(t(1),x1(1),'Marker','.','MarkerSize',30,'Color',[0.2 .63 0.6]);
pointer_x_2D(2) = line(t(1),x2(1),'Marker','.','MarkerSize',30,'Color',[0 0 .55]);
hold off
legend('x_{1}','x_{2}','FontSize',10,'FontName','Helvetica')
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
grid on

subplot(3,6,3:4)
plot(t,y1,'Color',[0.2 .63 0.6],'LineWidth',.95)
hold on
plot(t,y2,'Color',[0 0 .55],'LineWidth',.95)
pointer_y_2D(1) = line(t(1),y1(1),'Marker','.','MarkerSize',30,'Color',[0.2 .63 0.6]);
pointer_y_2D(2) = line(t(1),y2(1),'Marker','.','MarkerSize',30,'Color',[0 0 .55]);
hold off
legend('y_{1}','y_{2}','FontSize',10,'FontName','Helvetica')
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('y (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
grid on

subplot(3,6,5:6)
plot(t,z1,'Color',[0.2 .63 0.6],'LineWidth',.95)
hold on
plot(t,z2,'Color',[0 0 .55],'LineWidth',.95)
pointer_z_2D(1) = line(t(1),z1(1),'Marker','.','MarkerSize',30,'Color',[0.2 .63 0.6]);
pointer_z_2D(2) = line(t(1),z2(1),'Marker','.','MarkerSize',30,'Color',[0 0 .55]);
hold off
legend('z_{1}','z_{2}','FontSize',10,'FontName','Helvetica')
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('z (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
grid on

subplot(3,6,[7:9 13:15])
bottomview = plot([0 x1(1);x1(1) x2(1)],[0 z1(1);z1(1) z2(1)],...
    '.-','MarkerSize',20,'LineWidth',1.25,'Color','k','MarkerFaceColor','b','MarkerEdgeColor','b');
axis equal
grid on
xlabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('z (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
axis([-0.6 0.6 -0.6 0.6]);
bottomtitle = title(sprintf('Time: %0.2f s',full(t(1))),'FontName','Helvetica','FontSize',13);

subplot(3,6,[10:12 16:18])
frontview = plot([0 x1(1);x1(1) x2(1)],[0 y1(1);y1(1) y2(1)],...
    '.-','MarkerSize',20,'LineWidth',1.25,'Color','k','MarkerFaceColor','b','MarkerEdgeColor','b');
axis equal
grid on
xlabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('y (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
axis([-0.5 0.5 -0.5 0.1]);
fronttitle = title(sprintf('Time: %0.2f s',full(t)),'FontName','Helvetica','FontSize',13);

mov_2D(1,length(t)) = struct('cdata',[],'colormap',[]);


for i = 1:length(t)
    set(pointer_x_2D(1),'XData',t(i),'YData',x1(i));
    set(pointer_x_2D(2),'Xdata',t(i),'YData',x2(i));
    set(pointer_y_2D(1),'XData',t(i),'YData',y1(i));
    set(pointer_y_2D(2),'Xdata',t(i),'YData',y2(i));
    set(pointer_z_2D(1),'XData',t(i),'YData',z1(i));
    set(pointer_z_2D(2),'Xdata',t(i),'YData',z2(i));
    set(frontview(1),'XData',[0,x1(i)],'YData',[0,y1(i)]);
    set(frontview(2),'XData',[x1(i),x2(i)],'YData',[y1(i),y2(i)]);
    set(fronttitle,'String',sprintf('Time: %0.2f s',full(t(i))));
    set(bottomview(1),'XData',[0,x1(i)],'YData',[0,z1(i)]);
    set(bottomview(2),'XData',[x1(i),x2(i)],'YData',[z1(i),z2(i)]);
    set(bottomtitle,'String',sprintf('Time: %0.2f s',full(t(i))));
    pause(0.1);
    mov_2D(i) = getframe(gcf);
end

movie = VideoWriter('mov_2D.avi');
open(movie)
writeVideo(movie,mov_2D)
close(movie)

%% 3D Animation
%mov_3D.avi
figure('units','normalized','outerposition',[0 0 1 1],'Visible','off')
subplot(3,6,1:2)
plot(t,x1,'Color',[.2 .63 .6],'LineWidth',.95)
hold on
plot(t,x2,'Color',[0 0 .55],'LineWidth',.95)
pointer_x(1) = line(t(1),x1(1),'Marker','.','MarkerSize',30,'Color',[0.2 .63 0.6]);
pointer_x(2) = line(t(1),x2(1),'Marker','.','MarkerSize',30,'Color',[0 0 .55]);
hold off
legend('x_{1}','x_{2}','FontSize',10,'FontName','Helvetica')
xlabel('t (s)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
ylabel('x (m)','FontName','Helvetica','FontSize',11,'FontWeight','bold')
grid on

subplot(3,6,3:4)
plot(t,y1,'Color',[0.2 .63 0.6],'LineWidth',.95)
hold on
plot(t,y2,'Color',[0 0 .55],'LineWidth',.95)
pointer_y(1) = line(t(1),y1(1),'Marker','.','MarkerSize',30,'Color',[0.2 .63 0.6]);
pointer_y(2) = line(t(1),y2(1),'Marker','.','MarkerSize',30,'Color',[0 0 .55]);
hold off
legend('y_{1}','y_{2}','FontSize',10,'FontName','Helvetica')
xlabel('t (s)','FontName','Helvetica','FontSize',12,'FontWeight','bold')
ylabel('y (m)','FontName','Helvetica','FontSize',12,'FontWeight','bold')
grid on

subplot(3,6,5:6)
plot(t,z1,'Color',[0.2 .63 0.6],'LineWidth',.95)
hold on
plot(t,z2,'Color',[0 0 .55],'LineWidth',.95)
pointer_z(1) = line(t(1),z1(1),'Marker','.','MarkerSize',30,'Color',[0.2 .63 0.6]);
pointer_z(2) = line(t(1),z2(1),'Marker','.','MarkerSize',30,'Color',[0 0 .55]);
hold off
legend('z_{1}','z_{2}','FontSize',10,'FontName','Helvetica')
xlabel('t (s)','FontName','Helvetica','FontSize',12,'FontWeight','bold')
ylabel('z (m)','FontName','Helvetica','FontSize',12,'FontWeight','bold')
grid on

subplot(3,6,[9:10 15:16])
spatial = plot3([0 x1(1);x1(1) x2(1)],[0 y1(1);y1(1) y2(1)],[0 y1(1);z1(1) z2(1)],...
    '.-','MarkerSize',20,'LineWidth',1.25,'Color','k','MarkerFaceColor','b','MarkerEdgeColor','b');
view(30,30)
axis equal
grid on
axis([-.6 .6 -.6 .6 -.6 .1]);
spatial_title = title(sprintf('Time: %0.2f s',t(1)),'FontName','Helvetica','FontSize',11.5);

mov_3D(1,length(t)) = struct('cdata',[],'colormap',[]);

for i = 1:length(t)
    set(pointer_x(1),'XData',t(i),'YData',x1(i));
    set(pointer_x(2),'Xdata',t(i),'YData',x2(i));
    set(pointer_y(1),'XData',t(i),'YData',y1(i));
    set(pointer_y(2),'Xdata',t(i),'YData',y2(i));
    set(pointer_z(1),'XData',t(i),'YData',z1(i));
    set(pointer_z(2),'Xdata',t(i),'YData',z2(i));
    set(spatial(1),'XData',[0,x1(i)],'YData',[0,y1(i)],'ZData',[0,z1(i)]);
    set(spatial(2),'XData',[x1(i),x2(i)],'YData',[y1(i),y2(i)],'ZData',[z1(i),z2(i)]);
    set(spatial_title,'String',sprintf('Time: %0.2f s',t(i)));
    pause(0.1);
    mov_3D(i) = getframe(gcf);
end

movie = VideoWriter('mov_3D.avi');
open(movie)
writeVideo(movie,mov_3D)
close(movie)
