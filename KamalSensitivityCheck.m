% Teuwen
A1 = 5.5E5;
Ea1 = 6.86E4;
A2 = 5.5E5;
Ea2 = 5.94E4;
m = 1.55;
n = 1.21;
teuParam = [A1, Ea1, A2, Ea2, m, n];

% Mine
A1 = 1E5;
Ea1 = 7.5E4;%1.7E5;
A2 = 1E5;
Ea2 = 5.15E4;
m = 0.8242;
n = 0.7911;
myParam = [A1, Ea1, A2, Ea2, m, n];

B = [0.001:0.001:0.999]';
T = (150+273).*ones(length(B),1);
x = [B,T];


figure(1)
cla
% Teu vs Me
subplot(2,4,1)
hold on
plot(B,kamal6(teuParam,x,[]))
plot(B,kamal6(myParam,x,[]))

title('Teuwen vs My Parameter Set'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% m
subplot(2,4,2)
hold on
plot(B,kamal6(myParam,x,[]))
myParam(5) = m.*1.1;
plot(B,kamal6(myParam,x,[]),':')
myParam(5) = m.*0.9;
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('m'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% n 
subplot(2,4,3)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T]; %reset params, T
hold on
plot(B,kamal6(myParam,x,[]))
myParam(6) = n.*1.1;
plot(B,kamal6(myParam,x,[]),':')
myParam(6) = n.*0.9;
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('n'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% T
subplot(2,4,4)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
plot(B,kamal6(myParam,x,[]))
x(:,2) = T.*1.1.*ones(length(x(:,2)),1);
plot(B,kamal6(myParam,x,[]),':')
x(:,2) = T.*0.9.*ones(length(x(:,2)),1);
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('T'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% A1
subplot(2,4,5)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
plot(B,kamal6(myParam,x,[]))
myParam(1) = A1.*1.1;
plot(B,kamal6(myParam,x,[]),':')
myParam(1) = A1.*0.9;
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('A1'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% A2
subplot(2,4,6)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
plot(B,kamal6(myParam,x,[]))
myParam(3) = A2.*1.1;
plot(B,kamal6(myParam,x,[]),':')
myParam(3) = A2.*0.9;
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('A2'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% Ea1
subplot(2,4,7)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
plot(B,kamal6(myParam,x,[]))
myParam(2) = Ea1.*1.1;
plot(B,kamal6(myParam,x,[]),':')
myParam(2) = Ea1.*0.9;
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('Ea1'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% Ea2
subplot(2,4,8)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
plot(B,kamal6(myParam,x,[]))
myParam(4) = Ea2.*1.1;
plot(B,kamal6(myParam,x,[]),':')
myParam(4) = Ea2.*0.9;
plot(B,kamal6(myParam,x,[]),':')
legend('Baseline','+10%','-10%','Location','South')

title('Ea2'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')


%%% 
