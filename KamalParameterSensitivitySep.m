% Teuwen
A1 = 5.5E5;
Ea1 = 6.86E4;
A2 = 5.5E5;
Ea2 = 5.94E4;
m = 1.55;
n = 1.21;
teuParam = [A1, Ea1, A2, Ea2, m, n];

% % Mine
% A1 = 1E5;
% Ea1 = 7.5E4;%1.7E5;
% A2 = 1E5;
% Ea2 = 5.15E4;
% m = 0.8242;
% n = 0.7911;
myParam = [A1, Ea1, A2, Ea2, m, n];

B = [0.000:0.001:0.999]';
T = (150+273).*ones(length(B),1);
x = [B,T];




%%
% Override myParam
myParam = teuParam;
%%
figure(2)
cla
% Teu vs Me
subplot(2,4,1)
hold on
% plot(B,kamal6separate(teuParam,x))
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');

title('Unaltered'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

legendfigureStringArray ={'Baseline - nth','Baseline - Auto','+10%','-10%'};


% m
subplot(2,4,2)
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');

myParam(5) = m.*1.1;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');

myParam(5) = m.*0.9;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');

legend(legendfigureStringArray)

title('m'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% n 
subplot(2,4,3)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T]; %reset params, T
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');
myParam(6) = n.*1.1;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');

myParam(6) = n.*0.9;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');
% legend('Baseline','+10%','-10%','Location','South')

title('n'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% T
subplot(2,4,4)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');
x(:,2) = T.*1.1.*ones(length(x(:,2)),1);
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');
x(:,2) = T.*0.9.*ones(length(x(:,2)),1);
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');
% legend('Baseline','+10%','-10%','Location','NorthEast')

title('T'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% A1
subplot(2,4,5)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');
myParam(1) = A1.*1.1;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');
myParam(1) = A1.*0.9;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');
% legend('Baseline','+10%','-10%','Location','South')

title('A1'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% A2
subplot(2,4,6)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');
myParam(3) = A2.*1.1;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');
myParam(3) = A2.*0.9;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');
% legend('Baseline','+10%','-10%','Location','South')

title('A2'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% Ea1
subplot(2,4,7)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');
myParam(2) = Ea1.*1.1;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');
myParam(2) = Ea1.*0.9;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');
% legend('Baseline','+10%','-10%','Location','South')

title('Ea1'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')

% Ea2
subplot(2,4,8)

myParam = [A1, Ea1, A2, Ea2, m, n]; x =[B,T];
hold on
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r'); plot(B,b,'b');
myParam(4) = Ea2.*1.1;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r:'); plot(B,b,'b:');
myParam(4) = Ea2.*0.9;
[a,b] = kamal6separate(myParam,x);
plot(B,a,'r--'); plot(B,b,'b--');
% legend('Baseline','+10%','-10%','Location','South')

title('Ea2'); ylabel('dBdt (s^-1)'); xlabel('Conversion (B)')


%%% 
