load('XDataSample.mat');
close all

XDataR = XData(XData(:,2) == 160,:);

B = XDataR(:,1);
T = XDataR(1,2);
t = XDataR(:,3);

ParamBaseline = [250    % c
-200            % psi
3               % nc
0];             % th
% 969.75      -281.33            3       2.3202
figure(1)
% hold on

% Param is a 4x1 matrix of   c  psi  nc   th
alpha = isoAvrami(Param,XDataR);
Param = ParamBaseline;

%%
h(1) = subplot(3,2,1);
hold on
plot(t,isoAvrami(Param,XDataR),'b','LineWidth',1.5)

% 1000
Param(1) = 350;
plot(t,isoAvrami(Param,XDataR),'b:','LineWidth',1.5)

% 10
Param(1) = 150;
plot(t,isoAvrami(Param,XDataR),'b--','LineWidth',1.5)

xlabel('Time (min)')
ylabel('Relative Crystallinity')
grid 
legend('250','350','150','Location','SouthEast')

Param = ParamBaseline;

%%
h(2) = subplot(3,2,2);
hold on


% 1
Param(3) = 1;
plot(t,isoAvrami(Param,XDataR),'b','LineWidth',1.5)

% 2
Param(3) = 2;
plot(t,isoAvrami(Param,XDataR),'b:','LineWidth',1.5)

Param(3) = 3;
plot(t,isoAvrami(Param,XDataR),'r','LineWidth',1.5)
% 4
Param(3) = 4;
plot(t,isoAvrami(Param,XDataR),'r:','LineWidth',1.5)

xlabel('Time (min)')
ylabel('Relative Crystallinity')
grid 
legend('1','2','3','4','Location','SouthEast')

Param = ParamBaseline;

%%
h(3) = subplot(3,2,3);
hold on
plot(t,isoAvrami(Param,XDataR),'b','LineWidth',1.5)

Param(2) = -100;
plot(t,isoAvrami(Param,XDataR),'b--','LineWidth',1.5)

Param(2) = -300;
plot(t,isoAvrami(Param,XDataR),'b:','LineWidth',1.5)

xlabel('Time (min)')
ylabel('Relative Crystallinity')
grid
legend('-200','-100','-300','Location','SouthEast')

Param = ParamBaseline;

%%
h(4) = subplot(3,2,4);
hold on
plot(t,isoAvrami(Param,XDataR),'b','LineWidth',1.5)

Param(4) = 1;
plot(t,isoAvrami(Param,XDataR),'b--','LineWidth',1.5)

Param(4) = 5;
plot(t,isoAvrami(Param,XDataR),'b:','LineWidth',1.5)

xlabel('Time (min)')
ylabel('Relative Crystallinity')
grid
legend('0','1','5','Location','SouthEast')

Param = ParamBaseline;

%%
h(5) = subplot(3,2,5); % the last (odd) axes
hold on
XDataR(:,2) = 150*ones(length(XDataR(:,2)),1);
plot(t,isoAvrami(Param,XDataR),'b--','LineWidth',1.5)

XDataR(:,2) = 160*ones(length(XDataR(:,2)),1);
plot(t,isoAvrami(Param,XDataR),'b','LineWidth',1.5)

XDataR(:,2) = 170*ones(length(XDataR(:,2)),1);
plot(t,isoAvrami(Param,XDataR),'b:','LineWidth',1.5)

xlabel('Time (min)')
ylabel('Relative Crystallinity')
grid 
legend('150 C','160 C','170 C','Location','SouthEast')

Param = ParamBaseline;

%% reposition subplot 5
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)])