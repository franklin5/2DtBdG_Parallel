%%%%% plot the data.
%% ground state solution
% clear
% clc
% load superfluid.dat
% h = superfluid(:,1);
% Delta = superfluid(:,2);
% Mu = superfluid(:,3);
% Eg = superfluid(:,4);
% set(gca,'fontsize',16);
% figure(1)
% plot(h, Delta, 'r', h, Mu, '--', h, Eg, 'k','linewidth',2)
% xlabel('h/E_F')
% legend('\Delta','\mu','E_g')
%%
clear
filename = {'hi_0.3hf_1.2.dat',...
        'hi_2.1hf_0.dat'};
idata = 2;
data = load(filename{idata});
t = data(:,1);
Delta = data(:,2) + 1i* data(:,3);
plot(t,abs(Delta))
