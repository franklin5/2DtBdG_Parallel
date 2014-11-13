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
idata = 1;
data = load(filename{idata});
figure(idata)
t = data(:,1);
Delta = data(:,2) + 1i* data(:,3);
plot(t,abs(Delta))

%%
clear
Delta_K = load('hi_0.3hf_1.2_Delta_K_r.dat') +1i*load('hi_0.3hf_1.2_Delta_K_i.dat') ;
t = Delta_K(:,1);
nk = length(Delta_K(1,:));
kx = load('akx.OUT');nkx =length(kx);
ky = load('aky.OUT');nky =length(ky);
for nt = 1:10:length(t)
     temp = reshape(Delta_K(nt,:),nkx,nky);
     figure(nt)
     subplot(1,2,1)
     imagesc(kx,ky,angle(temp))
     axis([-5 5 -5 5])
     subplot(1,2,2)
     imagesc(kx,ky,abs(temp))
     axis([-5 5 -5 5])
end

