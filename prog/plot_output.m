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
close all
filename = {'hi_0.3hf_1.2.dat',...
        'hi_2.1hf_0.dat',...
        'hi_0.3hf_0.1.dat',...
        'hi_2.1hf_0.9.dat'};
file_r = {'hi_0.3hf_1.2_Delta_K_r.dat',...
            'hi_2.1hf_0_Delta_K_r.dat',...
            'hi_0.3hf_0.1_Delta_K_r.dat',...
            'hi_2.1hf_0.9_Delta_K_r.dat'};
file_i = {'hi_0.3hf_1.2_Delta_K_i.dat',...
            'hi_2.1hf_0_Delta_K_i.dat',...
            'hi_0.3hf_0.1_Delta_K_i.dat',...
            'hi_2.1hf_0.9_Delta_K_i.dat'};
idata = 1;
data = load(filename{idata});
figure(idata)
t = data(:,1);
Delta = data(:,2) + 1i* data(:,3);
plot(t,abs(Delta))
%%
Delta_K = load(file_r{idata}) +1i*load(file_i{idata}) ;
nk = length(Delta_K(1,:));
kx = load('akx.OUT');nkx =length(kx);
ky = load('aky.OUT');nky =length(ky);
for nt = 1:10:length(t)
     temp = reshape(Delta_K(nt,:),nkx,nky);
     figure(nt)
     subplot(2,1,1)
     imagesc(kx,ky,angle(temp))
     axis([-5 5 -5 5])
     subplot(2,1,2)
     imagesc(kx,ky,abs(temp))
     axis([-5 5 -5 5])
end
for nt = 1:10:length(t)
    temp = reshape(Delta_K(nt,:),nkx,nky);
    figure(nt)
    
    subplot(2,1,1)
    
    plot(kx,angle(temp(:,nky/2))/pi)
    ylabel('\theta(k)/\pi')
    
    subplot(2,1,2)
    plot(kx,abs(temp(:,nky/2)))
    xlabel('k_x/k_F')
    ylabel('|\Delta(k)|/\pi')
    title(['phase I, t=',num2str(data(nt,1))])
    saveas(figure(nt), [filename{idata},'t=',num2str(data(nt,1)),'.eps'],'epsc')
end


