init_seperations = [250]; 
Cooupling_ratio = 0.95; CavityLength = 1.45;
PsatTR =1.8; loss = 0.2+log(Cooupling_ratio)/CavityLength;
Omega = sqrt(12); del = 0.02; sigma = 0.01;
g0=0.65;    %Unsaturated gain
rng(1)
gam = 4; beta2 = -2;
PULSE_POWER_THREASHOLD = 0.1;
WIDTH_THREASHOLD = 50;
%Grid Parameters
Nt = 2^12; T = 2^10; dt = T/Nt;
t = (-Nt/2:1:Nt/2-1)'*dt;
dw = 2*pi/T; w = [0:Nt/2-1  -Nt/2:-1]'*dw;
Z = 4000000; h = 0.05; NumSteps = round(Z/h);
SaveDistance = 14500; SaveInterval = round(SaveDistance/h);
raundInterval = round(CavityLength/h);
NOISE_AMP = 0.01; gd = 0.00005; Tr = 100000;
%Operators
separation = [];
L = (1i*beta2*w.^2-loss)/2;
K = (1 - (w/Omega).^2)/2;
%Initial condition
u0 = 1*exp(-((t-init_seperations(1)/2)/2.5).^2);
u0 = 1*exp(-((t-init_seperations(1)/2)/1).^2) +1*exp(-((t+init_seperations(1)/2)/1).^2);
pulse_index_0 = find(t>-init_seperations(1)/2,1);
shift = [];
averaged_spectrum = u0;
Etot = norm(u0).^2*dt;
uf = fft(u0); uplot = abs(u0).'; uplot_w = abs(fftshift(uf)).'; gtplot= [-gd*cumsum(abs(u0).^2)*dt+Etot*gd/Tr*dt.*(1:length(u0))']';
zplot = 0; Psatf = PsatTR/dt*Nt;
% uplot_w_for_avg = zeros(SaveDistance/raundDistance, length(fftshift(uf)));
% uplot_w_avg_sq = [abs(fftshift(uf).').^2];
r=1;
for istep = 1:NumSteps
  uf = uf + sqrt((h/0.05)*(dt/(1/4)))*(1/sqrt(2)*(randn(1,length(uf))'*NOISE_AMP+1i*randn(1,length(uf))'*NOISE_AMP));
  if mod(istep,raundInterval)==0
     uf = uf*Cooupling_ratio;
%        if mod(istep,raundInterval)==0
% %            uplot_w_for_avg(r,:) = (fftshift(uf)).';
%        end
     r=r+1;
  end
  if istep == 25000000
      g0 = 0.58;
  end
  g1 = g0/(1+norm(uf)^2/Psatf);   %Saturated gain
  g2 = -2*(g1^2/g0/Psatf)*...
     real(dot(uf,(L + g1*K).*uf));    %frequency dependence
  u=ifft(exp(L*h/2+(g1*h/2+g2/8*h^2)*K).*uf);
  Etot =norm(u).^2*dt;
  Etot_w_noise = Etot;
  [pks,pulse_indexes,width,promince] = findpeaks(abs(u).^2,'MinPeakProminence',PULSE_POWER_THREASHOLD);
  if length(pulse_indexes)>=2
      current_separation = (abs(pulse_indexes(2)-pulse_indexes(1))*dt);
  end
  if length(pulse_indexes)>=2
    
      indeces = (1:length(u));
      dist_from_pulses = min(abs(indeces - pulse_indexes));
      max_width = max(width);
      noise_avg_pow = mean(abs(u(dist_from_pulses>(WIDTH_THREASHOLD*max_width))).^2);
      Etot_w_noise = Etot + noise_avg_pow*dt*length(u)*(Tr-T)/T;
%         Etot_w_noise = Etot;
      cavity_indexes = (Tr/2/dt-length(u)/2:Tr/2/dt+length(u)/2-1)';
      gt = -gd*cumsum(abs(u).^2)*dt - gd*noise_avg_pow*dt*(Tr/2/dt-length(u)/2)+...
        (Etot_w_noise)*gd/Tr*dt.*cavity_indexes;
  else
      current_separation = 0;
          cavity_indexes = (Tr/2/dt-length(u)/2:Tr/2/dt+length(u)/2-1)';
      gt = -gd*cumsum(abs(u).^2)*dt - Etot*gd/Tr*dt.*cavity_indexes;
  end
  gt = gt - mean(gt)*T/Tr ;
  uf = fft(exp(+(del+1i*gam)*h*abs(u).^2-sigma*h*abs(u).^4+gt.*h).*u);
  g1 = g0/(1+norm(uf)^2/Psatf);
  uf = exp(L*h/2 + (g1*h/2+g2/8*h^2)*K).*uf;
  if mod(istep,SaveInterval)==0
     r=1;
     uplot = [uplot; abs(ifft(uf)).'];
     zplot = [zplot, istep*h];
     uplot_w = [uplot_w; (fftshift(uf)).'];
     gtplot = [gtplot; gt'];
     separation = [separation, current_separation];
%        uplot_w_avg_sq = [uplot_w_avg_sq; mean(abs(uplot_w_for_avg).^2)];
     figure(1)
     subplot(4,2,1)
     mesh(t,zplot/(2e+5),gtplot);
     view(0,90); xlabel('Time (ps)'); ylabel('Time (s)');
     title('Gain')
     colormap jet
     subplot(4,2,3)
     plot(t,abs(uplot(size(uplot,1),:)).^2,t,abs(uplot(1,:)).^2) ; xlabel('Time (ps)'); ylabel('|A(z,t)|^2');
     grid on
     xlim([-20 20])
     subplot(4,2,2)
     mesh(ifftshift(w)/2/pi,zplot,abs(uplot_w).^2);
     xlim([-0.5 0.5])
     view(0,90); xlabel('Time (ps)'); ylabel('z (km)');
     title('|A(z,t)|^2')
     subplot(4,2,5)
     semilogy(t,abs(uplot(size(uplot,1),:)).^2); xlabel('Time (ps)'); ylabel('|A(z,t)|^2');
     ylim([-10,1])
     subplot(4,2,6)
     plot(ifftshift(w)/2/pi,mean(abs(uplot_w).^2));
     subplot(4,2,7)
     plot(t,exp(raundInterval/2*(-loss*h+g1*h+g2/2*h^2+gt*h))*Cooupling_ratio) ; xlabel('Time (ps)'); ylabel('|A(z,t)|^2');
     ylabel('loss ')
     xlabel('time (ps)')
     subplot(4,2,8)
     plot(separation,zplot(1:length(zplot)-1),'linewidth',2,'color','blue')
     ylabel('z (km)')
     xlabel('Separation (ps)')
     if length(pulse_indexes)>0
         uf = fft(circshift(u,pulse_index_0-pulse_indexes(1)));
         shift = [shift pulse_index_0-pulse_indexes(1)];
     end
     istep
  pause(0.01)
%     uplot_w_for_avg = []; uplot_for_avg = [];
  end
  if mod(istep,25e6) == 0
      filename = "NMI_simulation_data";
      save(filename + num2str(istep/25e6))
      uplot = [abs(ifft(uf)).'];
      zplot = istep*h;
      uplot = abs(u).'; uplot_w = abs(fftshift(uf)).';
      gtplot = [gt'];
      separation = [];
%         uplot_w_for_avg = zeros(SaveDistance/raundDistance, length(fftshift(uf)));
%         uplot_w_avg_sq = [abs(fftshift(uf).').^2];
      r=1;
    
  end
end
figure(5)
plot(separation,zplot(1:length(zplot)-1),'linewidth',2,'color','blue')
ylabel('z (km)')
xlabel('Separation (ps)')
filename = "NMI_Simulation_data_final";
save(filename)
