tic
 
dt = 1.4;% time step
n_t = 7337;% number of time steps
photo = importfile_hhg_data("myMPI.o107448");% data file
photo = photo';
t2=linspace(-2067.1*2.0,(n_t-1)*dt-2067.1*2.0,n_t);

photox = zeros(1,length(t2));
for i = 1:length(t2)
    photox(i) = photo(i);
end

t3=t2*24.1888*0.001;

%----------------------------------------%

photos_x1=photox;
window=(tanh((t2+2067.1*2.0)*0.01)-tanh((t2-2067.1*2.0)*0.01))/2;% window function (Phys. Rev. B 106, 205422)
photows_x1 = photos_x1 .* window;
w = linspace(-0.2,0.2,4000);
Iwx1 = zeros(1,length(w));
for j = 1:length(w)
for i = 1+1:length(t2)-1
Iwx1(j) = Iwx1(j)+(dt*photows_x1(i)*exp(1i*w(j)*t2(i)));
end
Iwx1(j) = Iwx1(j) + (dt*photows_x1(1)*exp(1i*w(j)*t2(1)))/2.0;
Iwx1(j) = Iwx1(j) + (dt*photows_x1(length(t2))*exp(1i*w(j)*t2(length(t2))))/2.0;
end
for i = 1:length(w)
Iwx1(i) = abs(Iwx1(i)).*abs(Iwx1(i)).*w(i).*w(i);
end
figure
semilogy(w*27.2114/0.26,Iwx1,'r','LineWidth',2)
xlim([0 10.0])
set(gca,'FontSize',16)
set(gca,'linewidth',2)

toc

