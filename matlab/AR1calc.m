E = 1500;
V = 300^2;
d = 1;
a_phys = 1;

phi = exp(-d*a_phys/sqrt(V))
c = (1-phi)*E
sigma2 = (1-phi^2)*V

x = 1:1:1000;
y = E;
for i=2: 1000
    y(i) = c + y(i-1) *  phi + normrnd(0,sqrt(sigma2));
end

close all
figure
plot(x,y)