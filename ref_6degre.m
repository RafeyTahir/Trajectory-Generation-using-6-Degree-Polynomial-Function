clc
clear all
% Define the polynomial function
f = @(t, a_6, a_5, a_4, a_3, a_2, a_1, a_0) a_6*t^6 + a_5*t^5 + a_4*t^4 + a_3*t^3 + a_2*t^2 + a_1*t + a_0;

% Set the time interval and step size
t_min = 0;
t_max = 1;
h = 0.01;
p=0.5;

% Set the polynomial coefficients
a_6 = (2*p-1)*(6*p^4 -12*p^3 + 4*p^2 +2*p +1)/(p^3 *(p -1)^3);
a_5 = -3*(4*p^6 -18*p^4 +16*p^3 -1)/(p^3 *(p -1)^3);
a_4 = 3*(10*p^6 -18*p^5 +10*p^3 -1)/(p^3 *(p -1)^3);
a_3 = -(20*p^6 -48*p^5 +30*p^4 -1)/(p^3 *(p-1)^3);
a_2 = 0;
a_1 = 0;
a_0 = 1;

% Initialize the time and reference trajectory arrays
t = t_min:h:t_max;
ref = zeros(size(t));

% Generate the reference trajectory
for i = 1:length(t)
ref(i) = f(t(i), a_6, a_5, a_4, a_3, a_2, a_1, a_0);
end

p=0.65;

% Set the polynomial coefficients
a_6 = (2*p-1)*(6*p^4 -12*p^3 + 4*p^2 +2*p +1)/(p^3 *(p -1)^3);
a_5 = -3*(4*p^6 -18*p^4 +16*p^3 -1)/(p^3 *(p -1)^3);
a_4 = 3*(10*p^6 -18*p^5 +10*p^3 -1)/(p^3 *(p -1)^3);
a_3 = -(20*p^6 -48*p^5 +30*p^4 -1)/(p^3 *(p-1)^3);
a_2 = 0;
a_1 = 0;
a_0 = 1;

refft = zeros(size(t));

% Generate the reference trajectory
for i = 1:length(t)
refft(i) = f(t(i), a_6, a_5, a_4, a_3, a_2, a_1, a_0);
end

p=0.58;

% Set the polynomial coefficients
a_6 = (2*p-1)*(6*p^4 -12*p^3 + 4*p^2 +2*p +1)/(p^3 *(p -1)^3);
a_5 = -3*(4*p^6 -18*p^4 +16*p^3 -1)/(p^3 *(p -1)^3);
a_4 = 3*(10*p^6 -18*p^5 +10*p^3 -1)/(p^3 *(p -1)^3);
a_3 = -(20*p^6 -48*p^5 +30*p^4 -1)/(p^3 *(p-1)^3);
a_2 = 0;
a_1 = 0;
a_0 = 1;

reff = zeros(size(t));

% Generate the reference trajectory
for i = 1:length(t)
reff(i) = f(t(i), a_6, a_5, a_4, a_3, a_2, a_1, a_0);
end

p=0.42;

% Set the polynomial coefficients
a_6 = (2*p-1)*(6*p^4 -12*p^3 + 4*p^2 +2*p +1)/(p^3 *(p -1)^3);
a_5 = -3*(4*p^6 -18*p^4 +16*p^3 -1)/(p^3 *(p -1)^3);
a_4 = 3*(10*p^6 -18*p^5 +10*p^3 -1)/(p^3 *(p -1)^3);
a_3 = -(20*p^6 -48*p^5 +30*p^4 -1)/(p^3 *(p-1)^3);
a_2 = 0;
a_1 = 0;
a_0 = 1;

refff = zeros(size(t));

% Generate the reference trajectory
for i = 1:length(t)
refff(i) = f(t(i), a_6, a_5, a_4, a_3, a_2, a_1, a_0);
end


p=0.35;

% Set the polynomial coefficients
a_6 = (2*p-1)*(6*p^4 -12*p^3 + 4*p^2 +2*p +1)/(p^3 *(p -1)^3);
a_5 = -3*(4*p^6 -18*p^4 +16*p^3 -1)/(p^3 *(p -1)^3);
a_4 = 3*(10*p^6 -18*p^5 +10*p^3 -1)/(p^3 *(p -1)^3);
a_3 = -(20*p^6 -48*p^5 +30*p^4 -1)/(p^3 *(p-1)^3);
a_2 = 0;
a_1 = 0;
a_0 = 1;

reffft = zeros(size(t));

% Generate the reference trajectory
for i = 1:length(t)
reffft(i) = f(t(i), a_6, a_5, a_4, a_3, a_2, a_1, a_0);
end





% Plot the reference trajectory
plot(t, ref,'--K')
hold on
plot(t,reff,'--b')
hold on
plot(t,refff,'--g')
plot(t,reffft,'.r')
plot(t,refft,'--y')
xlabel('time(s)')
ylabel('Amplitude')
legend({'$b=0.5$','$b=max$','$b=min$','$b<min$','$b>max$'},'interpreter','latex','Location','best')
print('-dpng', '-r300', 'poly');