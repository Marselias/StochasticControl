clear all
close all
clc

N=2000;
s=randn(1,N);

A = [1 -1.6 0.78 -0.18]; % A(q^-1)
B = [1 0.4 -0.45]; % B(q^-1)
C = [1 0.3]; % C(q^-1)

%z równañ diofantycznych MA

G = [1.9 -0.78 0.18]; % G(q^-1)
F = [1]; % F(q^-1)

y = filter(C,A,s);
plot(y)
title('bez sterowania')

u = filter(-G,B,s);
figure
plot(u)
title('sterowanie')

C_G = [-0.9 1.08 -0.18]; %C-G

y_mv = filter(C_G,A,s);
figure
plot(y_mv)
title('ze sterowaniem')

F_MA = [1 0.75];
varMA = 0;
for i=1:20
s1 = randn(1,N);
y_ma = filter(F_MA,1,s1);
varMA = varMA + var(y_ma); 
end
varMA = varMA/20

B_plus = [1 -0.5];
G_mv = [1.15 -0.615 0.15];
u_ma = filter(-G_mv,conv(B_plus,F_MA),y_ma);
figure
plot(u_ma)
title('Sterowanie MA')
figure
plot(y_ma)
title('Ze sterowaniem MA')

%Autocorrelation
r = zeros(1,N);
r(1) = 1;
r(2) = F_MA(2)/(1+F_MA(2)^2);

figure()
plot(r)
hold on
korelacja = xcorr(y_ma,y_ma);
korelacja =korelacja(ceil(length(korelacja)/2):end);
plot(korelacja./max(korelacja))
title('Porównanie autokorelacji')
legend('Teoretyczna','Estymowana')
hold off

[h1,w1] = freqz(F_MA,1,N/2+1);
figure
plot(w1/pi,20*log10(abs(h1)))
h_sr = zeros(N/2+1,1);
psdx_sr = zeros(1,N/2+1);
for rel=1:20
s1 = randn(1,N);
y_ma = filter(F_MA,1,s1);
korelacja = xcorr(y_ma,y_ma);
korelacja =korelacja(ceil(length(korelacja)/2):end);
order = 5; %rz¹d modelu.
a_hat = zeros(1,order+1);
a_hat(1) = 1;
K = -korelacja(2)/korelacja(1);
a_hat(2) = K;
Alpha = korelacja(1)*(1-K*K);
a_hat_new = a_hat;

for i=2:order
    suma = 0;
    for j=1:i-1
        suma = suma + korelacja(j+1)*a_hat(i-j+1);
    end
    suma = suma + korelacja(i+1);
    K = -suma/Alpha;
    for j=1:i-1
        a_hat_new(j+1)=a_hat(j+1) + K*a_hat(i-j+1);
    end
    a_hat_new(i+1)=K;
    a_hat = a_hat_new;
    Alpha = Alpha*(1-K*K);
end
a_hat_new;

[h,w] = freqz([1],a_hat_new,N/2+1);
 h_sr = h_sr + h;
 f = fft(y_ma);
f = f(1:(N)/2+1);
psdx = (1/(2*pi*N)) * abs(f).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
psdx_sr = psdx_sr + psdx;
end
h_sr = h_sr./20;
psdx_sr = psdx_sr./20;
figure
plot(w/pi,20*log10(abs(h_sr)));
title('Parametryczna estymata widmowej gêstoœci mocy Yula-Walkera');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');


figure
plot(w/pi,10*log10(psdx_sr))

figure
plot(w/pi,20*log10(abs(h1)),w/pi,20*log10(abs(h_sr)),w/pi,10*log10(psdx_sr))
legend('teoretyczna','parametryczna','nieparametryczna')
title('Porównanie metod estymacji widmowej gêstoœci mocy')
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');


K=100; %maksymalny rz¹d
K_=K;
FPE = zeros(1,K);
AIC  = zeros(1,K);
for l=1:K
order = l;
a_hat = zeros(1,order+1);
a_hat(1) = 1;
K = -korelacja(2)/korelacja(1);
a_hat(2) = K;
Alpha = korelacja(1)*(1-K*K);
a_hat_new = a_hat;

for i=2:order
    suma = 0;
    for j=1:i-1
        suma = suma + korelacja(j+1)*a_hat(i-j+1);
    end
    suma = suma + korelacja(i+1);
    K = -suma/Alpha;
    for j=1:i-1
        a_hat_new(j+1)=a_hat(j+1) + K*a_hat(i-j+1);
    end
    a_hat_new(i+1)=K;
    a_hat = a_hat_new;
    Alpha = Alpha*(1-K*K);
end    
    
    
variance = korelacja(1);
for m=2:order+1
    variance = variance - a_hat_new(m)*korelacja(m);
end

FPE(l) = variance*(1+order/N)/(1-order/N);
AIC(l) = N*log(variance) + 2*order;
end
figure
plot(1:K_,FPE)
title('FPE')
xlabel('Rz¹d modelu')
figure
plot(1:K_,AIC)
title('AIC')
xlabel('Rz¹d modelu')

