%freqz
clear all
close all 
clc
N=2048; %liczba próbek
a=[1 -2.2137 2.9403 -2.1697 0.9606]; % wspó³czynniki modelu AR
y = filter(1,a, (0+1.*randn(1,N)));
plot(y);
title('Próbki modelu AR');
xlabel('N');
ylabel('Amplituda');
axis([0 N -inf inf])
[h1,w1] = freqz(1,a,N);
figure('units','normalized','outerposition',[0 0 1 1])
plot(w1/pi,20*log10(abs(h1)))
title('Widmowa gêstoœæ mocy modelu AR');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');


K=2; %liczba segmentów
sum = zeros(1,(N/K)/2+1);
for i=1:K  %metoda uœrednionego perdiodogramu
    segment = y((i-1)*(N/K)+1:i*(N/K));
    f = fft(segment);
    f = f(1:(N/K)/2+1);
    psdx = (1/(2*pi*(N/K))) * abs(f).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    sum = sum + psdx;
end
sum = sum./K;
freq = 0:(2)/(N/K):1;
figure
plot(freq,10*log10(sum))
title(['Metoda uœrednionego periodogramu; K=', num2str(K)]);
axis([0 1 -inf inf])
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');

K_=K;




korelacja = zeros(1,N);
for i=0:N %estymata autokorelacji
   suma=0;
   for j=0:(N-i-1)
       suma = suma + conj(y(j+1))*y(j+i+1);  
   end
   korelacja(i+1)=suma/N;
end




%algorytm L-D - nie trzeba odwracaæ macierzy autokorelacji.
%Rozwi¹zuje r-a Yula-Walkera i zwraca estymowane wspó³czynniki procesu
order = 4; %rz¹d modelu.
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

[h,w] = freqz([1],a_hat_new,N);
figure
plot(w/pi,20*log10(abs(h)));
title('Parametryczna estymata widmowej gêstoœci mocy Yula-Walkera');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');

figure
plot(w1/pi,20*log10(abs(h1)),freq,10*log10(sum),w/pi,20*log10(abs(h)))
title('Porównanie metod');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');
legend('Rzeczywista',['nieparametryczna K=',num2str(K_)],['parametryczna rz¹d=',num2str(order)])


order = 2; %rz¹d modelu.
a_hat = zeros(1,order+1);
a_hat(1) = 1;
K = -korelacja(2)/korelacja(1);
a_hat(2) = K;
Alpha = korelacja(1)*(1-K*K);
a_hat_new = a_hat;
order_=order;
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

[h,w] = freqz([1],a_hat_new,N);
figure
plot(w/pi,20*log10(abs(h)));
title('Parametryczna estymata widmowej gêstoœci mocy - za ma³y rz¹d');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');

order = 50; %rz¹d modelu.
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

[h1,w1] = freqz([1],a_hat_new,N);
figure
plot(w1/pi,20*log10(abs(h1)));
title('Parametryczna estymata widmowej gêstoœci mocy - za du¿y rz¹d');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');

figure
plot(w1/pi,20*log10(abs(h1)),w/pi,20*log10(abs(h)))
title('Porównanie wp³ywu rzêdu na kszta³t estymaty');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');
legend(['rz¹d = ',num2str(order)],['rz¹d = ',num2str(order_)])

%Zbyt ma³y rz¹d sprawia, ¿e estymata staje siê niedok³adna. Zmniejsza siê 
%jej rozdzielczoœæ oraz wzrasta obci¹¿enie.
%Zbyt du¿y rz¹d estymaty wprowadza efekt zafalowania. Dla rzêdu 100 etymata
%parametryczna przypomina przebiegiem estymatê nieparametryczn¹. Jest
%nieregularna.


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

%Wartoœci FPE i AIC jednoznacznie pokazuj¹, ¿e zwiêkszanie rzêdu modelu
%daje najwiêcej korzyœci do momentu, gdy rz¹d estymatora nie przekracza
%rzêdu modelu estymowanego. Po tym progu wzrost wartoœci obu kryteriów jest
%wci¹¿ zauwa¿alny, ale tempo jego przyrostu jest znacznie mniejsze.




AR =[1 -2.2137 2.9403 -2.1697 0.9606]; %MODEL ARMA
MA = [-0.9 0.61 0.609 -0.899];

Y_arma = filter(MA,AR,(0+1.*randn(1,N)));
figure
plot(Y_arma)
title('Próbki modelu ARMA');
axis([0 N -inf inf])
xlabel('N');
ylabel('Amplituda');


korelacja = zeros(1,N);
for i=0:N
   suma=0;
   for j=0:(N-i-1)
       suma = suma + conj(Y_arma(j+1))*Y_arma(j+i+1);  
   end
   korelacja(i+1)=suma/N;
end

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
title('FPE ARMA')
xlabel('Rz¹d modelu')
figure
plot(1:K_,AIC)
title('AIC ARMA')
xlabel('Rz¹d modelu')


%algorytm L-D 
order = 100; %rz¹d modelu
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

[h4,w4] = freqz(MA,AR,N);
figure
plot(w4/pi,20*log10(abs(h4)));
title('Widmowa gêstoœæ mocy modelu ARMA');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');

[h3,w3] = freqz([1],a_hat_new,N);
figure
plot(w3/pi,20*log10(abs(h3)));
title('Parametryczna estymata widmowej gêstoœci mocy');
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');

%Mimo i¿ wartoœci AIC i FPE wskazuj¹, ¿e optymalny model ma rz¹d równy 2,
%to z obserwacji wykresów wynika, ¿e aby estymata by³a dok³adna nale¿y u¿yæ
%modelu o rzêdzie bliskim 100.
figure
plot(w4/pi,20*log10(abs(h4)),w3/pi,20*log10(abs(h3)))
title('Porównanie widma modelu ARMA z estymat¹ parametryczn¹')
xlabel('Znormalizowana czêstotliwoœæ');
ylabel('Moc/czêstotliwoœæ');
legend('Widmo',['Estymata, rz¹d =',num2str(order)]);


