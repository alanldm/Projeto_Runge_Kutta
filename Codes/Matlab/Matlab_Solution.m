%Limpando o console e o banco de variávies
clear all; clc;

%Definindo as condições iniciais do problema
y0 = -0.24593576; %y(1)
z0 = -0.55769344; %y'(1)

% ======= Calculando os valores exatos ======= %

%Definindo os intervalos para os valores exatos
x1 = 1:0.025:10*pi;
x2 = 1:0.25:10*pi;
x3 = 1:0.5:10*pi;

%Inicializando as variáveis que irão armazenar os valores exatos
exato_j025 = []; %Para h = 0.025
exato_j25 = [];  %Para h = 0.25
exato_j5 = [];   %Para h = 0.5

%Calculando os valores exatos com a função de bessel de primeira espécie e ordem zero
for x = x1
    
    j1 = sqrt(x) * besselj(0, 10*x); 
    exato_j025 = [exato_j025; j1]; %Armazenando o resultado
end
for x = x2
    
    j2 = sqrt(x) * besselj(0, 10*x);
    exato_j25 = [exato_j25; j2]; %Armazenando o resultado
end
for x = x3
    
    j3 = sqrt(x) * besselj(0, 10*x);
    exato_j5 = [exato_j5; j3]; %Armazenando o resultado
end

% ======= Calculando os valores com RK4 ======= %

%Definindo os valores de h
p = [0.025 0.25 0.5];

%Inicializando as variáveis que irão armazenar os valores para RK4
y_025_rk4 = [];
y_25_rk4 = [];
y_5_rk4 = [];


for h = p
    
   x = 1:h:10*pi;
   y = zeros(length(x),1);
   z = zeros(length(x),1);
   
   %Atribuindo as condições iniciais 
   y(1) = y0; 
   z(1) = z0; 
       
   for n = 2:length(x)
   [c1,c2] = coeficientes(h,x(n-1),y(n-1),z(n-1));
   z(n) = z(n-1) + (1/6) * c1;
   y(n) = y(n-1) + (1/6) * c2;
   end
   
    if h == 0.025
       y_025_rk4 = y;
    elseif h == 0.25
        y_25_rk4 = y;
    else
        y_5_rk4 = y;
    end    
 
end

% ======= Calculando os valores com RK6 (Método de Collatz) ======= %

%Inicializando as variáveis que irão armazenar os valores para RK6
y_025_rk6 = [];
y_25_rk6 = [];
y_5_rk6 = [];


for h = p
    
   x = 1:h:10*pi;
   y = zeros(length(x),1);
   z = zeros(length(x),1);
    
   y(1) = y0;
   z(1) = z0;
       
   for n = 2:length(x)
   [t1,t2] = termos(h,x(n-1),y(n-1),z(n-1));
   z(n) = z(n-1) + (1/(90*h))*t1;
   y(n) = y(n-1) + h*z(n-1) + (1/90)*t2;
   end
   
    if h == 0.025
       y_025_rk6 = y;
    elseif h == 0.25
        y_25_rk6 = y;
    else
        y_5_rk6 = y;
    end    
 
end

% ======= Visualizando os resultados obtidos  ======= %

t1 = table(y_025_rk4, y_025_rk6,exato_j025);
t2 = table(y_25_rk4, y_25_rk6,exato_j25);
t3 = table(y_5_rk4, y_5_rk6,exato_j5);
disp(t1);
disp(t2);
disp(t3);

%Função para calcular os coeficientes do RK4
function [c1,c2] = coeficientes(h, x, y, z)

  k1 = h*z;
  l1 = h*-(100+(1/x^2))*y;
  k2 = h*(z+l1/2);
  l2 = h*-(100+(1/(x+h/2)^2))*(y+k1/2); 
  k3 = h*(z+l2/2);
  l3 = h*-(100+(1/(x+h/2)^2))*(y+k2/2);
  k4 = h*(z+l3);
  l4 = h*-(100+(1/(x+h)^2))*(y+k3);

  c1 = l1+2*l2+2*l3+l4;
  c2 = k1+2*k2+2*k3+k4;

end

%Função para calcular os coeficientes do RK6
function [t1, t2] = termos(h,x,y,z)
  k0 = (h^2)*-(100+(1/x^2))*y;
  k1 = (h^2)*-(100+(1/(x+(h/4))^2))*(y+(h*z/4)+(k0/32));
  k2 = (h^2)*-(100+(1/(x+(h/2))^2))*(y+(h*z/2)-(k0/24)+(k1/6));
  k3 = (h^2)*-(100+(1/(x+(3*h/4))^2))*(y+(3*h*z/4)+(3*k0/32)+(k1/8)+(k2/16));
  k4 = (h^2)*-(100+(1/(x+h)^2))*(y+h*z+(3*k1/7)-(k2/14)+(k3/7));

  t1 = 7*k0+32*k1+12*k2+32*k3+7*k4;
  t2 = 7*k0+24*k1+6*k2+8*k3;

end
