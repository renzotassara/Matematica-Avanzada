function Proyecto_MATEMATICA_AVANZADA
clear all
clc
% --------------- DATOS ---------------------------
m=64.8; p=100; 
L=420;  EJ=9.45E11; ML=50000;
n=20;                           %Cantidad de grados de libertad
Dx=L/n;                         %Longitud de los n intervalos 
Dt=1e-3;                        %Incremento de tiempo
tf=20;                          %Tiempo de la respuesta
NDt=(tf/Dt)+1;                  %Cantidad de puntos en el intervalo tf
it=0;
for i=0:Dt:tf   %Dimensionando el vector t
    it=it+1;        
    t(it)=i;
end 
dim=it;

% =========================== MATRIZ DE MASA (dimensionamiento)============================
for i=1:n
  for j=1:n
    M(i,i)= m;
  end 
end 
% =========================== MATRIZ K (rigidez) ===========================
               % 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
K=(EJ/(Dx^4))*  [7,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %1 
                -4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %2
                 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %3
                 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %4
                 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %5
                 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %6
                 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %7
                 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0; %8
                 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0,0; %9
                 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0, 0,0; %10
                 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0, 0,0; %11
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0, 0,0; %12
                 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0, 0,0; %13
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0, 0,0; %14
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0, 0,0; %15
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1, 0,0; %16
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1,0; %17
                 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 6,-4, 1  %18
                 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-4, 5,-2; %19
                 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-4, 2];%20
                 

% =========================== IMPULSO UNITARIO ===========================
g=zeros(1,NDt); it=0;
g(1)=1/Dt;
% =========================== VECTOR fp ===========================
for i=1:n
  fp(i,1)=p;
end 
% =========================== VECTOR fm ===========================
fm(n-1,1)=(ML/(Dx^2))*-1;
fm(n,1)=(ML/(Dx^2))*2;
%============================ APLICANDO DIFERENCIA CENTRAL (impulso unitario en u(L,t), gc=0 y 20 gdl) ================
u=zeros(n,1);
u0=u;
for i=1:NDt-1
    u(:,i+1)=((Dt^2)*(inv(M)*(g(i)*fp)))+(inv(M)*(2*M-((Dt^2)*K))*u(:,i))-u0;
    u0=u(:,i);
end

figure(1)
plot(t,u(n,:),'r');
grid on

% =========================== AUTOVALORES Y AUTOVECTORES ===========================
A=inv(M)*K;
[V,D,W]=eig(A);
w2=min(D(n,n));
disp('El valor de la menor frecuencia natural w^2 es:');disp(w2);

phi=(V(:,n));                                    %Autovector asociado a la minima frecuencia
phin(1,:)=(abs(V(:,n)))./(norm(W(:,n),inf));     %Autovector (normalizado por izq) asociado a la minima frecuencia

%B=inv(M.')*(K.');  | Aqui compruebo que con la funcion eig(A) usando la
%[V,D,W]=eig(B);    | matriz W, obtengo el mismo resultado que como el calculo      
%V(:,20)            | con (Kt - w^2*Mt)

% ---------------- VALOR DE CONSTANTES ------------------------------

disp('Los valores de las constantes de los parámetros de 1 grado de libertad son:');
kn=phin*K*phi;Tx=[' kn = ',num2str(kn)];disp(Tx)
mn=phin*M*phi;Tx=[' mn = ',num2str(mn)];disp(Tx)
bp=(phin*fp)/mn;Tx=[' bp = ',num2str(bp)];disp(Tx)
bm=(phin*fm)/mn;Tx=[' bm = ',num2str(bm)];disp(Tx)
w2n=((kn)/(mn));Tx=['(wn)^2 = ',num2str(w2n)];disp(Tx)
T=(2*pi)/(w2n^(1/2));Tx=['Periodo T = ',num2str(T)];disp(Tx)

% ================= APLICANDO DIFERENCIA CENTRAL (impulso unitario en u(L,t), gc=0 y 1 gdl) ================
q=0;
q0=0;
for i=1:NDt-1
    q(i+1)=((Dt^2)*bp*g(i)+(2-(Dt^2)*w2n)*q(i)-q0); 
    q0=q(i);
end

figure(2)
plot(t,q,'r');
grid on


disp('');
qr=input('ingresar valor de qr:  ');


it=1;
for j=1:6 % ------------------------ EL ULTRA FOR -----------------------------

bm=(8.37/348.94); 
bp=(841.92/348.94);  %Cambio el valor de bm y bp para obtener datos coincidentes
disp(' ');

Kd=[30 70  100 140 300 600];
Kp=Kd(j)*3;
Ki=Kd(j)*2; 
C=1;

% ---------------- POLOS DE Q(s)-----------------------------
% Con g(t) impulso unitario

Ph=[1 0 w2n];
Phc=[1 Kp/Kd(j) Ki/Kd(j)];
Pq=[1 C*bm*Kd(j) (w2n)+(C*bm*Kd(j))*3 C*bm*Kd(j)*2];

polosq=roots(Pq);
polosh=roots(Ph);
poloshc=roots(Phc);
rp=real(poloshc);
ip=imag(poloshc);
disp(' ');




Tx=['Para Kd= ',num2str(Kd(j))];disp(Tx);
disp('Los polos de Q(s) son: ');disp(polosq);
disp('Los polos de H(s) son: ');disp(polosh);
disp('Los cero de Hc(s) son: ');disp(poloshc);
disp(' ');

% ======================== DETERMINACION DE q(t) ============

qf=0;
qf0=0;

%--------Posicion--------

for i=1:NDt-1
    qf(i+1)=(Dt^2*(qf(i)*(-w2n-bm*Kp*C)+(g(i)*bp)+bm*Kp*C*qr)+2*qf(i)+(Dt^2)*bm*C*Ki*(Dt*(qr-qf(i)))+((Dt/2)*bm*Kd(j)-1)*qf0)/(1+(Dt/2)*bm*Kd(j)); 
    qf0=qf(i);
end


%-------Fuerza de control con variacion de qr
Int=0;
for i=2:NDt-1
  Int=Int+(Dt*(qr-qf(i)));
  gc(i)=Kp*C*(qr-qf(i))-((Kd(j)*C)/(2*Dt))*(qf(i+1)-qf(i-1))-Ki*C*Int;
end
gc(1)=Kp*C*(qr-qf(1))-((Kd(j)*C)/(2*Dt))*(-3*qf(1)+4*qf(2)-qf(3))-Ki*C*((Dt/2)*(qr-qf(1)));   % Aplico derivada central en los ptos intermedios
gc(NDt)=Kp*C*(qr-qf(NDt))-((Kd(j)*C)/(2*Dt))*(3*qf(NDt)-4*qf(NDt-1)+qf(NDt-2))-Ki*C*(((Dt/2)*(qr-qf(NDt))+Int));  % y hacia delante o hacia atras para los ptos extremos

if it==1
  gcc=gc;
end


% ------------- U(L,t) EN 20 GRADOS DE LIBERTAD --------------
u=zeros(n,1);
u0=u;
qf(0)=0;
qf0=0;
for i=1:NDt-1
    qf(i+1)=(Dt^2*(qf(i)*(-w2n-bm*Kp*C)+(g(i)*bp)+bm*Kp*C*qr)+2*qf(i)+(Dt^2)*bm*C*Ki*(Dt*(qr-qf(i)))+((Dt/2)*bm*Kd(j)-1)*qf0)/(1+(Dt/2)*bm*Kd(j))
    Int=Int+(Dt*(qr-qf(i)));
    gc(1)=Kp*C*(qr-qf(1))-((Kd(j)*C)/(2*Dt))*(-3*qf(1)+4*qf(2)-qf(3))-Ki*C*((Dt/2)*(qr-qf(1)));
    gc(i)=Kp*C*(qr-qf(i))-((Kd(j)*C)/(2*Dt))*(qf(i+1)-qf(i-1))-Ki*C*Int;
    u(:,i+1)=((Dt^2)*inv(M))*((g(i)*fp)+gc(i)*fm-K*u(:,i))+2*u(:,i)-u0;
    qf0=qf(i);
    u0=u(:,i);
end


% -------Tiempo objetivo-----------

Tobj=20;
T=1;
while Tobj>1.4*(10^(-6))
    T=T+1;
    Tobj=abs((1/(2*Dt))*(qf(T+1)-qf(T-1)));
end

Tobjj(j)=T*Dt;
Normgc(j)=max(abs(gc));

% ----------------- GAFICOS -------------------





Kdd=Kd(j);

figure(3)
subplot(3,2,it);
plot(t,qf),title(['Modelo de un grado de libertad con qr=0 y Kd = ', num2str(Kdd)])
grid on

figure(4)
subplot(3,2,it);
plot(t,u(n,:)),title(['Modelo de 20 grados de libertad para Kd = ', num2str(Kdd)])
grid on

figure(5)
subplot(3,2,it);
plot(t,abs(gc)),title(['Fuerza de control con Kd = ', num2str(Kdd)])
grid on

figure(6)
subplot(3,2,it);
hold on
plot(polosq,'bo')
plot(polosh,'g*')
plot(rp,ip,'ks')
title(['Polos de Q(s) con Kd = ', num2str(Kdd)])
grid on

it=it+1;
end



figure(7)
plot(Tobjj,Kd),title('Tiempo objetivo en funcion de Kd')
grid on
figure(8)
plot(Normgc,Kd),title('Norma infinito de gc(t) en funcion de Kd')
grid on
TODO=[Kd;Tobjj;Normgc]
% -------------- MINIMOS CUADRADOS ------------------------------- 

P=3;
N=length(Kd);
  for i=1:N
    for j=1:P
      phii(i,j)=(1/Kd(i))^(j-1);
    end 
  end 

va=inv(phii'*phii)*(phii'*Tobjj');

alfa=va(1);Tx=[' Alfa = ',num2str(alfa)];disp(Tx)
beta=va(2);Tx=[' Beta = ',num2str(beta)];disp(Tx)
gamma=va(3);Tx=[' Gamma = ',num2str(gamma)];disp(Tx)




% %  ------------------ PARTE 3 -------------------------
%ACCION DADA POR UNA FUNCION ESCALON


%-------gc(t)=0
u=zeros(n,1);
u0=u;
for i=1:NDt-1
    u(:,i+1)=((Dt^2)*(inv(M)*(1*fp)))+(inv(M)*(2*M-((Dt^2)*K))*u(:,i))-u0;
    u0=u(:,i);
end

figure(9)
plot(t,u(n,:),'r');
grid on

%-------gc(t)no nula
u=zeros(n,1);
u0=u;

for i=1:NDt-1
    u(:,i+1)=((Dt^2)*(inv(M)*(1*fp)))+((Dt^2)*(inv(M)*(gcc(i))*fm))+(inv(M)*(2*M-((Dt^2)*K))*u(:,i))-u0;
    u0=u(:,i);
end
% Obtenemos el modelo de 20 grados de libertad respondiendo ante el impulso unitario
% para gc(t) con qr=0

figure (10)
plot(t,u(n,:)),title('20 grados de libertad para g(t) como funcion escalon')
grid on

%ACCION DADA POR UN IMPULSO TRIANGULAR
gt=zeros(1,NDt);
gt(1,500)=1;
gt(1,1000)=0;
for i=1:1000
    if i<500
     gt(1,i)=(2*i)/1000;
    end
    if i>500
     gt(1,i)=abs((2*i-2000)/1000);
    end
end

%-------gc(t)=0
u=zeros(n,1);
u0=u;
for i=1:NDt-1
    u(:,i+1)=((Dt^2)*(inv(M)*(gt(1,i)*fp)))+(inv(M)*(2*M-((Dt^2)*K))*u(:,i))-u0;
    u0=u(:,i);
end

figure(11)
plot(t,u(n,:),'r');
grid on

%-------gc(t)no nula
u=zeros(n,1);
u0=u;

for i=1:NDt-1
    u(:,i+1)=((Dt^2)*(inv(M)*(gt(1,i)*fp)))+((Dt^2)*(inv(M)*(gc(i))*fm))+(inv(M)*(2*M-((Dt^2)*K))*u(:,i))-u0;
    u0=u(:,i);
end
% Obtenemos el modelo de 20 grados de libertad respondiendo ante el impulso unitario
% para gc(t) con qr=0

figure (12)
plot(t,u(n,:))
grid on






  
  
  end 