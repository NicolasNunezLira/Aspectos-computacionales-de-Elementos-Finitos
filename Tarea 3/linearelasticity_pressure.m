% Resolucion del problema de valores de contorno de elasticidad lineal
%  /  \nabla\cdot \sigma(u)=-f en \Omega
% <   u=u_D en \Gamma_D
%  \  \sigma\cdot n=g en \Gamma_N
% En su formulacion mixta desplazamiento lineal-presion con Mef estabilizado P1-P1

clear all;close all;clc
%% Configuracion
d=3; %Dimension
Neu=[1]; %Identificador Newmann
Dir=[2]; %Identificador Dirichlet

%% Datos del problema 

%   Calculo de los coeficientes de Lame:
%   (el material que se considerara sera el acero estructural
% Modulo elastico (E) y modulo de poisson (nu)
E=2.1e5; nu=0.26;
% Coeficientes de Lame (lambda, mu)
lambda=(E*nu)/((1+nu)*(1-2*nu));
mu=E/(2*(1+nu));
eps=1/lambda;

%Parametro estabilizacion
alf=1;

% Funciones involucradas
f = @(x,y,z) [(x-x);(y-y);-9.8*x.^0.*y.^0.*z.^0];
ud = @(x,y,z) [(x-x);(y-y);(z-z)];
g = {@(x,y,z) [(x-x);(y-y);(z-z)], @(x,y,z) [(x-x);(y-y);1000*x.^0.*y.^0.*z.^0]};

%% Creacion del mallado
[archivo,ruta]=uigetfile('*.poly','Seleccione archivo de malla');%Lectura archivo de malla
name=archivo(1:end-5);%Rescatar el nombre del .poly

mmax=.01;mmax=num2str(mmax);%Medida maxima como string
c=strcat('./tetgen -pqa',mmax,{' '},archivo);%Comando para crear mallado
system(c{1});%Ejecucion comando anterior

fid=fopen(strcat(ruta,[name,'.1.node']),'r');%Lectura nodos
nid = fscanf(fid,'%i',4);
mesh.p = fscanf(fid,'%f',[nid(3)+(nid(4)+1+d) nid(1)]);
mesh.p = (mesh.p(2:end,:))';
fclose(fid);

fid=fopen(strcat(ruta,[name,'.1.ele']),'r');%Lectura elementos
nele = fscanf(fid,'%i',3);
mesh.t = fscanf(fid,'%f',[nele(3)+(2+d) nele(1)]);
mesh.t = (mesh.t(2:end,:))';
fclose(fid);

fid=fopen(strcat(ruta,[name,'.1.face']),'r');%Lectura caras
nface = fscanf(fid,'%i',2);
mesh.f = fscanf(fid,'%f',[1+d+nface(2) nface(1)]);
mesh.f = (mesh.f(2:end,:))';
fclose(fid);

%% Base de espacio de funciones discreto tetraedro de referencia
%Base polinomios v=[P]*[v]
p=@(x,y,z) [1-x-y-z,x,y,z];
P=@(x,y,z) blkdiag(p(x,y,z),p(x,y,z),p(x,y,z));
%Base gradiente \nabla u=[DP]*[u]
dp=[-1 -1 -1;1 0 0;0 1 0;0 0 1]';
DP=blkdiag(dp,dp,dp);
z=[1,0,0,0,1,0,0,0,1];%Divergencia \nabla\cdot u=[z]*[DP]*[u]
%Base tensorial \varepsilon(u)=[D]*[DP]*[u]
D=[1 0 0 0 0 0 0 0 0; %\partial_1 u_1
    0 0 0 0 1 0 0 0 0; %\partial_2 u_2
    0 0 0 0 0 0 0 0 1; %\partial_3 u_3
    0 1 0 1 0 0 0 0 0; %2*(\partial_2 u_1 + \partial_1 u_2)
    0 0 1 0 0 0 1 0 0; %2*(\partial_3 u_1 + \partial_1 u_3)
    0 0 0 0 0 1 0 1 0]; %2*(\partial_3 u_2 + \partial_2 u_3)

%% Ensamble global 
nnodos=size(mesh.p,1);
nelem=size(mesh.t,1);
K=sparse((d+1)*nnodos,(d+1)*nnodos);
F=sparse((d+1)*nnodos,1);
I=speye(9);

% Matrizes de masa
% \int P'*P
m=@(x,y,z) P(x,y,z)'*P(x,y,z);
M=gauss(m,d);%[M]_{9\times 9}
% \int p'*p
m=@(x,y,z)p(x,y,z)'*p(x,y,z);
m=gauss(m,d);%[m]_{3\times 3}
% \int p
pp=gauss(p,3);
% \int P
PP=gauss(P,d);

% Ensamble global
for i=1:nelem
    %Identificar elemento i
    ind=mesh.t(i,:);
    xx=(mesh.p(ind,1:d))';
    %Cambio de variable
    b=xx(:,2:end)-xx(:,1);
    Jf=abs(det(b));
    B=blkdiag(b,b,b);
    bb=xx(:,1);
    %Diametro elemento i
    h=0;
    for j=1:d+1
        for k=2:d+1
            if norm(xx(:,j)-xx(:,k))>h
                h=norm(xx(:,j)-xx(:,k));
            end
        end
    end
    h=h^2;
    
     %% Matriz de rigidez
    
    %Matriz \int \varepsilon(u):\varepsilon(v) 
    K1=2*mu*(Jf/6)*DP'*(B\I)*D'*D*(B'\DP);%VOl(elem referencia)=1/6;
    % Matriz \int (\nabla\cdot v)p
    K2=-Jf*DP'*(B\I)*z'*pp;
    % Matriz \int (\nabla\cdot u)q
    K3=-Jf*pp'*z*(B'\DP);
    % Matriz \int pq
    K4=-eps*Jf*m;
    % Matriz \int \nabla p \nabla q ---- Termino estabilizacion
    K5= -alf*h*(Jf/6)*dp'*(b\I(1:3,1:3))*(b'\dp);
    %Ensamble K
    Ind=[mesh.t(i,:),mesh.t(i,:)+nnodos,mesh.t(i,:)+2*nnodos,mesh.t(i,:)+3*nnodos];
    K(Ind,Ind)=K(Ind,Ind)+[K1,K2;K3,K4+K5];
    
    %% Vector de carga
    
    %Vector \int f\cdot v
    Ind=Ind(1:end-4);
    F(Ind)=F(Ind)+Jf*M*f(xx(1,:)',xx(2,:)',xx(3,:)');
    %Vector \int \varepsilon(u_d):\varepsilon(v)
    F(Ind)=F(Ind)-2*mu*(Jf/6)*DP'*(B\I)*D'*D*(B'\DP)*ud(xx(1,:)',xx(2,:)',xx(3,:)');
    %Vector \int q(\nabla\cdot u_D)
    ind=mesh.t(i,:)+3*nnodos;
    F(ind)=F(ind)-Jf*pp'*z*(B'\DP)*ud(xx(1,:)',xx(2,:)',xx(3,:)');
    %Vector \int f\cdot (\nabla p) ----- Termino de estabilizacion
    F(ind)=F(ind)-alf*h*dp'*(b\I(1:3,1:3))*PP*f(xx(1,:)',xx(2,:)',xx(3,:)');
end    

%% Imposicion de condiciones de frontera

%Condiciones Dirichlet
u=zeros((d+1)*nnodos,1);
nodirich=1:(d+1)*nnodos;
for i=1:length(Dir)
    dirich=unique(mesh.f(mesh.f(:,end) == Dir(i),1:3));
    Dirich=[dirich ; dirich+nnodos ;  dirich+2*nnodos];
    u(Dirich)=ud(mesh.p(dirich,1),mesh.p(dirich,2),mesh.p(dirich,3));
    nodirich=setdiff(nodirich,Dirich);
end

%Condiciones de Neumann
for i=1:length(Neu)
    gn=g{i};
    neu=mesh.f(find(mesh.f(:,end)==Neu(i)),1:3);
    for j=1:length(neu)
        ind=find(sum((ismember(mesh.t,neu(j,:)))')==3);%Indice elemento de la cara actual
        ind2=find(ismember(mesh.t(ind,:),neu(j,:)));
        ind3=mesh.t(ind,ind2);%Indices nodo de la cara actual
        ind4=find(~ismember(mesh.t(ind,:),neu(j,:)));%Nodo interior elemento actual
        %Nodos cara actual
        a=mesh.p(ind3,1:end-1)';
        %Vector normal cara actual
        n=norm(cross(a(:,2)-a(:,1),a(:,3)-a(:,1)));
        %Parametrizacion
        switch ind4%Composicion [P\circ r]
            case 1
                R=@(s,t) P(s,t,1-s-t)'*P(s,t,1-s-t);
                R=gauss(R,d-1);
            case 2
                R=@(s,t) P(0,s,t)'*P(0,s,t);
                R=gauss(R,d-1);
            case 3
                R=@(s,t) P(s,0,t)'*P(s,0,t);
                R=gauss(R,d-1);
            case 4
                R=@(s,t) P(s,t,0)'*P(s,t,0);
                R=gauss(R,d-1);
        end
        %Nodos elemento actual
        w=(mesh.p(mesh.t(ind,:),1:d))';
        %Interpolacion g en elemento actual
        G=gn(w(1,:)',w(2,:)',w(3,:)');
        %Producto interior de traza local: 
        %\int_{\Gamma_N\cap T} g\cdot v
        F2=n*R*G;
        %Ensamble vector de carga
        Ind=[mesh.t(ind,:) , mesh.t(ind,:)+nnodos , mesh.t(ind,:)+2*nnodos];
        F(Ind)=F(Ind)+F2;
    end
end

%% Resolucion sistema

u(nodirich) = K(nodirich,nodirich)\(F(nodirich)-K(nodirich,Dirich)*u(Dirich));

%% Exportar VTK y visualizacion

%Parte vectorial---Defleccion
sol.type='vector';
sol.name='Deflection';
sol.data=[u(1:nnodos),u(nnodos+1:2*nnodos),u(2*nnodos+1:3*nnodos)];
%Parte escalar---Presion
sol(2).type='scalar';
sol(2).name='Pressure';
sol(2).data=u(3*nnodos+1:end);

name=strcat(name,'.vtk');
vtk_write_tetrahedral_grid_and_data(name,'Linear_elasticity',...
    mesh.p(:,1:end-1),mesh.t,sol,false)

%Visualizar en paraview
c=strcat('paraview --data=''',name,'''');
system(c);
