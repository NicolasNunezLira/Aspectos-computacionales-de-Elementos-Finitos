%Aprox de la solicion del problema de valores en la frontera
% -\epsilon u+a\cdot grad(u)+bu=f en \Omega\subset R^2
%  u_0 en \Gamma_D
%  grad(u)*n=0 en \gamma_N
clear all;close all;clc;
d=2;%Dimension del espacio
SUPG=0;%0: Elementos finitos clasico,1:SUPG
%% Creaci�n mallado
%[archivo,ruta]=uigetfile('*.poly','Archivo de malla');
%name = archivo(1:end-5);
%filename = strcat('triangle -pqa.008'," ", archivo);
%system(sprintf('%s',filename)) 

%fid=fopen(strcat(ruta,[name,'.1.node']),'r');
%nid = fscanf(fid,'%i',4);
%mesh.p = fscanf(fid,'%f',[nid(3)+4 nid(1)]);
%mesh.p = (mesh.p(2:end,:))';
%fclose(fid);

%fid=fopen(strcat(ruta,[name,'.1.ele']),'r');
%nele = fscanf(fid,'%i',3);
%mesh.t = fscanf(fid,'%f',[nele(3)+4 nele(1)]);
%mesh.t = (mesh.t(2:end,:))';
%fclose(fid);

%% Creacion malla ubuntu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Código para cambiar el mallado del dominio:
%"    !./triangle -pqea0.1 Test.poly     "
%triángulos con area a lo mas 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

file = fopen('hem.1.node');
num=fscanf(file ,'%i',4); nverf=num(1);
mesh.p=fscanf(file ,'%f',[4,nverf]);
mesh.p = mesh.p';
fclose(file);
mesh.p=mesh.p(:,2:4);

file = fopen('hem.1.ele');
num=fscanf(file,'%i',3); nelf=num(1);
mesh.t=fscanf(file,'%i',[4,nelf]);
mesh.t = mesh.t';
fclose(file);
mesh.t=mesh.t(:,2:4);


%% Visualizar malla
%trisurf(mesh.t(:,1:3),mesh.p(:,1),mesh.p(:,2),rand(size(mesh.p(:,1))))
%view(2)

%% Datos del problema
a = [ 1 ; 0 ]' ;
b= 0 ;
epsilon=10^(-4);
f=@(x,y) x-x;%Termino fuente
u_0= {@(x,y) x-x+y-y,@(x,y) x.^0.*y.^0};%Cond dirichlet

%% Variables a utilizar
P=@(x,y)[1-x-y,x,y];
DP=[-1 1 0;-1 0 1];
n1=size(mesh.t,1);
n2=size(mesh.p,1);
K=sparse(n2,n2);
F=sparse(n2,1);

%% Ensamble
for i=1:n1
 %% Identificacion del elemento i
 ind=mesh.t(i,:);
 xx=mesh.p(ind,1:2)';
 %% Cambio de variable
 B=xx(:,2:end)-xx(:,1);
 Jf=abs(det(B));
 bb=xx(:,1);
 %% Matriz de rigidez del elemento i
 aT= @(x,y) Jf*(epsilon*DP'*(B\(B'\DP)) + (P(x,y)'*a)*(B'\DP) + b*P(x,y)'*P(x,y));
 KT=gauss(aT,d);
 switch SUPG
     case 0
        K(ind,ind)=K(ind,ind)+KT;
     case 1
         chi=@(x) coth(x)-1./x;
         hk = 2*norm(a)/(sum(abs(DP'*(B\a'))));
         Pe = (norm(a)*hk)/(2*1*epsilon);
         dk = hk*chi(Pe)/(2*norm(a));
         r=@(x,y,w)(dk*(B\DP)'*a')*(a*(B\DP) + b*P(x,y) - f(w(1),w(2)))*Jf;
         R = gauss(r,d,B,bb);
         K(ind,ind)=K(ind,ind)+KT+R;
 end
 %% Ensamblaje de vector de traccion
 t = @(x,y) Jf*(P(x,y)'*P(x,y)*(f(xx(1,:),xx(2,:)))');
 FT=gauss(t,d);
 F(ind)=F(ind)+FT;
end
%% Imponer Condiciones Dirichlet
uu=zeros(n2,1);
id=[2,3];
nodirich=1:n2;
for i=1:length(id)
    dirich=find(mesh.p(:,3)==id(i));
    dirich = unique(dirich);
    uu(dirich)=u_0{i}(mesh.p(dirich,1),mesh.p(dirich,2));
    nodirich=setdiff(nodirich,dirich);
end
%% Calculo de los nodos restantes
uu(nodirich) = K(nodirich,nodirich)\(F(nodirich)-K(nodirich,setdiff(1:n2,nodirich))*uu(setdiff(1:n2,nodirich)));

%% Graficos
figure
trisurf(mesh.t(:,1:3), mesh.p(:,1),mesh.p(:,2),uu)
shading(gca,'interp');
colormap(jet);
colorbar;
title('Ecuacion de Hemker con MEF clasico')


exportTriangulation2VTK('Caso2D',[mesh.p(:,1:2) uu],mesh.t);
