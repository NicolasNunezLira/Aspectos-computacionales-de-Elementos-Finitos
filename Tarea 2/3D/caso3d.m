%Aprox de la solicion del problema de valores en la frontera
% -\epsilon u+a\cdot grad(u)+bu=f en \Omega\subset R^3
%  u_0 en \Gamma_D
%  grad(u)*n=0 en \gamma_N
clear all;close all;clc;
d=3;%Dimension
SUPG=1;%0: Elementos finitos clasico,1:SUPG
%% Creación mallado
[archivo,ruta]=uigetfile('*.poly','Archivo de malla');
name = archivo(1:end-5);
filename = strcat('tetgen -pqea.008k'," ", archivo);
system(sprintf('%s',filename))

fid=fopen(strcat(ruta,[name,'.1.node']),'r');
nid = fscanf(fid,'%i',4);
mesh.p = fscanf(fid,'%f',[nid(3)+(nid(4)+1+d) nid(1)]);
mesh.p = (mesh.p(2:end,:))';
fclose(fid);

fid=fopen(strcat(ruta,[name,'.1.ele']),'r');
nele = fscanf(fid,'%i',3);
mesh.t = fscanf(fid,'%f',[nele(3)+(2+d) nele(1)]);
mesh.t = (mesh.t(2:end,:))';
fclose(fid);

fid=fopen(strcat(ruta,[name,'.1.face']),'r');
nface = fscanf(fid,'%i',2);
mesh.f = fscanf(fid,'%f',[1+d+nface(2) nface(1)]);
mesh.f = (mesh.f(2:end,:))';
fclose(fid);


%% Visualizar malla
%figure
%trisurf(mesh.t(:,1:3),mesh.p(:,1),mesh.p(:,2),rand(size(mesh.p(:,1))))
%view(2)

%% Datos del problema
a = [ 1 ; 0 ;0]' ;
b= 0 ;
epsilon=10^(-4);
f=@(x,y,z) x-x;%Termino fuente
u_0= {@(x,y,z) x-x,@(x,y,z) x.^0.*y.^0.*z.^0};%Cond dirichlet

%% Variables a utilizar
P=@(x,y,z) [1-x-y-z;x;y;z]';
DP=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
n1=size(mesh.t,1);
n2=size(mesh.p,1);
K=sparse(n2,n2);
F=sparse(n2,1);

%% Ensamble
 
for i=1:n1
 %% Identificacion del elemento i
 ind=mesh.t(i,:);
 xx=mesh.p(ind,1:3)';
 %% Cambio de variable
 B=xx(:,2:end)-xx(:,1);
 Jf=abs(det(B));
 bb=xx(:,1);
 %% Matriz de rigidez del elemento i
 aT= @(x,y,z) Jf*(epsilon*DP'*(B\(B'\DP)) + (P(x,y,z)'*a)*(B'\DP) + b*P(x,y,z)'*P(x,y,z));
 KT=gauss(aT,3);
  switch SUPG
     case 0
        K(ind,ind)=K(ind,ind)+KT;
     case 1
         chi=@(x) coth(x)-1./x;
         hk = 2*norm(a)/(sum(abs(DP'*(B\a'))));
         Pe = (norm(a)*hk)/(2*epsilon);
         delta = hk*chi(Pe)/(2*norm(a));
         r=@(x,y,z,w)(delta*(B\DP)'*a')*(a*(B\DP) + b*P(x,y,z) - f(w(1),w(2),w(3)))*Jf;
         R = gauss(r,d,B,b);
         K(ind,ind)=K(ind,ind)+KT+R;
 end
 %% Ensamblaje de vector de traccion
 t = @(x,y,z) Jf*(P(x,y,z)'*P(x,y,z)*(f(xx(1,:),xx(2,:),xx(3,:)))');
 FT=gauss(t,3);
 F(ind)=F(ind)+FT;
end
%% Imponer Condiciones Dirichlet
uu=zeros(n2,1);
id=[1,3];
nodirich=1:n2;
for i=1:length(id)
    dirich=find(mesh.p(:,4)==id(i));
    dirich = unique(dirich);
    uu(dirich)=u_0{i}(mesh.p(dirich,1),mesh.p(dirich,2),mesh.p(dirich,3));
    nodirich=setdiff(nodirich,dirich);
end
%% Calculo de los nodos restantes
uu(nodirich) = K(nodirich,nodirich)\(F(nodirich)-K(nodirich,setdiff(1:n2,nodirich))*uu(setdiff(1:n2,nodirich)));

%% Graficos
figure
trisurf(mesh.t(:,1:4), mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),uu)
title('Ecuacion de Hemker con MEF SUPG')
%shading(gca,'interp');
colormap(jet);
colorbar;

writeVTK('Caso3D',mesh.t(:,1:4),mesh.p(:,1:end-1),uu)
