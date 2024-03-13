function [L2,H1,sH1]=err(f,df,p,t,uu,d)
%Calculo de los errores en normas L^2 y H^1 y seminorma H^1 dados los datos de ingreso:
%   f:Funcion solucion
%   df:Gradiente de la solucion
%   p:Nodos del mallado generado por Triangle/Tengen
%   t:Elementos del mallado generado por Triangle/Tengen
%   uu:Solucion numerica
%   d:Dimension del espacio(d=2,3)
L2=0;sH1=0;
for i=1:size(t,1)
    switch d
        case 2
            P=@(x,y)[1-x-y,x,y];
            DP=[-1 1 0;-1 0 1];
            ind=t(i,:);
            xx=p(ind,1:2)';
            B=xx(:,2:end)-xx(:,1);
            Jf=abs(det(B));
            b=xx(:,1);
            l=@(x,y,w) (f(w(1),w(2))-P(x,y)*uu(ind)).^2*Jf;
            L2=L2+gauss(l,d,B,b);
            sn=@(x,y,w) (df(w(1),w(2))-(B'\DP)*u(ind)).^2*Jf;
            sH1=sH1+gauss(sn,d,B,b);
        case 3
            P=@(x,y,z) [1-x-y-z;x;y;z]';
            DP=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
            ind=t(i,:);
            xx=p(ind,1:3)';
            B=xx(:,2:end)-xx(:,1);
            Jf=abs(det(B));
            b=xx(:,1);
            l=@(x,y,w) (f(w(1),w(2),w(3))-P(x,y,z)*uu(ind)).^2*Jf;
            L2=L2+gauss(l,d,B,b);
            sn=@(x,y,w) (df(w(1),w(2),w(3))-(B'\DP)*u(ind)).^2*Jf;
            sH1=sH1+gauss(sn,d,B,b);
    end
end
H1=sqrt(L1+sH1);
sH1=sqrt(sH1);
L2=sqrt(L2);
