function int = gauss(f,d,B,b)
%Aproximacion de integrales con cuadratura de Gauss en 2D y 3D, en triangulos y tetraedros respectivamente,
%considerando un cambio de variable lineal, afin e invertible de la forma
%  F:\hat{K}\rightarrow K
%  F(\hat{x})=B*\hat{x}=x
%  donde \hat{K} es el triangulo o el tetraedro de referencia con nodos \theta,\e_i,i=1:d
% Parametros de entrada:
% f:Integrando
% d:Dimension del espacio
% B,b:Constantes de F
%----------------------
%Si nargin=2, se considera que f solo se debe evaluar \hat{K}.
%Si nargin=4,se considera que f depende de vectores en \hat{K} y K
if nargin==2
switch d
    case 2
        ww = [1/6, 1/6, 1/6];
        xx = {2/3,1/6;1/6,1/6;1/6,2/3};
        int = 0;
        if nargin==2
            for j = 1:d+1
                int = int + ww(j)*f(xx{j,:});
            end
        end
        if nargin==4
            for j = 1:d+1
                int = int + ww(j)*f(xx{j,:},[xx{j,:}]*B'+b);
            end
        end
    case 3
        ww = [1/24,1/24,1/24,1/24];
        xx = {0.5854101966249685  0.1381966011250105  0.1381966011250105;
            0.1381966011250105  0.1381966011250105  0.1381966011250105;
            0.1381966011250105  0.1381966011250105  0.5854101966249685;
            0.1381966011250105  0.5854101966249685  0.1381966011250105};
        int = 0;
        if nargin==2
            for j = 1:d+1
                int = int + ww(j)*f(xx{j,:});
            end
        end
        if nargin==4
            for j = 1:d+1
                int = int + ww(j)*f(xx{j,:},[xx{j,:}]*B'+b);
            end
        end
    end
end