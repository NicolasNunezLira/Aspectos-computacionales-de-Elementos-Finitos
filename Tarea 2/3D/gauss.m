function int = gauss(f,d,B,b)
%Calculo de integrales en dimension d=2,3 por cuadratura de
%Gauss en el elemento de referencia \hat{T} correspondiente para los datos
%   f:Funcion anonima a integrar
%   d:Dimension del espacio
%   B,b: Parametros de la aplicacion afin F(\hat{T})=T
if nargin==4
switch d
    case 2
        
        ww = [1/6, 1/6, 1/6];
        xx = {2/3,1/6;1/6,1/6;1/6,2/3};
        int = 0;
        for j = 1:d+1
            int = int + ww(j)*f(xx{j,:},[xx{j,:}]*B'+b);
        end
    case 3
        ww = [1/24,1/24,1/24,1/24];
        xx = {0.5854101966249685  0.1381966011250105  0.1381966011250105;
            0.1381966011250105  0.1381966011250105  0.1381966011250105;
            0.1381966011250105  0.1381966011250105  0.5854101966249685;
            0.1381966011250105  0.5854101966249685  0.1381966011250105};
        int = 0;
        for j = 1:d+1
            int = int + ww(j)*f(xx{j,:},[xx{j,:}]*B'+b);
        end
end
elseif nargin==2
switch d
    case 2
        
        ww = [1/6, 1/6, 1/6];
        xx = {2/3,1/6;1/6,1/6;1/6,2/3};
        int = 0;
        for j = 1:d+1
            int = int + ww(j)*f(xx{j,:});
        end
    case 3
        ww = [1/24,1/24,1/24,1/24];
        xx = {0.5854101966249685  0.1381966011250105  0.1381966011250105;
            0.1381966011250105  0.1381966011250105  0.1381966011250105;
            0.1381966011250105  0.1381966011250105  0.5854101966249685;
            0.1381966011250105  0.5854101966249685  0.1381966011250105};
        int = 0;
        for j = 1:d+1
            int = int + ww(j)*f(xx{j,:});
        end
    end
end
    
