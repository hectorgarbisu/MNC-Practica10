
function [x,y] = miMDF(a,b,N,ya,yb,func,params)
% Metodo de Diferencias Finitas para problemas de contornos
% ULPGC, EII, MNC

    % espaciado constante
    x = linspace(a,b,N);
    h = (b-a)/(N-1);
    
    % función de cálculo de los coeficientes
    P = func(x,params);
    A = P(1,:);
    B = P(2,:);
    C = P(3,:);
    D = P(4,:);
    
    % factores E,F, G y H
    E = A;
    F = h*h*C - h*B - 2*A;
    G = h*B + A;
    H = h*h*D;
     
    % Ensamblado de RHS
    RHS = H';
    RHS(1)= ya;
    RHS(N) = yb;
    
    % Ensamblado de LHS
    LHS = zeros(N,N);
    LHS(1,1)=1;
    LHS(N,N) = 1;
    for a=2:N-1
        LHS(a,a-1)= E(a);
        LHS(a,a) = F(a);
        LHS(a,a+1) = G(a);
    end
    
    %disp(LHS);
    %disp(RHS);
    
    % si el determinante es nulo o muy pequeño,
    % el sistema es no compatible -> no hay solución
    deter = det(LHS);
    if abs(deter) < 1.0e-6 
        error('El sistema no tiene solución');
    end
    X = LHS\RHS;
    y = X';
end
