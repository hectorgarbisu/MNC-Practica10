function [x,y] = inicioMDF(a,b,N,ya,y1a,func,params)
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
    RHS = zeros(N,1);
    RHS(1)   = ya;
    RHS(2)   = h*y1a;
    RHS(2:N) = H(1:N-1)';
    
    % Ensamblado de LHS
    LHS = zeros(N,N);
    LHS(1,1) = 1;
    LHS(2,1) = -1;
    LHS(2,2) = 1;
    for a=3:N
        LHS(a,a-2)= E(a-2);
        LHS(a,a-1) = F(a-2);
        LHS(a,a) = G(a-2);
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
