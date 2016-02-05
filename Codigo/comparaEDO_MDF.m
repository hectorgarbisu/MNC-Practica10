function comparaEDO_MDF
    T=5;
    %muchos puntos para dibujo preciso de la solucion exacta
    tplot = linspace(0,T,200);
    yplot = SolucionExacta(tplot);
    
    figure;
    hold on;
    plot(tplot,yplot,'r-');
    
    %Solucion con ODE45()
    [tode, Yode] = ode45(@sistemaODE,[0,T],[0,0]);
    yode = Yode(:,1);
    plot(tode,yode,'b-');
    y1 = SolucionExacta(tode);
    %solucion exacta en los puntos de ODE45
    %se evitan los primeros terminos dado que on casi nulos
    %al ser error relativo implica dividir por casi cero
    errorODE = 100*mean(abs(y1(5:end)-yode(5:end))./y1(5:end));
    
    %solucion con MDF con 6 veces el numero de puntos
    N=length(tode);
    [tMDF,yMDF] = inicioMDF(0,T,6*N,0,0,@terminosMDF,[]);
    plot(tMDF,yMDF,'k-');
    %solucion exacta para los puntos de MDF
    y2 = SolucionExacta(tMDF);
    %Se evitan 6 veces los primeros terminos
    errorMDF = 100*mean(abs(y2(5*6:end)-yMDF(5*6:end))./y2(5*6:end));
    
    grid on;
    xlabel('t');
    title(sprintf('Error ODE(%%): %g, Error MDF(%%): %g',errorODE,errorMDF));
    hold off;
end

%solucion analitica de la EDO
function y = SolucionExacta(t)
    y = 1 -exp(-t/2).*(cos((sqrt(3)/2)*t)+(1/sqrt(3))*sin((sqrt(3)/2)*t));
end


function dy = sistemaODE(t,y)
    dy = zeros(2,1);
    dy(1,1) = y(2);
    dy(2,1) = 1-y(1)-y(2);
end

function P = terminosMDF(x,params)
    P = ones(4,length(x));
end
    