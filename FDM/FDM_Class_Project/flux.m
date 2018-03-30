function F = flux(u,v)
    if u>=v && (u+v)/2 > 0
        F = .5*u^2;
    elseif u>=v && (u+v)/2<=0;
        F = .5*v^2;
    elseif u<v && u>0;
        F = .5*u^2;
    elseif u<v && v<0;
        F = .5*v^2;
    else u<v && (u<= 0 <=v);
        F = 0;
    end
end