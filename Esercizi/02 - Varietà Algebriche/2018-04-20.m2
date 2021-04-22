----- Esercizi risultante -----
R = QQ[x,y]
f = random(4,R)
g = random(5,R)
sylv = sylvesterMatrix(f,g,x)
det(sylv)
resultant(f,g,x)
discriminant(f,x)
-----
restart
R = QQ[a_0..a_5,x]
for n from 2 to 5 do (
    f = sum(n+1, i->a_i*x^i),
    print(f, factor(discriminant(f,x)))
)

----- Equazione implicita di una
----- superficie razionale.
----- Procedimento errato:
restart
R = QQ[u,v,x,y,z]
I = ideal(x*v-u^2, y*u-v^2, z-u)
eliminate({u,v},I)
----- Procedimento corretto:
restart
R = QQ[t,u,v,x,y,z]
I = ideal(x*v-u^2, y*u-v^2, z-u, 1-t*u*v)
eliminate({t,u,v},I)

----- Folium di Cartesio -----
restart
R = QQ[t,x,y]
I = ideal((1+t^3)*x-3*t, (1+t^3)*y-3*t^2)
eliminate(t,I)
