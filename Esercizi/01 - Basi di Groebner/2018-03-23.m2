----- Esercizio 0.1.5 -----
R = QQ[x]
f = x^16-4*x^12+6*x^8-4*x^4+1
g = sum(10, i -> (-1)^i * x^i)
unbox = (m) -> m_(0,0)
-- Usando gcd
d1 = gcd(f,g)
-- Usando la principalità degli ideali
d2 = unbox(gens(gb(ideal(f,g))))
-- Intersecando gli ideali
m3 = unbox(gens(intersect(ideal(f),ideal(g))))
d3 = sub(f*g/m3, R)
-- Eliminando t a mano
RT = QQ[t,x,MonomialOrder=>Lex]
I = ideal(t*sub(f,RT), (1-t)*sub(g,RT))
m4 = unbox(gens(gb(I)))
d4 = sub((sub(f,RT)*sub(g,RT))/m4, R)
-- Eliminando t in modo automatico
m5 = unbox(gens(eliminate({t},I)))
d5 = sub((sub(f,RT)*sub(g,RT))/m5, R)

----- Esercizio 0.1.6 -----
R = QQ[x,y]
I = ideal(x^4+y^4-1, x^2-x*y+y^2-1)
f = x^8+y^8-1
f%I

----- Esercizio 0.2.1 -----
R = QQ[a..f, MonomialOrder=>Lex]
I = ideal(sum(6,i->(R_i)^3), sum(6,i->(R_i)^2))
E1 = eliminate({a},I)
E2 = ideal((gens(gb(I)))_(0,0))
E1 == E2

----- Esercizio 0.2.2 -----
R = QQ[t,x,y,z,MonomialOrder=>Lex]
f = x^3+2*x*y*z-z^2
g = x^2+y^2+z^2-1
xyz = matrix{{x,y,z}}
I = ideal(diff(xyz,f)-t*diff(xyz,g),g)
G = gens(gb(I))
factor(G_(0,0))
-- Trovate le radici del polinomio nella sola z,
-- possiamo sostituirle negli altri polinomi di G
-- e determinare le x e y dei punti critici
sub(G,z=>0)
sub(G,z=>1)
sub(G,z=>-1)
xyz%sub(G,z=>2/3)
xyz%sub(G,z=>-2/3)
sub(G,R/ideal(128*z^2-11))

----- Esercizio 0.2.3 -----
R = QQ[x,y]
I = ideal(x^2-x+1, y-1)
J = ideal(x^3*y+x*y^3-2, x+y-1)
K1 = intersect(I,J)

RT = QQ[t,x,y,MonomialOrder=>Lex]
G2 = gens(gb(t*sub(I,RT)+(1-t)*sub(J,RT)))
K2 = sub(ideal(submatrix(G2, {0}, 0..2)), R)

K1 == K2 --true

----- Esercizio 0.2.4 -----
R = QQ[x,y,z]
f = x^3+y^3+z^3-3*x*y*z
g = x+y+z
d = gcd(f,g)
m = lcm(f,g)
