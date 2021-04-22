-- Esercizio 2
R = QQ[x,y,z,a..f]
M = matrix({{a,b,d},
	    {b,c,e},
	    {d,e,f}})
I1 = minors(2,M)
I1 = radical(I1)
isHomogeneous(I1)
parametrizzazione = (x,y,z) -> {x^2, x*y, y^2, x*z, y*z, z^2}
I2 = eliminate({x,y,z}, ideal(toList(a..f) - parametrizzazione(x,y,z)))
I1 == I2

restart

--Esercizio 3
R = QQ[x,y,a]    -- è fondamentale dare a x e y la precedenza su a
f = (x+y)^2+2*x*y-y
g = x^2+y^2-a*x+1
I = ideal(f,g)
leadTerm(I)      -- I è zero-dimensionale in x,y. Bene che a sia assente.
base = matrix{{1,x,y,y^2}}    -- Q-base di R/I per ogni a fissato
d = 4

compx = sub(contract(transpose(base), base*x%I), {x=>0_QQ,y=>0_QQ})
compy = sub(contract(transpose(base), base*y%I), {x=>0_QQ,y=>0_QQ})
m = (transpose base) * base
b = mutableMatrix(R,d,d)
for i from 0 to d-1 do (
    for j from 0 to d-1 do (
        expx = (exponents(m_(i,j)))#0#0;
	expy = (exponents(m_(i,j)))#0#1;
        b_(i,j) = trace(compx^expx * compy^expy)
    )
)
bez = matrix(b)

det(sub(bez,a=>-1)-x*id_(R^d)) -- 2-2 = 0 punti reali
det(sub(bez,a=>0)-x*id_(R^d))  -- 2-2 = 0 punti reali
det(sub(bez,a=>3)-x*id_(R^d))  -- 3-1 = 2 punti reali

h = det(bez)
J = ideal(x,y,h)
baseh = sub(basis(R/J), R)
comph = sub(contract(transpose(baseh), baseh*a%J), a=>0_QQ)
eigenvalues(sub(comph,CC))
b = mutableMatrix(R,6,6)
for i from 0 to 5 do
    for j from 0 to 5 do
	b_(i,j) = trace(comph^(i+j))
bezh = matrix(b)
det(bezh - x*id_(R^6))         -- 5-1 = 4 radici reali di h

det(sub(bez,a=>29)-x*id_(R^d))
