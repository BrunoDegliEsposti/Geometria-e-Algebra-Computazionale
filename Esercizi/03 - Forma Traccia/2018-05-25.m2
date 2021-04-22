----- Esercitazione
-- f = x^5+15*x^4+7*x^3+a*x^2-x+3
-- Per quali valori di a f ha radici doppie?
-- Come varia il numero di radici reali di f per a = 1,3,5?
restart
R = QQ[x,a]
f = x^5+15*x^4+7*x^3+a*x^2-x+3
df = diff(x,f)
gcd(f,df)
d = discriminant(f,x)

R1 = QQ[a]
d = sub(d,R1)
I = ideal(d)
mm = mutableMatrix(R,5,5)
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = sub(contract(a^i,a^j*a%I),a=>0)
comp = matrix(mm)
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = trace(comp^(i+j))
bez = matrix(mm)
eigenvalues(sub(comp,CC)) -- tre radici reali distinte
det(bez) -- essendo negativo, gli autovalori negativi di bez
-- sono dispari. non possono essere 3 o 5 altrimenti
-- ci sarebbero troppi autovalori in totale (dato che
-- ogni coppia di radici complesse coniugate contribuisce
-- un autovalore positivo e uno negativo).
-- Quindi abbiamo conferma delle tre radici reali distinte.

-- Radici positive e negative?
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = trace(comp^(i+j+1))
beza = matrix(mm)
det(beza) -- essendo positivo, ci sono o 0 o 2 radici positive
eigenvalues(sub(beza,CC),Hermitian=>true) -- 0 radici positive

-- Radici dentro (-500,-400) o fuori?
id5 = id_(R1^5)
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = trace(-(comp+500*id5)*(comp+400*id5)*comp^(i+j))
bezh = matrix(mm)
det(bezh)
eigenvalues(sub(bezh,CC),Hermitian=>true)

-- Osserviamo che gli unici cambiamenti del numero di
-- radici reali di f avviene in corrispondenza di radici doppie
-- Quindi basta calcolarne il numero nei quattro
-- intervalli delimitati dalle tre radici

R = QQ[x,a]
f = x^5+15*x^4+7*x^3+a*x^2-x+3
I = ideal(f)
mm = mutableMatrix(R,5,5)
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = sub(contract(x^i,x^j*x%I),x=>0)
comp = matrix(mm)
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = trace(comp^(i+j))
bez = matrix(mm)
p = det(bez-bez^0*x)
sub(p,a=>-500) -- 4 variazioni di segno -> 3 radici reali
sub(p,a=>-400) -- 5 variazioni di segno -> 5 radici reali
sub(p,a=>-15)  -- 4 variazioni di segno -> 3 radici reali
sub(p,a=>-5)   -- 3 variazioni di segno -> 1 radice reale


----- Esercitazione 2 -----
restart
R=QQ[x,y]
I=ideal((x+y)^4+2*x*y^2,x^2+y^2-x-3) --nota il -3
bb=sub(basis (R/I),R)
d=numcols basis(R/I)

--prima colonna della matrice compagna rispetto a x
compx=sub(contract(transpose bb,(bb_(0,0))*x%I),{x=>0_QQ,y=>0_QQ})
--aggiungo le colonne e trovo compx
for i from 1 to (numcols bb-1) do 
compx=compx|sub(contract(transpose bb,(bb_(0,i))*x%I),{x=>0_QQ,y=>0_QQ})
compx
det(compx-x*id_(R^{8:0}))
eigenvalues(compx)

compy=sub(contract(transpose bb,(bb_(0,0))*y%I),{x=>0_QQ,y=>0_QQ})
for i from 1 to (numcols bb-1) do 
compy=compy|sub(contract(transpose bb,(bb_(0,i))*y%I),{x=>0_QQ,y=>0_QQ})
compy
det(compy-y*id_(R^{8:0}))
eigenvalues(compy)

b=mutableMatrix(R,d,d)
for i from 0 to d-1 do for j from 0 to d-1 do
b_(i,j)=trace(compx^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,0))*compy^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,1)))
bez=matrix b
rank(bez), det(bez)
det(bez-x*id_(R^8))
eigenvalues(sub(bez,CC))
