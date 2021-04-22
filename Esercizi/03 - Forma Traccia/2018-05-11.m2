----- Polinomio Minimo di una matrice -----
restart
n = 3
R = QQ[a_0..a_n,t]
m = sub(matrix{
    {1,1,0},
    {0,1,0},
    {0,0,3}}, R)
ide = m^0

polcar = det(m - t*ide)
factor(polcar)

im = minors(1,a_0*m^0+a_1*m^1+a_2*m^2+a_3*m^3)
coeff = matrix{{a_0..a_n}}%im
factor((coeff*matrix{{1},{t},{t^2},{t^3}})_(0,0))

load("minpoly.m2")

----- Matrice Compagna -----
restart
R = QQ[x]
f = x^5+x+1
I = ideal(f)
mm = mutableMatrix(R,5,5)
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = sub(contract(x^i,x^j*x%I),x=>0)
comp = matrix(mm)
eigenvalues(sub(comp,CC))
eigenvectors(sub(comp,CC))

h = x+5*x^3-x^7
for i from 0 to 4 do
    for j from 0 to 4 do
    	mm_(i,j) = sub(contract(x^i,x^j*h%I),x=>0)
comp2 = matrix(mm)
comp2 == (comp + 5*comp^3 - comp^7)

comp^0 + comp + comp^5 == 0

----- Forma Traccia di Killing
restart
R = QQ[x]
f = (x^2+1)*(x-1)*x*(x+1)
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
det(bez)
det(bez-bez^0*x)
-- 4 variazioni di segno su 5, quindi per la regola di Cartesio
-- bez ha 4 autovalori positivi e 1 negativo (il prossimo
-- comando lo verifica numericamente).
-- 4-1 = 3, infatti per costruzione f ha 3 radici reali.
eigenvalues(sub(bez,CC))
for i from 0 to 4 do
    print(i, det(submatrix(bez,{0..i},{0..i})))

----- Interpolazione di Hermite -----
restart
R = QQ[x]
a1 = 5+(x-1)-13*(x-1)^2
a2 = 2018+5*(x-2)
f = (x-1)^3 * (x-2)^2
g1 = (x-2)^2
g2 = (x-1)^3
q = (quotientRemainder(matrix{{1_R}},matrix{{g1,g2}}))_0
b1 = q_(0,0)
b2 = q_(1,0)
polint = (g1*b1*a1+g2*b2*a2)%f
sub(polint, x=>1)
sub(polint, x=>2)

-- PER CASA: Trovare un polinomio che per x=i vale al primo ordine
-- 2*i-(x-i) per i=1,2,3. Calcolare il numero di radici reali di
-- 2018*x^5+5*x^4+11*x^3+yy*x^2+mm*x+dd.