restart
R = QQ[x,y]
I = ideal((x+y)^4+2*x*y^2,x^2+y^2-x)
bb = sub(basis(R/I),R)
d = numcols(bb) --numero di elementi della base (#colonne della matrice)
--prima colonna della matrice compagna rispetto a x
compx=sub(contract(transpose(bb),(bb_(0,0))*x%I),{x=>0_QQ,y=>0_QQ}) ---contrarre significa calcolare lo sviluppo di Taylor
---compx sta per matrice compagna di x
--aggiungo le colonne e trovo compx
for i from 1 to d-1 do
    compx=compx|sub(contract(transpose bb,(bb_(0,i))*x%I),{x=>0_QQ,y=>0_QQ})

compy=sub(contract(transpose bb,(bb_(0,0))*y%I),{x=>0_QQ,y=>0_QQ})
for i from 1 to d-1 do 
    compy=compy|sub(contract(transpose bb,(bb_(0,i))*y%I),{x=>0_QQ,y=>0_QQ})
det(compx-x*id_(R^8))
det(compy-x*id_(R^8))
eigenvalues(compx)
eigenvalues(compy)
eigenvalues(compx+compy)
---valutare numericamente le otto soluzioni del sistema sapendo gia che quattro coincidono con l`origine
compy
bb
sols=(eigenvectors(transpose (compx+compy)))_1

b=mutableMatrix(R,d,d)
for i from 0 to d-1 do for j from 0 to d-1 do
b_(i,j)=trace(compx^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,0))*compy^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,1)))
---il comando precedente calcola gli esponenti degli elementi della base
---in generale exponents prende in ingresso un monomio e restituisce il multiindice che gli fa da esponente
bez=matrix b
rank bez

contract(matrix{{2},{x},{3*x*x},{4*x*x*x}},x^2)
diff(1/2*x, x^2)
diff(matrix{{x},{y}},x*y)
