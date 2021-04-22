restart
R=QQ[x,y]
I=ideal((x+y)^4+2*x*y^2,x^2+y^2-x)
bb=sub(basis (R/I),R)
d=numcols basis(R/I)

--prima colonna della matrice compagna rispetto a x
compx=sub(contract(transpose bb,(bb_(0,0))*x%I),{x=>0_QQ,y=>0_QQ})
--aggiungo le colonne e trovo compx
for i from 1 to (numcols bb-1) do 
compx=compx|sub(contract(transpose bb,(bb_(0,i))*x%I),{x=>0_QQ,y=>0_QQ})
det(compx-x*id_(R^{8:0}))
eigenvalues(compx)

compy=sub(contract(transpose bb,(bb_(0,0))*y%I),{x=>0_QQ,y=>0_QQ})
for i from 1 to (numcols bb-1) do 
compy=compy|sub(contract(transpose bb,(bb_(0,i))*y%I),{x=>0_QQ,y=>0_QQ})
compy
det(compy-y*id_(R^{8:0}))
eigenvalues(compy)

eigenvalues(compx+compy) --si possono dedurre gli accoppiamenti
-- cercando due valori in compx e compy che sommati danno quelli
-- in questa lista 
sols=(eigenvectors(transpose (compx+compy)))_1

b=mutableMatrix(R,d,d)
for i from 0 to d-1 do for j from 0 to d-1 do
b_(i,j)=trace(compx^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,0))*compy^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,1)))
bez=matrix b
rank bez