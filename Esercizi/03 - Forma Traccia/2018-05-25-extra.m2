----esercizio in due variabili
restart
R=QQ[x,y]
I=ideal((x+y)^4+2*x*y^2+7,x^2+y^2-x-3)
bb=sub(basis (R/I),R)
d=numcols basis(R/I)
--prima colonna della matrice compagna rispetto a x
compx=sub(contract(transpose bb,(bb_(0,0))*x%I),{x=>0_QQ,y=>0_QQ})
--aggiungo le colonne e trovo compx
for i from 1 to (numcols bb-1) do 
compx=compx|sub(contract(transpose bb,(bb_(0,i))*x%I),{x=>0_QQ,y=>0_QQ})
compx
compy=sub(contract(transpose bb,(bb_(0,0))*y%I),{x=>0_QQ,y=>0_QQ})
for i from 1 to (numcols bb-1) do 
compy=compy|sub(contract(transpose bb,(bb_(0,i))*y%I),{x=>0_QQ,y=>0_QQ})
----esercizio, valutare numericamente le 8 soluzioni del sistema
compy
bb
sols=(eigenvectors(transpose (compx)))_1
for i from 0 to 7 do print(i,(sols_i)_1/(sols_i)_0,(sols_i)_5/(sols_i)_0)
---verifica in un caso che le soluzioni annullano numericamente il sistema
i=3
a=(sols_i)_1/(sols_i)_0,b=(sols_i)_5/(sols_i)_0
sub(gens I,{x=>a,y=>b})

----bezoutiante
b=mutableMatrix(R,d,d)
for i from 0 to d-1 do for j from 0 to d-1 do
b_(i,j)=trace(compx^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,0))*compy^((matrix exponents(bb_(0,i)*bb_(0,j)))_(0,1)))
bez=matrix b
rank bez, det bez
det(bez-x*id_(R^{8:0}))---4 variazioni, ci sono 0 radici reali per I e 4 coppie complesse coniugate
eigenvalues sub(bez,RR)---sono 4 e 4
