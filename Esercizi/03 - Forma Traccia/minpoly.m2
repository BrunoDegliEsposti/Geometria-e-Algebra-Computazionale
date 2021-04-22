--INPUT: matrice A nxn
--OUTPUT: polinomio minimo di A
n=3
R=QQ[a_0..a_n,t]
--il comando seguente definisce una matrice A random con coefficienti razionali,
--puo' essere sostituito da una matrice A a piacere con la forma A=matrix{{...}}
A=random(R^{n:0},R^{n:0})
A=matrix{{0_R,0,-1},{1,0,-3},{0,1,-3}}
--il ciclo seguente mostra se le matrici A^0..A^(k-1) sono indipendenti calcolando il nucleo
for k from 1 to (n+1) do 
print(k,gens kernel transpose diff(transpose matrix {{a_0..a_(k-1)}},gens minors(1,sum(k,i->a_i*A^i))))
--il ciclo seguente automatizza il calcolo del primo indice k tale che A^0..A^(k-1) sono dipendenti
nk=0,k=0,while nk==0 do (k=k+1,nk=numcols gens kernel transpose diff(transpose matrix {{a_0..a_(k-1)}},gens minors(1,sum(k,i->a_i*A^i))))
--il seguente p e' il polinomio minimo
p=sum(k,i->t^i*(gens kernel transpose diff(transpose matrix {{a_0..a_(k-1)}},gens minors(1,sum(k,i->a_i*A^i))))_(i,0))
if gcd(p,diff(t,p))==1 then print "A è diagonalizzabile su C" else print "A non è diagonalizzabile su C"

--matrice compagna