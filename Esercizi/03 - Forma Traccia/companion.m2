--struttura esercitazione 11/5/2018
R=QQ[x]
f=x^5+x+1
d=degree I
--costruiamo la matrice compagna
mm=mutableMatrix(R,d,d)
for i from 0 to d-1 do for j from 0 to d-1 do
mm_(i,j)=sub(contract(x^i,x*x^j%I),x=>0_QQ)
comp=matrix mm
comp^5+compx+1
---i due comandi seguenti danno la stessa matrice
for i from 0 to d-1 do for j from 0 to d-1 do mm_(i,j)=sub(contract(x^i,x*x^j%I),x=>0_QQ)
comp2=matrix mm
comp^2
comp^2-comp2---matrice nulla
---------------
eig=eigenvalues(sub(comp,CC))
eig_0
sum(5,i->eig_i)
sum(5,i->(eig_i)^2)
trace(comp),trace(comp^2)
comp^2
------------------------
----costruzione della forma traccia
b=mutableMatrix(R,d,d)
for i from 0 to d-1 do for j from 0 to d-1 do b_(i,j)=trace(comp^(i+j))
bez=matrix b
rank bez---siccome il rango è 5, ci sono 5 radici distinte 
eigenvalues(sub(bez,CC))
det(bez-x*id_(R^{d:0}))----dalle variazioni di bez si calcola il numero di autovalori positivi e quindi la segnatura di bez

restart
R=QQ[x]
----cerco un polinomio che per x=1  sia pol. 2grado assegnato
---per x=2 sia pol. 1grado assegnato
q=(quotientRemainder(matrix{{1_R}},matrix{{(x-2)^2,(x-1)^3}}))_0
b1=q_(0,0)
b2=q_(1,0)
g1=(x-2)^2
g2=(x-1)^3
b1*g1+b2*g2
h=((2+3*(x-1)+5*(x-1)^2)*b1*g1+(7+11*(x-2))*b2*g2)%ideal((x-2)^2*(x-1)^3)
sub(h,x=>1)
sub(h,x=>2)
----esercizi per casa

----trovare un polinomio che per x=i vale al primo ordine 2*i-(x-i)
---per i=1,2,3

---trovare quante radici reali ha il polinomio 2018*x^5+5*x^4+11*x^3+yyyy*x^2+mm*x+dd
---dove  dd/mm/yyyy è la propria data di nascita
