R = QQ[x_0..x_3];	--Di default M2 usa l'ordine lessicografico graduato inverso
I = ideal(sum(4,i->x_i^2),sum(4,i->x_i^3),sum(4,i->x_i^4));
time gens gb I;		--Quanto tempo ci mette a calcolare la base di groebner?
G = gens gb I; 		--Il ; sopprime l'output
betti G			--Quanti polinomi abbiamo per ciascun grado?

R = QQ[x_0..x_3, MonomialOrder=>Lex]
I = ideal(sum(4,i->x_i^2),sum(4,i->x_i^3),sum(4,i->x_i^4))
time G = gens gb I; {* Stavolta la base è più complicata e il suo calcolo più lento *}
leadTerm(I)
f=sum(4,i->x_i^5)
f%I

-- disomogenizzo l'ideale aggiungendo delle costanti
R = QQ[x_0..x_3, MonomialOrder=>Lex]
I = ideal(sum(4,i->x_i^2)+1,sum(4,i->x_i^3)+1,sum(4,i->x_i^4)+1)
time G = gens gb I;
betti G
leadTerm(I)
-- Morale: siamo in presenza di una forte instabilità algebrica

----- Esercizio 0.1.1 -----
R = QQ[x,y,MonomialOrder=>Lex];
I = ideal(x^2+y^2,x*y);
gg = gens gb I;
quotientRemainder(gg, gens I);

R = QQ[x,y,MonomialOrder=>GLex];
I = ideal(x^2+y^2,x*y);
gg = gens gb I;
quotientRemainder(gg, gens I);

R = QQ[x,y,MonomialOrder=>GRevLex];
I = ideal(x^2+y^2,x*y);
gg = gens gb I;
quotientRemainder(gg, gens I);

----- Esercizio 0.1.2 -----
-- Sintassi per le matrici in M2: matrix{{1,2},{3,4}}
n = 123456789;
m = 987654321;
d = gcd(n,m);
quotientRemainder(matrix{{d}},matrix{{n,m}})

----- Esercizio 0.1.3 -----
R = QQ[x_1..x_4]
A = matrix{{1,0,1,0},{0,1,0,-1},{1,2,3,4}}
I = minors(1,A*transpose basis(1,R))
gens gb I
basis(1,R)%I

----- Esercizio 0.1.4 -----
R = QQ[a_0..a_4,x]
aa = i->sum(5,j->(a_j)*i^j)
I = ideal(aa(1))
for i from 2 to 5 do I = I +ideal(aa(i)-((-1)^i*i^2+1))
--apply(5,i->(-1)^(i+1)*(i+1)^2+1)
matrix{{a_0..a_4}}%I
--A = matrix{aa(1),aa(2),aa(3),aa(4),aa(5)}





