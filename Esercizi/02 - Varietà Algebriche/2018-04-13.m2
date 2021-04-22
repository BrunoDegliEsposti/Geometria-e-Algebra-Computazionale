----- Esercizio in classe -----
R = QQ[x,y,MonomialOrder=>Lex];
I = ideal(x^2+2*y^2-2,x^2+x*y+y^2-2);
gens(gb(I))
R = QQ[y,x,MonomialOrder=>Lex];
I = ideal(x^2+2*y^2-2,x^2+x*y+y^2-2);
gens(gb(I))

----- Esercizio in classe -----
R = QQ[x,y,z,MonomialOrder=>Lex];
I = ideal(x^2+y^2+z^2-4, x^2+2*y^2-5, x*z-1);
gens(gb(I))

----- Esercizio in classe -----
R = QQ[x,y,z,MonomialOrder=>Lex];
I = ideal(x^10+1-y*x^5, x^2+1-z*x);
gens(gb(I))

----- Esercizio evoluta
R = QQ[t,x,y,MonomialOrder=>Lex];
I = ideal(2*y*(-t^4+1)-8*x*t*(t^2+1)-16*t*(t^2-1), 8*y*t^3-8*x*(3*t^2+1)-20*(3*t^2-1))
G = gens(gb(I))
G_(0,0)