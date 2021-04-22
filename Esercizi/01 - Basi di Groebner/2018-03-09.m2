R = QQ[x,y,MonomialOrder=>Lex]
I = ideal(x^2+y^2-1,x*y)
for i from 1 to 10 do print((x^i+y^i)%I)
f = random(8,R)
f%I
