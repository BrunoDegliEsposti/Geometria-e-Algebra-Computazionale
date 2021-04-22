---attempt for Sudoku
restart
K2=toField (QQ[q]/ideal(1+q^3+q^6))
R=K2[x_(0,0)..x_(8,8)]---81 variables
V=x->x^9-1
E=(x,y)->x^8+x^7*y+x^6*y^2+x^5*y^3+x^4*y^4+x^3*y^5+x^2*y^6+x*y^7+y^8
I=ideal(0_R)
for i from 0 to 8 do for j from 0 to 8 do I=I+ideal(V(x_(i,j)))
for i from 0 to 8 do for j from 0 to 7 do for k from j+1 to 8 do I=I+ideal(E(x_(i,j),x_(i,k)))---rows
for i from 0 to 8 do for j from 0 to 7 do for k from j+1 to 8 do I=I+ideal(E(x_(j,i),x_(k,i)))---columns
for i from 0 to 2 do for j from 0 to 2 do for k from 0 to 7 do for p from k+1 to 8 do
I=I+ideal(E(x_(3*i+floor(k/3),3*j+(k%3)),x_(3*i+floor(p/3),3*j+(p%3))))--blocks
I=I+ideal(x_(0,0)-q^8,x_(0,2)-q^6,x_(0,6)-q^7,x_(0,8)-1,x_(1,2)-q^4,x_(1,6)-q^2,
    x_(2,1)-q^7,x_(2,3)-q^4,x_(2,4)-1,x_(2,7)-q,x_(3,5)-q^3,x_(3,6)-q^5,
    x_(4,2)-q^2,x_(4,4)-q^6,x_(4,6)-q^8,x_(5,2)-q^1,x_(5,3)-1,
    x_(6,1)-q^4,x_(6,4)-q^2,x_(6,5)-q^7,x_(6,7)-q^5,x_(7,2)-1,x_(7,6)-q^6,
    x_(8,0)-q^5,x_(8,2)-q^7,x_(8,6)-q^1,x_(8,8)-q^3);
betti I
time gens gb I
toString oo
----this is the solution, after 1620 seconds = 27 minutes
 matrix {{x_(8,8)-q^3, x_(8,7)-q^2, x_(8,6)-q, x_(8,5)-1, x_(8,4)-q^4, x_(8,3)+q^3+1, x_(8,2)+q^4+q, x_(8,1)+q^5+q^2, x_(8,0)-q^5,
      x_(7,8)-q^4, x_(7,7)+q^4+q, x_(7,6)+q^3+1, x_(7,5)-q^5, x_(7,4)-q^3, x_(7,3)+q^5+q^2, x_(7,2)-1, x_(7,1)-q^2, x_(7,0)-q, x_(6,8)+q^5+q^2,
      x_(6,7)-q^5, x_(6,6)-1, x_(6,5)+q^4+q, x_(6,4)-q^2, x_(6,3)-q, x_(6,2)-q^3, x_(6,1)-q^4, x_(6,0)+q^3+1, x_(5,8)+q^4+q, x_(5,7)+q^3+1,
      x_(5,6)-q^4, x_(5,5)-q^2, x_(5,4)+q^5+q^2, x_(5,3)-1, x_(5,2)-q, x_(5,1)-q^5, x_(5,0)-q^3, x_(4,8)-q, x_(4,7)-q^3, x_(4,6)+q^5+q^2,
      x_(4,5)-q^4, x_(4,4)+q^3+1, x_(4,3)-q^5, x_(4,2)-q^2, x_(4,1)-1, x_(4,0)+q^4+q, x_(3,8)-q^2, x_(3,7)-1, x_(3,6)-q^5, x_(3,5)-q^3, x_(3,4)-q,
      x_(3,3)+q^4+q, x_(3,2)+q^5+q^2, x_(3,1)+q^3+1, x_(3,0)-q^4, x_(2,8)+q^3+1, x_(2,7)-q, x_(2,6)-q^3, x_(2,5)+q^5+q^2, x_(2,4)-1, x_(2,3)-q^4,
      x_(2,2)-q^5, x_(2,1)+q^4+q, x_(2,0)-q^2, x_(1,8)-q^5, x_(1,7)+q^5+q^2, x_(1,6)-q^2, x_(1,5)+q^3+1, x_(1,4)+q^4+q, x_(1,3)-q^3, x_(1,2)-q^4,
      x_(1,1)-q, x_(1,0)-1, x_(0,8)-1, x_(0,7)-q^4, x_(0,6)+q^4+q, x_(0,5)-q, x_(0,4)-q^5, x_(0,3)-q^2, x_(0,2)+q^3+1, x_(0,1)-q^3,
      x_(0,0)+q^5+q^2}}
----for example x_(8,3)=-1-q^3=q^6
---x_(2,1)=-q-q^4=q^7...
codim I, degree I

------------------------------
----attempt removing entries
restart
K2=toField (QQ[q]/ideal(1+q^3+q^6))
R=K2[x_(0,0)..x_(8,8)]---81 variables
V=x->x^9-1
E=(x,y)->x^8+x^7*y+x^6*y^2+x^5*y^3+x^4*y^4+x^3*y^5+x^2*y^6+x*y^7+y^8
I=ideal(0_R)
for i from 0 to 8 do for j from 0 to 8 do I=I+ideal(V(x_(i,j)))
for i from 0 to 8 do for j from 0 to 7 do for k from j+1 to 8 do I=I+ideal(E(x_(i,j),x_(i,k)))---rows
for i from 0 to 8 do for j from 0 to 7 do for k from j+1 to 8 do I=I+ideal(E(x_(j,i),x_(k,i)))---columns
for i from 0 to 2 do for j from 0 to 2 do for k from 0 to 7 do for p from k+1 to 8 do
I=I+ideal(E(x_(3*i+floor(k/3),3*j+(k%3)),x_(3*i+floor(p/3),3*j+(p%3))))--blocks
---remove x_(0,2)
I=I+ideal(x_(0,0)-q^8,x_(0,6)-q^7,x_(0,8)-1,x_(1,2)-q^4,x_(1,6)-q^2,
    x_(2,1)-q^7,x_(2,3)-q^4,x_(2,4)-1,x_(2,7)-q,x_(3,5)-q^3,x_(3,6)-q^5,
    x_(4,2)-q^2,x_(4,4)-q^6,x_(4,6)-q^8,x_(5,2)-q^1,x_(5,3)-1,
    x_(6,1)-q^4,x_(6,4)-q^2,x_(6,5)-q^7,x_(6,7)-q^5,x_(7,2)-1,x_(7,6)-q^6,
    x_(8,0)-q^5,x_(8,2)-q^7,x_(8,6)-q^1,x_(8,8)-q^3);
betti I
---following is now long
time gens gb I
