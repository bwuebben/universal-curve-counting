"""Independent quotient-algebra/resultant verification of Paper 2."""
import sympy as s
r,x,y,w=s.symbols('r x y w')
a=-2*(r-1)/(r+1)
b=-6*(r*r-1)/(r*r+r+1)
z=(x*x-a*x+3)/3
p=s.Poly(s.expand((x-b)*(x*x-a*x+3)+6*x),x,domain=s.QQ.frac_field(r))
c=-r**3/(r-1)**3-r**3*x/((r-1)*(r*r-1))-r**3*z/(r**3-1)
e=-18+2*(a+b)*x+(9+a*b-b*b)*z
c=s.Poly(c,x,domain=s.QQ.frac_field(r))
e=s.Poly(e,x,domain=s.QQ.frac_field(r))
beta=(r**3*e*c.invert(p)**2).rem(p)
print('p:',s.factor(p.as_expr()))
print('beta:',[s.factor(v) for v in beta.all_coeffs()])
q=s.Poly(s.resultant(p.as_expr(),y-beta.as_expr(),x),y)
# substitution r+1/r=2-w without importing the inspected script
quad=s.Poly(r*r-(2-w)*r+1,r,domain=s.QQ.frac_field(w))
def symmetrize(v):
    n,d=s.fraction(s.factor(v))
    np=s.Poly(n,r,domain=s.QQ.frac_field(w)).rem(quad)
    dp=s.Poly(d,r,domain=s.QQ.frac_field(w)).rem(quad)
    out=(np*dp.invert(quad)).rem(quad)
    assert out.degree()<=0
    return s.factor(out.as_expr())
co=[symmetrize(v) for v in q.all_coeffs()]
print('charpoly:',co)
tr2=s.expand(co[1]**2-2*co[2])
print('trace beta squared:',tr2)
print('trace beta squared factored:',s.factor(tr2))
print('r=2 charpoly:',[s.factor(v.subs(r,2)) for v in q.all_coeffs()])
print('r=2 tr2:',s.factor((q.all_coeffs()[1]**2-2*q.all_coeffs()[2]).subs(r,2)))
print('r coefficient at -6:',s.expand(tr2).coeff(w,6))
