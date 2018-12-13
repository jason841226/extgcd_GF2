#############################

p = 761; q61 = 765; q = 6*q61+1; w = 250
Zx.<x> = ZZ[]; R.<xp> = Zx.quotient(x^p-x-1)
Fq = GF(q); Fqx.<xq> = Fq[]; Rq.<xqp> = Fqx.quotient(x^p-x-1)
F3 = GF(3); F3x.<x3> = F3[]; R3.<x3p> = F3x.quotient(x^p-x-1)


##############################


# functions for lists

def conditional_swap( condition , a , b ) :
  if 0 == condition:
    return a , b
  else :
    return b , a


def list_shift(l,shift) :
  src_idx = max(-shift,0)
  dst_idx = max(shift,0)
  length = max(len(l)-abs(shift),0)
  new_l = [ l[0].parent()(0) ]* len(l)
  new_l[dst_idx:dst_idx+length] = l[src_idx:src_idx+length]
  return new_l


def list_pad_zero( al , len ):
  al_ref = [ al[0].parent(0) ]*(len)
  al = map( lambda x,y: y if x is None else x, al , al_ref )
  return al


def list_pad_front(al,nlen):
  assert nlen >= len(al)
  ral = list(al)
  ral[:0] = [al[0].parent(0)]*(nlen-len(al))
  return ral


def list_mul( l , c ):
  return [ x*c for x in l ]

def list_add( l1 , l2 ):
  return [ x+y for x,y in zip(l1,l2) ]


###################################

# functions for polynomials(lists)

def poly_in_deg( a , deg ):
  return list_pad_zero( a.list() , deg+1 )


def poly_rev( a , deg ):
  r = poly_in_deg(a,deg)
  rev_r = r[::-1]
  return rev_r


def poly_mul( a , b ):
  r = [a[0].parent(0)]*(len(a)+len(b)-1)
  for i in range(len(a)):
    r[i:i+len(b)] = list_add( r[i:i+len(b)] , list_mul(b,a[i]) )
  return r


###################################
##
##  This is not for a normal matrix. 
##  It is for the "transition matrix".
##
###################################

def mat_mul_vec( mat , vec ) :
  m_l = len(mat[0][0])
  f_m00 = poly_mul(mat[0][0],vec[0])
  g_m01 = poly_mul(mat[0][1],vec[1])
  f = list_add( f_m00 , g_m01 )
  f_m10 = poly_mul(mat[1][0],vec[0])
  g_m11 = poly_mul(mat[1][1],vec[1])
  g = list_add( f_m10 , g_m11 )
  ## XXX: assert 0 for truncated part.
  rf = f[m_l-1:]
  rg = list_pad_zero( g[m_l:] , len(vec[1]) )
  return  rf , rg 


def mat_mul( mat1 , mat2 ) :
  m1_len = len(mat1[0][0])
  m2_len = len(mat2[0][0])
  r_len = m1_len+m2_len
  r00_u = poly_mul(mat1[0][0],mat2[0][0])
  r00_d = poly_mul(mat1[0][1],mat2[1][0])
  r00 = list_add(list_pad_front(r00_u,r_len),list_pad_zero(r00_d,r_len))
  r01_u = poly_mul(mat1[0][0],mat2[0][1])
  r01_d = poly_mul(mat1[0][1],mat2[1][1])
  r01 = list_add(list_pad_front(r01_u,r_len),list_pad_zero(r01_d,r_len))
  r10_u = poly_mul(mat1[1][0],mat2[0][0])
  r10_d = poly_mul(mat1[1][1],mat2[1][0])
  r10 = list_add(list_pad_front(r10_u,r_len),list_pad_zero(r10_d,r_len))
  r11_u = poly_mul(mat1[1][0],mat2[0][1])
  r11_d = poly_mul(mat1[1][1],mat2[1][1])
  r11 = list_add(list_pad_front(r11_u,r_len),list_pad_zero(r11_d,r_len))
  return [ [r00,r01] , [ r10 , r11 ] ]

def deg1mat_mul( deg1mat , mat2 ) :
  inp_poly_len = len(mat2[0][0])
  np_len = inp_poly_len+1
  mat2_r0 = [list_pad_front(mat2[0][0],np_len),list_pad_front(mat2[0][1],np_len)]
  mat2_r1 = [list_pad_zero(mat2[1][0],np_len),list_pad_zero(mat2[1][1],np_len)]
  r0 , r0_ = conditional_swap( deg1mat[0][0][0] , mat2_r1 , mat2_r0 )   # use AND
  m1_10 = deg1mat[1][0][0]
  m1_11 = deg1mat[1][1][0]
  r10_u = list_mul(mat2[0][0],m1_10)
  r10_d = list_mul(mat2[1][0],m1_11)
  r11_u = list_mul(mat2[0][1],m1_10)
  r11_d = list_mul(mat2[1][1],m1_11)
  r10 = list_add(list_pad_front(r10_u,np_len),list_pad_zero(r10_d,np_len))
  r11 = list_add(list_pad_front(r11_u,np_len),list_pad_zero(r11_d,np_len))
  return [ r0 , [ r10 , r11 ] ]


##############################
#
#  Main Content Starts.
#
##############################



# Some variables for testing

def random_poly( deg ):
  r = Fqx.random_element()
  r = Fqx(r-r + Fq.random_element())
  for i in range(deg):
    r = Fqx( r*xq + Fq.random_element() )
  return r

a1 = random_poly(1)
b1 = random_poly(1)
c1 = random_poly(1)
a2 = random_poly(2)
b2 = random_poly(2)
c2 = random_poly(2)
a3 = random_poly(3)
b3 = random_poly(3)
c3 = random_poly(3)


##############################

# The standard GCD in the textbook

#output: gcd = s*f + t*g
def ext_gcd( f , g ) :
  if f.degree() < g.degree() :  f,g = g,f  #swap
  T = f.parent()
  lr = [ f , g ]
  ls = [ T(1) , T(0) ]
  lt = [ T(0) , T(1) ]
  #print "lr: [f,g] = ", lr,'\nls:',ls, "\nlt:", lt, '\n'
  print "lr[0]: " , lr[0] , "\tls[0]:", ls[0] , "\tlt[0]:", lt[0] 
  print "lr[1]: " , lr[1] , "\tls[1]:", ls[1] , "\tlt[0]:", lt[1] 
  print ""
  while lr[-1] != 0 :
    qq = lr[-2] // lr[-1]
    rr = lr[-2] % lr[-1]
    lr.append( rr )
    ls.append( ls[-2] - qq*ls[-1] )
    lt.append( lt[-2] - qq*lt[-1] )
    print "lr[-1]:",lr[-1] , "\tls[-1]:", ls[-1], "\tlt[-1]" , lt[-1]
  print "\nlr:", lr,'\nls:',ls, "\nlt:", lt, '\n'
  #normalize
  inv_n = lr[-2].leading_coefficient()^-1
  return lr[-2]*inv_n,ls[-2]*inv_n,lt[-2]*inv_n


# demo code:
# ext_gcd( a2, b2 )


#first variant: Performing operations as matrix-multiplication of 2x2 matrices.

def ext_gcd1( f , g ) :
  if f.degree() < g.degree() :  f,g = g,f  #swap
  Ta = f.parent()
  x = Ta.gen()
  lr = [ f , g ]
  ls = [Ta(1),Ta(0)]
  lt = [Ta(0),Ta(1)]
  print "[f,g]:" , lr,"\n",ls,"\n",lt , "\n"
  while lr[-1] != 0 :
    if lr[0].degree() > lr[1].degree() : 
       lr = lr[::-1]
       ls , lt = lt , ls
       print "op0: [0,1] <-- swap, choose g"
    else:
       print "op0: [1,0] <-- choose f"
    delta = lr[1].degree() - lr[0].degree()
    lc_f = lr[0].leading_coefficient()
    lc_g = lr[1].leading_coefficient()
    lr[1] = lr[1]*lc_f - lr[0]*lc_g*(x^delta)
    lt[0] = lt[0]*lc_f-ls[0]*lc_g*(x^delta)
    lt[1] = lt[1]*lc_f-ls[1]*lc_g*(x^delta)
    print "op1: [" , lc_f , "," , lc_g , "x^" , delta , "]"
    print "-->\n[f,g]:" , lr,"\n",ls,"\n",lt , "\n"
  #normalize
  inv_n = (lr[-2].leading_coefficient())^-1
  return lr[-2]*inv_n,map(lambda x:x*inv_n, ls )

# demo code:
# ext_gcd1( a2, b2 )


############################################################


# One 'division step' in the GCD.
# Trying to make it time-constant by applying 'conditional swap'

def divstep( delta , f , g ):
  coef_k = f[0].parent()
  kx.<x> = coef_k[] #kx = T_[0][0].parent()
  #x = kx.gen()
  #M2kx = T_.parent() 
  M2kx = MatrixSpace(kx.fraction_field(),2)   
  f = list_pad_zero( f , max(len(f),len(g)))
  g = list_pad_zero( g , max(len(f),len(g)))
  T0 = M2kx( (1,0,-g[0]/x,f[0]/x) )
  T1 = M2kx( (0,1,g[0]/x,-f[0]/x) )
  #conditional move
  c = 1 if delta > 0 and g[0] != 0 else 0
  rc , r1 = conditional_swap( c , [delta,f,g] , [-delta,g,f] )
  T0 , T1 = conditional_swap( c , T0 , T1 )
  # elimination
  rc[0] = 1 + rc[0]
  f = rc[1]
  g = rc[2]
  rg = [ f[0] * gi - g[0] * fi for fi,gi in zip(f,g) ]
  #rg.pop(0)
  rg = list_shift(rg,-1)
  return rc[0] , rc[1] , rg , T0

# demo code:
# Warning: The lists are in the reverse order of polynomials
#
# f = poly_rev( a2 , 2 )
# g = poly_rev( b2 , 2 )
# print "[f,g] = ", [ f, g ] , "\ndivstep(0,f,g): -->"
# divstep( 0 , f , g )
#


# The reversed Constant-time GCD

def ext_gcd2( a , b ):
  deg_a = a.degree()
  deg_b = b.degree()
  if deg_a < deg_b :
    a , b = b , a
    deg_a , deg_b = deg_b , deg_a
  if deg_a < 0 : return a
  #if deg_a == deg_b : 
  #  b = b*a.leading_coefficient() - a*b.leading_coefficient()
  #deg_b = deg_a - 1
  delta = deg_a - deg_b
  r1 = poly_rev( a , deg_a )
  r2 = poly_rev( b , deg_b )
  M2kx = MatrixSpace( a.parent().fraction_field() , 2 )
  T0 = M2kx(1)
  print delta , r1 , r2 , "\n", T0 , "\n"
  delta , r1 , r2 , T1 = divstep( delta , r1 , r2 )
  T0 = M2kx(T1)*T0
  print delta , r1 , r2 , '\n' , T1 , '\n' , T0 , '\n'
  for i in range(deg_a+deg_b-1) :
    delta , r1 , r2 , T1 = divstep( delta , r1 , r2 )
    T0 = M2kx(T1)*T0
    print delta , r1 , r2 , '\n' , T1 , '\n' , T0 , '\n'
  #normalize
  inv = r1[0]^-1
  r1 = map( lambda x: x*inv , r1 )
  T0[0] = T0[0]*inv
  print "output:"
  print delta , r1 , r2 , '\n' , T0 , '\n'
  return r1 , T0


# demo code:
# ext_gcd2( a2, b2 )
# ext_gcd2( a2*b1 , a1*b1 )
# ext_gcd1( a2*b1 , a1*b1 )


########################################################


# The "list" version of previous ext_ext2()


def divstep_list( delta , f , g ):
  coef_k = f[0].parent()
  T00 = [ [coef_k(1)] , [coef_k(0)] ]
  T10 = [ [-g[0]] , [f[0]]  ]
  T01 = [ [coef_k(0)] , [coef_k(1)] ]
  T11 = [ [g[0]]  , [-f[0]] ]
  #conditional move
  c = 1 if delta > 0 and g[0] != 0 else 0
  rd , rd_ = conditional_swap( c , delta , -delta )
  f  , g   = conditional_swap( c , f , g )
  T0 , T0_ = conditional_swap( c , [ T00, T10 ] , [ T01 , T11 ] )
  # elimination
  rd = rd + 1
  rg = [ f[0] * gi - g[0] * fi for fi,gi in zip(f,g) ]
  rg = list_shift( rg , -1 )
  return rd , f , rg , T0


# f = poly_rev( a2 , 2 )
# g = poly_rev( b2 , 2 )
# print "[f,g] = ", [ f, g ] , "\ndivstep(0,f,g): -->"
# divstep_list( 0 , f , g )


def extgcd_list( a , b ):
  deg_a = a.degree()
  deg_b = b.degree()
  if deg_a < deg_b :
    a , b = b , a
    deg_a , deg_b = deg_b , deg_a
  if deg_a < 0 : return a
  if deg_a == deg_b : 
    b = b*a.leading_coefficient() - a*b.leading_coefficient()
  deg_b = deg_a - 1
  delta = deg_a - deg_b
  r1 = poly_rev( a , deg_a )
  r2 = poly_rev( b , deg_b )
  lenlist = deg_a + 1
  r2 = list_pad_zero( r2 , lenlist )
  coef_k = r1[0].parent()
  print delta , r1 , r2 , "\n"
  delta , r1 , r2 , T1 = divstep_list( delta , r1 , r2 )
  T0 = T1    # T1*Identity = T1
  print delta , r1 , r2 , '\n' , T1 , '\n'
  for i in range(deg_a+deg_b-1) :
    delta , r1 , r2 , T1 = divstep_list( delta , r1 , r2 )
    T0 = deg1mat_mul( T1 , T0 )
    print delta , r1 , r2 , '\n' , T1 , '\n' , T0 , '\n'
  #normalize
  inv = r1[0]^-1
  r1 = map( lambda x: x*inv , r1 )
  T0[0][0] = list_mul( T0[0][0] , inv )
  T0[0][1] = list_mul( T0[0][1] , inv )
  print "output:"
  print delta , r1 , r2 , '\n' , T0 , '\n'
  return r1 , T0


# demo code:
# extgcd_list( a2*b1 , a1*b1 )


  

##########################################3


def jumpdivsteps(n,t,delta,f,g):
  print "jumpdivsteps(",n,t,delta,f,g,")"
  assert t>= n and n>=0
  f,g = f.truncate(t), g.truncate(t)
  print "truncated: f=",f , "  g=", g
  kx=f.parent()
  x=kx.gen()
  M2kx=MatrixSpace(kx.fraction_field(),2)
  if n == 0: return delta,f,g,(),M2kx(1)
  if n==1:
    c = 1 if delta > 0 and g[0] != 0 else 0
    Tc, Td = conditional_swap( c , M2kx((1,0,-g[0]/x,f[0]/x)), M2kx( (0,1,g[0]/x,-f[0]/x) ) )
    delta , d2 = conditional_swap( c , delta , -delta )
    fc , f2 = conditional_swap( c , f , g )
    gc , g2 = conditional_swap( c , f[0]*g-g[0]*f , g[0]*f-f[0]*g )
    return 1+delta,fc,kx(gc/x),(Tc,),Tc
  # other cases of n
  j = n//2
  print "j=",j
  delta,f1,g1,T1,P1 = jumpdivsteps(j,j,delta,f,g)
  print "P1:\n" , P1 
  f,g = P1*vector( (f,g) )
  print "P1x(f,g): f=" , f , "  g=", g
  f,g = kx(f).truncate(t-j),kx(g).truncate(t-j)
  print "truncated: f=",f , "  g=", g
  delta,f2,g2,T2,P2 = jumpdivsteps(n-j,n-j,delta,f,g)
  print "P2:\n" , P2
  f,g = P2*vector( (f,g) )
  print "P2x(f,g): f=" , f , "  g=", g
  f,g = kx(f).truncate(t-n+1),kx(g).truncate(t-n)
  print "truncated: f=",f , "  g=", g
  return delta,f,g,T2+T1,P2*P1

# demo code:
# f = a2.reverse( 2 )
# g = b2.reverse( 2 )
# delta, f , g , T2 , T1 = jumpdivsteps( 4 , 5 , 0 , f , g )
# print delta, "[f,g] = " , [f,g]
# print T2
# print T1


def jumpgcd(R0,R1):
  d = R0.degree()
  assert d>0 and d>R1.degree()
  f,g = R0.reverse(d), R1.reverse(d-1)
  delta,f,g,T,P = jumpdivsteps(2*d-1,3*d-1,1,f,g)
  return f.reverse(delta//2)/f[0]

# demo code:
# f = a2.reverse( 2 )
# g = b1.reverse( 1 )
# jumpgcd(f,g)


#########################################################

# The list version of previous JumpGCD

def jumpdivstep_mat(n,delta,f,g):
  print "jumpdivstep_mat(",n,delta,f,g,")"
  assert n>0
  ### truncate 
  f,g = f[:n], g[:n]
  print "truncated(",n,"): f=",f , "  g=", g
  if n==1:
    rd,f,g,T = divstep_list( delta , f , g )
    return rd,T
  # other cases of n
  j = n//2
  print "j=",j
  delta,P1 = jumpdivstep_mat(j,delta,f,g)
  print "P1:\n" , P1 
  f,g = mat_mul_vec( P1 , (f,g) )
  print "P1x(f,g): f=" , f , "  g=", g
  delta,P2 = jumpdivstep_mat(n-j,delta,f,g)
  print "P2:\n" , P2
  #f,g = mat_mul_vec( P2 , (f,g) )
  #print "P2x(f,g): f=" , f , "  g=", g
  #f,g = kx(f).truncate(t-n+1),kx(g).truncate(t-n)
  #print "truncated: f=",f , "  g=", g
  rP = mat_mul( P2 , P1 )
  print "ret: ", delta , rP
  return delta, rP

def jumpdivsteps_list(n,delta,f,g):
  print "jumpdivsteps_list(",n,delta,f,g,")"
  assert n>0
  if n==1:
    return divstep_list( delta , f , g )
  # other cases of n
  j = n//2
  print "j=",j
  delta,P1 = jumpdivstep_mat(j,delta,f,g)
  print "P1:\n" , P1 
  f,g = mat_mul_vec( P1 , (f,g) )
  print "P1x(f,g): f=" , f , "  g=", g
  delta,P2 = jumpdivstep_mat(n-j,delta,f,g)
  print "P2:\n" , P2
  f,g = mat_mul_vec( P2 , (f,g) )
  print "P2x(f,g): f=" , f , "  g=", g
  #f,g = kx(f).truncate(t-n+1),kx(g).truncate(t-n)
  print "truncated ???: f=",f , "  g=", g
  rP = mat_mul( P2 , P1 )
  return delta, f , g , rP


# demo code:
# f = poly_rev( a2 , 2 )
# g = poly_rev( b2 , 2 )
# jumpdivsteps_list( 5 , 0 , f , g )



def jumpgcd_list(R0,R1):
  d = R0.degree()
  d1 = R1.degree()
  assert d>0 and d>d1
  f = poly_rev(R0,d)
  g = poly_rev(R1,d1)
  g = list_pad_zero(g,len(f))
  delta,f,g,T = jumpdivsteps_list(2*d-1,1,f,g)
  #return delta,f,g,T
  gcd_deg = delta>>1
  gcd_len = gcd_deg + 1
  inv = f[0]^-1
  r_gcd = list_mul( f , inv )
  r_s = list_mul( T[0][0] , inv )
  r_t = list_mul( T[0][1] , inv )
  deg_s = d1 - gcd_deg -1
  deg_t = d - gcd_deg -1
  return gcd_len , r_gcd , r_s[gcd_deg*2:gcd_deg*2+deg_s+1] , r_t[gcd_deg*2:gcd_deg*2+deg_t+1]

# demo code:
# ff, gg = a1*b1*c1 , b2
# ext_gcd( ff , gg  )
# jumpgcd_list( ff , gg )