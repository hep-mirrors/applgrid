#ifndef __AMP_MACROS_H__
#define __AMP_MACROS_H__ 1



#define S(i, j) ((_M_s->ptr)[p##i][p##j])
#define A(i, j) ((_M_a->ptr)[p##i][p##j])
#define B(i, j) ((_M_b->ptr)[p##i][p##j])

#define Sij(i, j) ((_M_s->ptr)[i][j])
#define Aij(i, j) ((_M_a->ptr)[i][j])
#define Bij(i, j) ((_M_b->ptr)[i][j])

#define S3(i, j, k) (S(i,j)+S(i,k)+S(j, k))

#define A3(i, j, k) (A(i,j)*B(j,k))
#define B3(i, j, k) (B(i,j)*A(j,k))

#define A4(i, j, k, l) (A(i,j)*B(j,k)*A(k,l))
#define B4(i, j, k, l) (B(i,j)*A(j,k)*B(k,l))

#define A5(i, j, k, l, m) (A(i,j)*B(j,k)*A(k,l)*B(l,m))
#define B5(i, j, k, l, m) (B(i,j)*A(j,k)*B(k,l)*A(l,m))

#define A6(i, j, k, l, m, n) (A(i,j)*B(j,k)*A(k,l)*B(l,m)*A(m,n))
#define B6(i, j, k, l, m, n) (B(i,j)*A(j,k)*B(k,l)*A(l,m)*B(m,n))

#define C(i, j, k, l) (A(i,j)*B(j,l)+A(i,k)*B(k,l))
#define D(i, j, k, l) (B(i,j)*A(j,l)+B(i,k)*A(k,l))


#define D12  (S(1,2)-S(3,4)-S(5,6))
#define D34  (S(3,4)-S(1,2)-S(5,6))
#define D56  (S(5,6)-S(1,2)-S(3,4))
#define Dt23 (S(2,3)-S(5,6)-S(1,4))
#define Dt14 (S(1,4)-S(5,6)-S(2,3))
#define Dt56 (S(5,6)-S(2,3)-S(1,4))


#define G3  (pow(D56, 2)  - 4.0*S(1,2)*S(3,4))
#define Gt3 (pow(Dt56, 2) - 4.0*S(2,3)*S(1,4))


#define MASS_OF_TOP  175.0
#define mtsq         30625



#define R(i, j, k, l)   (A(i,k)*A(j,l)/(A(i,l)*A(j,k)))

#endif
