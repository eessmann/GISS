#define NOT_ZERO 1.e-30

/*-----------------------------------------------------* 
 *MYC - Mixed Youngs and Central Scheme                *
 *-----------------------------------------------------*/
/* 

Known problems: the index [1][1][1], i.e. the central cell
in the block, never occurs: neither in the central scheme
nor in Youngs' method. Therefore an isolated droplet will have
a normal with all components to zero. I took care of the
division-by-zero issue, but not of this one.

Ruben

*/
static void mycs(double c[3][3][3],double mxyz[3])
{ 
  double m1,m2,m[4][3],t0,t1,t2;
  int cn;

  /* write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
                             m00 X = m01 Y + m02 Z + alpha */
  m1 = c[0][1][0] + c[0][1][2] + c[0][0][1] + c[0][2][1] + 
       c[0][1][1];
  m2 = c[2][1][0] + c[2][1][2] + c[2][0][1] + c[2][2][1] + 
       c[2][1][1];
  t0 = fabs(m1-m2) + NOT_ZERO; 
  m[0][0] = (m1-m2)/t0;

  m1 = c[0][0][1]+ c[2][0][1]+ c[1][0][1];
  m2 = c[0][2][1]+ c[2][2][1]+ c[1][2][1];
  m[0][1] = 0.5*(m1-m2);

  m1 = c[0][1][0]+ c[2][1][0]+ c[1][1][0];
  m2 = c[0][1][2]+ c[2][1][2]+ c[1][1][2];
  m[0][2] = 0.5*(m1-m2);

  /* write the plane as: sgn(my) Y =  mx X +  mz Z + alpha 
                             m11 Y = m10 X + m12 Z + alpha */
  m1 = c[0][0][1] + c[0][2][1] + c[0][1][1];
  m2 = c[2][0][1] + c[2][2][1] + c[2][1][1];
  m[1][0] = 0.5*(m1-m2);

  m1 = c[1][0][0] + c[1][0][2] + c[2][0][1] + c[0][0][1] +
       c[1][0][1];
  m2 = c[1][2][0] + c[1][2][2] + c[2][2][1] + c[0][2][1] +
       c[1][2][1];
  t0 = fabs(m1-m2) + NOT_ZERO; 
  m[1][1] = (m1-m2)/t0;

  m1 = c[1][0][0]+ c[1][1][0]+ c[1][2][0];
  m2 = c[1][0][2]+ c[1][1][2]+ c[1][2][2];
  m[1][2] = 0.5*(m1-m2);

  /* write the plane as: sgn(mz) Z =  mx X +  my Y + alpha 
                             m22 Z = m20 X + m21 Y + alpha */

  m1 = c[0][1][0]+ c[0][1][2]+ c[0][1][1];
  m2 = c[2][1][0]+ c[2][1][2]+ c[2][1][1];
  m[2][0] = 0.5*(m1-m2);

  m1 = c[1][0][0]+ c[1][0][2]+ c[1][0][1];
  m2 = c[1][2][0]+ c[1][2][2]+ c[1][2][1];
  m[2][1] = 0.5*(m1-m2);

  m1 = c[0][1][0] + c[2][1][0] + c[1][0][0] + c[1][2][0] +
       c[1][1][0];
  m2 = c[0][1][2] + c[2][1][2] + c[1][0][2] + c[1][2][2] +
       c[1][1][2];
  t0 = fabs(m1-m2) + NOT_ZERO; 
  m[2][2] = (m1-m2)/t0;

  /* normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1 */
  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]) + NOT_ZERO;
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]) + NOT_ZERO;
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]) + NOT_ZERO;
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;

  /* choose among the three central scheme */ 
  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0 && t1 > t2)
    cn = 1;
  if (t2 > t0 && t2 > t1)
    cn = 2;

  /* Youngs-CIAM scheme */  
  m1 = c[0][0][0] + c[0][2][0] + c[0][0][2] + c[0][2][2] +
       2.*(c[0][0][1] + c[0][2][1] + c[0][1][0] + c[0][1][2]) +
       4.*c[0][1][1];
  m2 = c[2][0][0] + c[2][2][0] + c[2][0][2] + c[2][2][2] +
       2.*(c[2][0][1] + c[2][2][1] + c[2][1][0] + c[2][1][2]) +
       4.*c[2][1][1];
  m[3][0] = m1-m2;

  m1 = c[0][0][0] + c[0][0][2] + c[2][0][0] + c[2][0][2] +
       2.*( c[0][0][1] + c[2][0][1] + c[1][0][0] + c[1][0][2]) +
       4.*c[1][0][1];
  m2 = c[0][2][0] + c[0][2][2] + c[2][2][0] + c[2][2][2] +
       2.*(c[0][2][1] + c[2][2][1] + c[1][2][0] + c[1][2][2]) +
       4.*c[1][2][1];
  m[3][1] = m1-m2;

  m1 = c[0][0][0] + c[0][2][0] + c[2][0][0] + c[2][2][0] +
       2.*(c[0][1][0] + c[2][1][0] + c[1][0][0] + c[1][2][0]) +
       4.*c[1][1][0];
  m2 = c[0][0][2] + c[0][2][2] + c[2][0][2] + c[2][2][2] +
       2.*(c[0][1][2] + c[2][1][2] + c[1][0][2] + c[1][2][2]) +
       4.*c[1][1][2];
  m[3][2] = m1-m2;

  /* normalize the set (mx,my,mz): |mx|+|my|+|mz| = 1 */
  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]) + NOT_ZERO;
  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;

  /* choose between the previous choice and Youngs-CIAM */
  t0 = MAX(fabs(m[3][0]),fabs(m[3][1]));
  t0 = MAX(t0,fabs(m[3][2]));

  if (fabs(m[cn][cn]) > t0) 
    cn = 3;

  /* components of the normal vector */
  mxyz[0] = m[cn][0];
  mxyz[1] = m[cn][1];
  mxyz[2] = m[cn][2];

  return; 
}
