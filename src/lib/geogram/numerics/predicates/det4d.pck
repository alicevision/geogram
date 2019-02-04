#include "kernel.pckh"

#ifdef DIM4

Sign predicate(det)(
    point(p0), point(p1), point(p2), point(p3) DIM
) {

   double m12 = p1_0*p0_1 - p0_0*p1_1;
   double m13 = p2_0*p0_1 - p0_0*p2_1;
   double m14 = p3_0*p0_1 - p0_0*p3_1;
   double m23 = p2_0*p1_1 - p1_0*p2_1;
   double m24 = p3_0*p1_1 - p1_0*p3_1;
   double m34 = p3_0*p2_1 - p2_0*p3_1;

   double m123 = m23*p0_2 - m13*p1_2 + m12*p2_2;
   double m124 = m24*p0_2 - m14*p1_2 + m12*p3_2;
   double m134 = m34*p0_2 - m14*p2_2 + m13*p3_2;
   double m234 = m34*p1_2 - m24*p2_2 + m23*p3_2;
        
   scalar Delta = (m234*p0_3 - m134*p1_3 + m124*p2_3 - m123*p3_3);
   return sign(Delta);
}

#endif
