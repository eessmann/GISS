BEGIN {
  x1 = -3./4.;
  x2 = 1.;
  x = 1./4.;

  a1 = 1./(x1*x1 - x1*x2);
  a2 = - 1./(x1*x2 - x2*x2);
  a0 = - a1 - a2;

  b1 = - x2/(x1*x1 - x1*x2);
  b2 = x1/(x1*x2 - x2*x2);
  b0 = - b1 - b2;

  printf ("f(x)  = %g*p0 + %g*p1 + %g*p2\n",
	  x*x*a0 + x*b0 + 1,
	  x*x*a1 + x*b1,
	  x*x*a2 + x*b2);
  printf ("f'(x) = %g*p0 + %g*p1 + %g*p2\n",
	  2.*a0*x + b0,
	  2.*a1*x + b1,
	  2.*a2*x + b2);
}
