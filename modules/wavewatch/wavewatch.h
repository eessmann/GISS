/*
 * These types must match the lengths of wavewatch's REAL and INTEGER types
 */
typedef double REAL;
typedef long long INTEGER;

/**
 * Initialisation of wavewatch
 */
extern void __gfsw3init__gfsw3_init (void);

/**
 * imported from w3srcemd.f90:W3SRCE
 * @IX,@IY: Discrete grid point counters (hopefully unused)
 * @IMOD: Model number (always 1)
 * @SPEC: Spectrum (action) in 1-D form (NK*NTH elements): INPUT/OUTPUT
 * @ALPHA: Nondimensional 1-D spectrum (NK)                OUTPUT
 * @WN1: Discrete wavenumbers (NK)                         INPUT
 * @CG1: Discrete group velocities (NK)                    INPUT
 * @DEPTH: Depth                                           INPUT
 * @U10ABS: Wind speed at reference height                 INPUT
 * @U10DIR: Wind direction at reference height             INPUT
 * @USTAR: Friction velocity                               INPUT/OUTPUT
 * @USTDIR: Friction velocity direction                    OUTPUT (maybe model dependent?)
 * @EMEAN: Mean energy                                     OUTPUT (maybe model dependent?)
 * @FMEAN: Mean frequency                                  OUTPUT (maybe model dependent?)
 * @WMEAN: Mean wavenumber                                 OUTPUT (maybe model dependent?)
 * @AMAX: Maximum energy                                   OUTPUT
 * @FPI: Peak-input frequency                              INPUT/OUTPUT
 * @CD: Drag coefficient                                   OUTPUT (maybe model dependent?)
 * @Z0: Roughness length                                   OUTPUT (maybe model dependent?)
 * @DTDYN: Average dynamic time step                       OUTPUT
 * @FCUT: Cut-off frequency for tail                       OUTPUT
 * @DTG: Global time step                                  INPUT
 */
extern void __w3srcemd__w3srce (INTEGER * IX, INTEGER * IY, INTEGER * IMOD,
				REAL * SPEC,
				REAL * ALPHA,
				REAL * WN1, 
				REAL * CG1, 
				REAL * DEPTH, 
				REAL * U10ABS, 
				REAL * U10DIR, 
				REAL * USTAR,
				REAL * USTDIR,
				REAL * EMEAN,
				REAL * FMEAN,
				REAL * WMEAN,
				REAL * AMAX,
				REAL * FPI,
				REAL * CD,
				REAL * Z0,
				REAL * DTDYN,
				REAL * FCUT,
				REAL * DTG);
