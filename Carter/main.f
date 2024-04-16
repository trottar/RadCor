	Program StrucFunc_Test
C	Program to generate f1(x,Q2), f2(x,Q2) etc. using proper low-Q2 behavior
C	
C	4/2/97 SEK
C
C	Use new NMC parametrization for F2 and
C	and SLAC 1998 fit for R
C	
C	Rewritten to call subroutines in "models.f" for all
C	calculations
c
C       Final H version 10/26/1999 SEK
C
C	Added necessary stuff for deuteron 12/2/1999 SEK
C
C 	Changed to allow calling of new Hall C parametrization for F1, Fl 23-Feb-2004 SEK
C
C	Updated to refer only to the newest models 15-Jun-2009 SEK

	implicit none

	include 'instruct.inc'

	real*8 amp/0.939/ ! xz
        character*200 fname  ! xz
	
Csk	REAL*8 VSTRT(37),STP(37),MY_LO(37),MY_HIGH(37),XVAL(37)
Csk	INTEGER NPRM(37)
Csk	CHARACTER*10 PNAM(37)

	real*8 Q2, Nu, X, R, F2, F1, qsq, W2, DR,E
	real*8 TH, G1, G2, A1, A2
	real*8 fn,fnerr,xval1(41),xvall(41),w2inc,temp(4)
        real*8 mp,mp2,pi,alpha

	character*1 Answer, T, nT /'N'/, pT /'P'/ !Target Type
	logical dummy, suppress !, newhallc
	integer model/12/, i, imax, npts, sf, YoniIndex
	integer pSFChoice, nSFChoice, pAsymChoice, nAsymChoice
	integer pIPOL, nIPOL, pIA1, nIA1
	
        REAL*8 Z_targ, A_targ

!     xz: new variables needed to calculate cross sections
        REAL*8 GAM2,GAMMA2,EP,Y,ROOTQ2
        REAL*8 SIN2,TAN2
        REAL*8 EPSILON,ETA,ZETA,ucD,lcD
        REAL*8 nbarn
        REAL*8 Apar,Aperp
        REAL*8 sig_mott,sig_unpol,sig_par,sig_perp
!     xz: end of new variables
       
!     flags for SF choices
	common /YoniDeut/ pSFChoice, nSFChoice, pAsymChoice, nAsymChoice, 
     *       pIPOL, nIPOL, pIA1, nIA1

        integer ifile
        
	NMC = .TRUE.
	suppress = .FALSE.
	BODEK1 = .FALSE.
	NORES = .FALSE.
	ERROR = .FALSE.
	IPOL = 1
csk	IPOL = 1: standard version of A2 in the DIS region
csk	IPOL = 2: extra term g2tw3 in A2 in DIS
csk	IPOL = 4: g2 = 0 in DIS
csk	IPOL = 5: A2 = 0 in DIS
csk	IPOL = 7: A2 = Soffer bound in DIS
	IPOLRES = 1
	IA1 = 4 ! Use to determine specific A1/A2 DIS model used
csk     IA1 = 4: sek_me (2006/7 combined Doaa/sek fit to world data including EG1)
csk     IA1 = 5: sek_me; add correlated fit error to A1
csk     IA1 = 6: sek_me; subtract correlated fit error to A1
csk     IA1 = 7: sek_2021 (New version by Pushpa)
csk     IA1 = 8: sek_2021 + dA1
csk     IA1 = 9: sek_2021 - dA1

	AsymChoice = 11 ! use to determine specific A1/A2 models used in the RR
csk    11: New Standard Resonance Model 2008-9 (c) Guler/Kuhn
csk    12: Preliminary version v1 of A1 model - 2008-9 (c) Guler/Kuhn 
csk    13: Old Standard for A1 2008-9 (c) Guler/Kuhn
csk    14: Alternative A2 resonance model: A2_MAID
csk    15: Old version of A2 resonance model 

	SFChoice = 23
csk    10: 2007 version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2007, HERMES
csk    11: 2009 version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2009, HERMES
csk    12: Same version as 11, but with errors added to F2 (and proportionally F1)
csk    13: Same version as 11, but with errors subtracted from R (F2 unchanged)
csk    14: Same version as 11, but without tabulated R values in RR substituted
csk    15-16: Analog to 12-13, but without tabulated R values in RR substituted
csk    17: Newest kludge for Rd, F2d from Eric Christy December 2011. Should be only for D
csk    18-19: Analog to 15-16 for new kludge
csk	20: New version by Christie/Kalantarians 2014
csk	21-22 corresponding error estimates
csk 23-25: New version F1F221
	
Cxz	open (unit=22, file='f1.in',status='old')
Cxz	open (unit=11, file='f1.out',status='replace')
	
1009	write(6,*) ' Enter Target type (P,N,3 or D - please capitalize):'
	read(5,111) T
111	format(a1)

	 if ((T .ne. 'P').and.(T .ne. 'D').and.(T. ne. 'N')
     *       .and.(T .ne. '3')) goto 1009


       if(T .eq. 'P') then
         Z_targ = 1.0D0
         A_targ = 1.0D0
       elseif(T .eq. 'D') then
         Z_targ = 1.0D0
         A_targ = 2.0D0
       elseif(T .eq. 'N') then
         Z_targ = 0.0D0
         A_targ = 1.0D0
       elseif(T .eq. '3') then
         Z_targ = 2.0D0
         A_targ = 3.0D0
       else
         write(6,*) ' Wrong Target ', T
         goto 1009
       endif
         
1523	write(6,*) ' Enter smearing preference: ',
     >      '0 = none, 1 = quasielastic only, 2 = all, 3 = inelastic only.'
      	read(5,112) YoniIndex
112	format(I1)
	if (YoniIndex.gt.3 .or. YoniIndex.lt.0) goto 1523

Cxz  Yoni code is only implemented for D
        if ((T.eq.'3').and.(YoniIndex.ne.0)) then
          print *,'smearing is implemented only for D. ',
     >         'Using D-smearing for 3He!'
!          YoniIndex=0
        endif
Cxz  end
	
	if ((T.eq.pT).or.(T.eq.nT).or.(YoniIndex.lt.1)) then
1245	  write(6,*) ' What value of IPOL do you want use? 1-7'
	  read(5,*) IPOL
	  if (IPOL.lt.1 .or. IPOL.gt.7) goto 1245
	
1246	  write(6,*) ' What value of IA1 do you want use? 4-9'
	  read(5,*) IA1
	  if (IA1.lt.4 .or. IA1.gt.9) goto 1246
	
1247	  continue
	  if (T.eq. 'D' .and. (YoniIndex.lt.1)) then
	    write(6,*) ' What value for SFChoice? 10-19'
	    read(5,*) SFChoice
	    if (SFChoice.lt.10 .or. SFChoice .gt. 19) goto 1247
	  else
	    write(6,*) ' What value for SFChoice? 11-25'
	    read(5,*) SFChoice
	    if (SFChoice.lt.10 .or. SFChoice .gt. 25) goto 1247
	  endif
	  
	
1248	  write(6,*) ' What value for AsymChoice? 11-13'
	  read(5,*) AsymChoice
	  if (AsymChoice.lt.11 .or. AsymChoice .gt. 13) goto 1248
	  
	else
	
	  write(6,*) ' Enter model parameters for proton'
2245	  write(6,*) ' What value of IPOL do you want use? 1-7'
	  read(5,*) pIPOL
	  if (pIPOL.lt.1 .or. pIPOL.gt.7) goto 2245
	
2246	  write(6,*) ' What value of IA1 do you want use? 4-9'
	  read(5,*) pIA1
	  if (pIA1.lt.4 .or. pIA1.gt.9) goto 2246
	
2247	  write(6,*) ' What value for SFChoice? 10-25'
	  read(5,*) pSFChoice
	  if (pSFChoice.lt.10 .or. pSFChoice .gt. 25) goto 2247
	
2248	  write(6,*) ' What value for AsymChoice? 11-13'
	  read(5,*) pAsymChoice
	  if (pAsymChoice.lt.11 .or. pAsymChoice .gt. 13) goto 2248
	  
	  write(6,*) ' Enter model parameters for neutron'
3245	  write(6,*) ' What value of IPOL do you want use? 1-7'
	  read(5,*) nIPOL
	  if (nIPOL.lt.1 .or. nIPOL.gt.7) goto 3245
	
3246	  write(6,*) ' What value of IA1 do you want use? 4-6'
	  read(5,*) nIA1
	  if (nIA1.lt.4 .or. nIA1.gt.6) goto 3246
	
3247	  write(6,*) ' What value for SFChoice? 10-25'
	  read(5,*) nSFChoice
	  if (nSFChoice.lt.10 .or. nSFChoice .gt. 25) goto 3247
	
3248	  write(6,*) ' What value for AsymChoice? 11-13'
	  read(5,*) nAsymChoice
	  if (nAsymChoice.lt.11 .or. nAsymChoice .gt. 13) goto 3248
	  
	endif
	
	write(6,*) ' Include Resonances? [Y]/n'
	read(5,111) Answer
	if ((Answer.eq.'N').or.(Answer.eq.'n')) then
	  NORES = .TRUE.
C	this makes h2mod use only the nonres. = background terms for sigma
	  suppress = .TRUE.
c       this skips h2mod entirely, using a DIS extrapolation instead.
	endif

C
C     xz: here is the main loop, set up for calculating structure functions
C     vs. x for fixed Q2. This can be changed to other loops or single settings
C
C     an input of E is added to calculate kinematic factors needed by Apar and Aperp
C     if only calculating F1, F2, G1, G2, A1, A2, then E is not needed
C
        E = 10.384
        print *,'Enter E beam energy value in GeV'
        read *,E
        
 3249   print *,'Enter Q2 value in GeV^2'
	read *,Q2
	
	if (Q2.lt.10) then
	  write(fname,'("outfiles/sig",a1,"_sm",I1,"_ipol",I1,"IA1",I1,"_SF",
     *      I2,"_AC",I2,"_Q2_",F5.3,".out")')
     *      T,YoniIndex,IPOL,IA1,SFChoice,AsymChoice,Q2
	else if (Q2.lt.100) then
	  write(fname,'("outfiles/sig",a1,"_sm",I1,"_ipol",I1,"IA1",I1,"_SF",
     *      I2,"_AC",I2,"_Q2_",F6.3,".out")')
     *         T,YoniIndex,IPOL,IA1,SFChoice,AsymChoice,Q2
        else
          print *,'Q2 too high, please check'
          goto 3249
	endif
        ifile=NextUn()
	open (unit=ifile, file=fname,status='replace')
	
	print *,'writing output to ',fname
        
Cxz: uncomment the following if reading W2, Q2 grid from an input file
C	do i = 1, 200000
C     read(22,*,end=99) W2, Q2
C	  Nu = (W2-amp**2+Q2)/2.0/amp
C	  qsq = Q2 + Nu*Nu
C     X = Q2/2.0/amp/Nu

Cxz: uncomment the following if reading X,Q2 grid from an input file
C	do i = 1, 200000        
C	read(22,*,end=99) X, Q2

C     now use loop for X at constant Q2
        X=0.01
        
        do i = 1, 200000
          X=X+0.001

        Nu = Q2/2.0/amp/X
        qsq = Q2 + Nu*Nu
        W2=amp**2. + 2.0*amp*Nu - Q2

        if (W2.lt.1.07**2.) goto 99 ! quit X loop

        if (YoniIndex .eq. 0)then                         
          call F1F2new (x, Q2, T, F1, F2, R, suppress)
C          print *,'check:',X,Q2,T,F1,F2,R,
C     >         (4.D0*amp*amp*X*X/Q2+1.D0)*F2/2.0D0/X/F1-1
C	R is here only used as placeholder for sign/sigp
          TH = 25.0 D00         ! Arbitrary, not really needed
          ! xz: note: TH is used in G1G2new, weird
	    call G1G2new(x, Q2, TH, T, G1, G2, A1, A2)
	    if (F1 .gt. 0. .and. Nu .gt. 0. .and. Q2 .gt. 0.) then
	      R = F2/F1/2.0/x*(1+Q2/Nu/Nu)-1.0
	    else
	      R = 0.0
	    endif          
C	  else if (T .ne. 'P' .and. T .ne. 'N') then ! intended for any non-zero smearing of deuterium
        else if (T .eq. 'D') then ! DSF only works for deuteron, not 3He

          call DSFs(X,Q2,YoniIndex,F1,F2,R,G1,G2,A1,A2,suppress)

C     xz's temporary fix
        elseif (T .eq. '3') then ! Yoni works only for smearing of P, N in deuterium, not 3He
          call DSFs3(X,Q2,YoniIndex,F1,F2,R,G1,G2,A1,A2,suppress)
Cxz end of XZ's fixes
          
        else ! P or N, call Yoni directly for D-smeared SF
          call Yoni(X,Q2,T,YoniIndex,F1,F2,R,G1,G2,A1,A2,suppress)
        endif
!        write(ifile,1234) Q2,W2,X,F1,F2,R,A1,A2,G1,G2
! 1234  format(' ',10f12.6)
        
! xz: now calculate Apar and Aperp from A1 and A1
! some of the varaibles below are already defined, redefine just to be sure
       GAM2 = 4.D0*amp*amp/Q2
       GAMMA2 = GAM2*X*X ! gamma^2
       NU = Q2/(2.D0*amp*X)
       EP = E - NU
       Y = NU/E
       if (abs(1-Q2/(2.*E*EP)).le.1.0) then
         TH = acos(1-Q2/(2.*E*EP)) ! Q2=2EE'(1-costh)
       
         print *,'calculating for E,EP,Q2,X=',E,EP,Q2,X
         
       SIN2 = (DSIN(TH/2.D0))**2
       TAN2 = SIN2/(1.D0 - SIN2)
       ROOTQ2 = DSQRT(Q2)
         
       EPSILON = 1.D0/(1.D0 + 2.D0*TAN2*(1.D0 + NU*NU/Q2)) ! Eq.(1.50) xz thesis
       ETA = EPSILON*ROOTQ2/(E - EP*EPSILON)  ! Eq.(1.51)
       ZETA = ETA*(1.D0 + EPSILON)/(2.D0*EPSILON) ! Eq.(1.52)

       R = F2/F1/2.0/X*(1+Q2/NU/NU)-1.0 ! sig_L/sig_T Eq.(1.20) and inverting (1.18)
       
       ucD = (1-(1-Y)*EPSILON)/(1+EPSILON*R) ! Eq.(1.49) uppercase D
       lcD = ucD * SQRT(2*EPSILON/(1+EPSILON)) ! Eq.(1.53) lowercase d

       Apar = ucD * (A1 + ETA * A2)
       Aperp = lcD * (A2 - ZETA * A1)

       alpha=7.3e-3             !  alpha = 1/137
       nbarn=0.389e6            ! barn: (1 GeV)**-2 = 0.389e-3 barn
c     Calculate Mott cross section in nbarn/(sr)
       sig_mott = (((alpha * DCOS(TH/2.)/2.*E*SIN2))**2)*nbarn
       sig_unpol = sig_mott * (F2/NU + 2*F1*TAN2/amp) ! unpol cross section in nbar/sr/GeV
       sig_par = sig_unpol * Apar
       sig_perp = sig_unpol * Aperp

c     xz: here can write whatever needed to the output file
       
C       write(ifile,1234) Q2,W2,X,F1,F2,R,A1,A2,G1,G2 ! structure functions
C 1234  format(' ',12f12.6)
       write(11,1234) Q2,W2,X,A1,A2,Apar,Aperp,sig_unpol,sig_par,sig_perp ! observables
       
       else
         print *,'invalid setting for E,EP,Q2,X=',E,EP,Q2,X
       endif
       
 1234  format(' ',9f12.6,3G12.4)
	enddo
99	continue
CXZ	close (22)
	close (ifile)
	Stop
	end

