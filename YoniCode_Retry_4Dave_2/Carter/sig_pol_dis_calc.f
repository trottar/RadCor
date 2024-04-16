c---------------------------------------------------------------
      subroutine sig_pol_dis_calc(iA, iZ, X, QSQ, G1F1, G2F1)
      implicit none
C      INCLUDE 'radcon.inc'
      include 'instruct.inc'
      real*8 QSQ,x,mp
      real*8 F1p,F1n,F2p,F2n,g1p,g2p,g1n,g2n,g1f1p,g1f1n,g2f1p,g2f1n
      real*8 R,F1,F2,G1,G2,A1,A2,A1p,A1n,A2p,A2n
      real*8 THdum
      real*8 g1f1,g2f1
      real*8 gam, gam2
      real*8 epsilon,rootq2,eta,zeta,ucD,lcD
      CHARACTER*1 TARG, TARGP /'P'/, TARGN /'N'/
      integer iA,iZ,yoni_index
      integer pIPOL,pIA1,pAsymchoice,pSFchoice
      integer nIPOL,nIA1,nAsymchoice,nSFchoice
      logical suppress
      common /YoniDeut/ pSFChoice, nSFChoice, pAsymChoice, nAsymChoice, 
     *     pIPOL, nIPOL, pIA1, nIA1

      parameter(mp=0.938272)

      logical DEBUG/.false./

C      print *,'enable printing by setting DEBUG to .true.'
      if (DEBUG) write(6,'("check kine:",2F8.4)')X,QSQ
      
      suppress=.false.
      yoni_index=2

      IPOL=1
      IA1=4
      Asymchoice=11
      SFchoice=23

      pIPOL = 1
      pIA1 = 4
      pAsymchoice=11
      pSFchoice=23
      
      nIPOL = 1
      nIA1 = 4
      nAsymChoice=11
      nSFchoice=23

      if (yoni_index.eq.0) then 
        
C     XZ calling free proton and free neutron SF directly
C     this would differ only slightly from Yoni output for most kinematic
C     points, except very high x, where the two are very different
C     
        TARG=TARGP
        call F1F2new(X, QSQ, TARG, F1p, F2p, R, suppress)
C     R is here only used as placeholder for sign/sigp
        THdum = 25.0 D00        ! Arbitrary, not really needed
        call G1G2new(X, QSQ, THdum, TARG, G1p, G2p, A1p, A2p)
        if (DEBUG) write(6,'("check free proton newSF:",2G12.4,F10.6)')F1p,G1p,G1p/F1p
          
        TARG=TARGN
        call F1F2new(X, QSQ, TARG, F1n, F2n, R, suppress)
C     R is here only used as placeholder for sign/sigp
        THdum = 25.0 D00        ! Arbitrary, not really needed
        call G1G2new(X, QSQ, THdum, TARG, G1n, G2n, A1n, A2n)
        if (DEBUG) write(6,'("check free neutron newSF:",2G12.4,F10.6)')F1n,G1n,G1n/F1n
        F1 = 2.*F1p + F1n
        F2 = 2.*F2p + F2n
        G1 = -2.*0.027*G1p + 0.87*G1n
        G2 = -2.*0.027*G2p + 0.87*G2n

        G1F1=0
        G2F1=0
        if (F1.ne.0) then
          G1F1 = G1/F1
          G2F1 = G2/F1
        endif
      
        if (DEBUG) write(6,'("check free-N he3 newSF:",4F8.4)')F1,G1,G1F1,G2F1
        goto 999

      else ! call Yoni for smearing. It will apply deuteron smearing to free p and n SF
        TARG=TARGP
        call Yoni(X,QSQ,TARG,yoni_index,F1p,F2p,R,G1p,G2p,A1p,A2p,
     >       suppress)
        if (DEBUG) write(6,'("check proton:",2G12.4,F10.6)')F1p,G1p,G1p/F1p

        TARG=TARGN
        call Yoni(X,QSQ,TARG,yoni_index,F1n,F2n,R,G1n,G2n,A1n,A2n,
     >       suppress)
        
        if (DEBUG) write(6,'("check neutron:",2G12.4,F10.6)')F1n,G1n,G1n/F1n

        F1 = 2.*F1p + F1n
        F2 = 2.*F2p + F2n
        G1 = -2.*0.027*G1p + 0.87*G1n
        G2 = -2.*0.027*G2p + 0.87*G2n
        
        G1F1=0
        G2F1=0
        if (F1.ne.0) then
          G1F1 = G1/F1
          G2F1 = G2/F1
        endif
      
        if (DEBUG) write(6,'("check he3 Yoni:",4F8.4)')F1,G1,G1F1,G2F1

        goto 998

C     XZ calling DSFs3. This is equivalent to above unless there are bugs in the code
      
C      print *,'check all flags:',IPOL,IA1,AsymChoice,SFchoice
C      print *,'check all flags:',pIPOL,pIA1,pAsymChoice,pSFchoice
C     print *,'check all flags:',nIPOL,nIA1,nAsymChoice,nSFchoice
        
        call DSFs3(X,QSQ,yoni_index,F1,F2,R,G1,G2,A1,A2,suppress)
        if (DEBUG) write(6,'("check he3 DSF3:",4F8.4)')F1,G1,G1/F1,G2/F1
      
        G1F1=0
        G2F1=0
        if (F1.ne.0) then
          G1F1 = G1/F1
          G2F1 = G2/F1
        endif
 998    continue
      endif
 999  continue
      return
      end
