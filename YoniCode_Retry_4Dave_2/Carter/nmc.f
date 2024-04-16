C
c--------------------------------------------------------
c
c     UPDATE Jan 18, 2012
c     The definition of structure functions here (as in 
c     the quark parton model, or LO in QCD), is inconsistent
c     with the order of the PDF used, which is usually NLO 
c     or NNLO.
c
c
c   function to compute DIS cross sections using 
c   MRST or CTEQ PDFs for F1
c
c   for A>=3, EMC effect has not been included yet
c
c   inputs:
c       A: target atomic number (double precision number!)
c       e: beam energy in GeV
c       th: scattering angle in rad
c       ef: scattered electron momentum in GeV
c       IPDF: 1=CTEQ/2=MMHT
c
c   output:
c       DIS cross section in nbarn/sr/GeV, per nuclei
c--------------------------------------------------------
        subroutine get_xsec_pdf(Z,A,e,th,ef,IPDF,xsec)

        implicit none

        double precision Z,A,alpha,e,ef
        double precision mn,mp,mott,nbarn,nu
        double precision s2,t2,theta,th
        double precision xbj,qmu2,wmm
        double precision f1p,f2p,f1n,f2n
        double precision f1,f2,w1,w2,R
        double precision xsec
        double precision r1998,dr1998
        integer IPDF
        logical DEBUG/.false./

        double precision Rs,dRs,Rv,dRv,Rc,dRc,Ad,dAd,Ap,dAp,a2,da2,
     >       b2,db2,c2,dc2
        integer STAT

C        print *,'warning! inconsistency between struc. func.'
C        print *,'definition and PDF'
C        print *,'see nmc.f head notes'

        if (th.le.0) then
           theta=-th
        else
           theta=th
        endif
        mp=0.9382727
        mn=0.9395653

        alpha=7.3e-3 !  alpha = 1/137
        nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn
c
        s2=sin(abs(theta)/2.)**2

c       Calculate kinematic variables
        nu = e - ef             !energy transfer
        qmu2 = 4*e*ef*s2        !Q2
        xbj=qmu2/2/nu/mp        !Bjorken x

        wmm=mp*mp+2*mp*nu-qmu2  !invariant mass

c       Calculate Mott cross section in nbarn/(sr)
        mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn

        if (DEBUG) print *,'in xsec pdf using IPDF=',IPDF
        if (IPDF.eq.9) then
           call F1F2IN09(Z, A, qmu2, wmm, f1, f2, R)
        elseif (IPDF.eq.21) then
           call F1F2IN21(Z, A, qmu2, wmm, f1, f2)
        else
           call GETF2PDF(e,xbj,qmu2,f2p,f2n,IPDF)
           call GETF1PDF(e,xbj,qmu2,f1p,f1n,IPDF)
           f1 = Z*f1p+(A-Z)*f1n
           f2 = Z*f2p+(A-Z)*f2n
        endif
      
        w1 = f1/mp
        w2 = f2/nu
c
        t2 = tan(abs(theta)/2.)**2
c
c       DIS cross section in nbarn/sr-GeV
c
        xsec = mott * (w2 + 2 * w1 * t2)

        if (DEBUG) write(*,'("in xsec PDF: Z=",F2.0," A=",F2.0,
     >" Eb,th,Ep=",3F6.2," Mott=",F10.5," F1,F2=",2F8.4," xsec=",
     >F10.5)')Z,A,e,theta*180./3.1416,ef,mott,f1,f2,xsec
        
 103    continue
        return
        end

c--------------------------------------------------------
        subroutine get_xsec_f1f2_pdf(Z0,A0,e0,th0,ef0,IPDF,f1,f2,xsec)

        implicit none

        double precision alpha
        double precision Z0,A0,e0,th0,ef0
        double precision Z,A,e,th,ef
        double precision mn,mp,mott,nbarn,nu
        double precision s2,t2,theta
        double precision xbj,qmu2,wmm
        double precision f1p,f2p,f1n,f2n
        double precision f1,f2,df1,df2,w1,w2,dw1,dw2,R
        double precision xsec,dxsec
        double precision r1998,dr1998
        integer IPDF
        logical debug/.false./

        double precision Rs,dRs,Rv,dRv,Rc,dRc,Ad,dAd,Ap,dAp,a2,da2,
     >       b2,db2,c2,dc2
        integer STAT

        Z=Z0
        A=A0
        e=e0
        th=th0
        ef=ef0
        
C        print *,'warning! inconsistency between struc. func.'
C        print *,'definition and PDF'
C        print *,'see nmc.f head notes'

        if (th.le.0) then
           theta=-th
        else
           theta=th
        endif
        mp=0.9382727
        mn=0.9395653

        alpha=7.3e-3 !  alpha = 1/137
        nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn
c
        s2=sin(abs(theta)/2.)**2

c       Calculate kinematic variables
        nu = e - ef             !energy transfer
        qmu2 = 4*e*ef*s2        !Q2
        xbj=qmu2/2/nu/mp        !Bjorken x

        wmm=mp*mp+2*mp*nu-qmu2  !invariant mass

c       Calculate Mott cross section in nbarn/(sr)
        mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn

        if (DEBUG) print *,'in xsec_f1f2_pdf using IPDF=',IPDF

        if (IPDF.eq.9) then
           call F1F2IN09(Z, A, qmu2, wmm, f1, f2, R)
        elseif (IPDF.eq.21) then
           call F1F2IN21(Z, A, qmu2, wmm, f1, f2)
        else
           if (DEBUG) print *,'calling PDF',IPDF
           call GETF2PDF(e,xbj,qmu2,f2p,f2n,IPDF)
           call GETF1PDF(e,xbj,qmu2,f1p,f1n,IPDF)
           f1 = Z*f1p+(A-Z)*f1n
           f2 = Z*f2p+(A-Z)*f2n
        endif
      
        w1 = f1/mp
        w2 = f2/nu
c
        t2 = tan(abs(theta)/2.)**2
c
c       DIS cross section in nbarn/sr-GeV
c
        xsec = mott * (w2 + 2 * w1 * t2)

        if (DEBUG) write(*,'("in xsec PDF: Z=",F2.0," A=",F2.0,
     >" Eb,th,Ep=",3F6.2," Mott=",F10.5," F1,F2=",2F8.4," xsec=",
     >F10.5)')Z,A,e,theta*180./3.1416,ef,mott,f1,f2,xsec
        
 103    continue
        return
        end


C------------------------------------------------------------
C     subroutine to calculate F1 from PDF fits
C     IPDF=1: CT18
C          2: MMHT2014
C          
      SUBROUTINE GETF1PDF(e0,xbj0,qmu20,f1p,f1n,IPDF)
      DOUBLE PRECISION e0,xbj0,qmu20
      DOUBLE PRECISION e,xbj,qmu2,f2p,f2n
      DOUBLE PRECISION f1p,df1p,f1n,df1n,gamma2,AMP
      double precision r1998,dr1998
      INTEGER IPDF

      logical DEBUG

      DEBUG=.false.

c$$$        print *,'warning! inconsistency between struc. func.'
c$$$        print *,'definition and PDF'

      e=e0
      xbj=xbj0
      qmu2=qmu20
      
      AMP=0.93827
      
      call GETF2PDF(e,xbj,qmu2,f2p,f2n,IPDF)
      
      gamma2=4.*AMP**2.*xbj**2./qmu2

      if (DEBUG) print *,'gamma=',r1998(xbj,qmu2),gamma2,
     >     (1+r1998(xbj,qmu2))/(1+gamma2),
     >     (1+gamma2)/(1+r1998(xbj,qmu2))
        
      f1p=f2p*(1+gamma2)/(2*xbj*(1+r1998(xbj,qmu2)))
      f1n=f2n*(1+gamma2)/(2*xbj*(1+r1998(xbj,qmu2)))

C      print *,'prev F1p=',f1p/((1+gamma2)/(1+r1998(xbj,qmu2)))
C      print *,'prev F1n=',f1n/((1+gamma2)/(1+r1998(xbj,qmu2)))

      return
      end

C------------------------------------------------------------
C     subroutine to calculate F2 from PDF fits
C     IPDF=1: CT10
C          2: MSTW2008 - LO hardwired here
C          3: MMHT2014 - LO hardwired here
C          14: CT14
C
C     check R too!
C------------------------------------------------------------

      SUBROUTINE GETF2PDF(e0,xbj0,qmu20,f2p,f2n,IPDF)
      DOUBLE PRECISION e0,xbj0,qmu20
      DOUBLE PRECISION e,xbj,qmu2,f2p,f2n
      DOUBLE PRECISION uv,duv,dv,ddv,u,du,d,dd,s,ds,c,dc,
     >     ubar,dubar,dbar,ddbar,sbar,dsbar,cbar,dcbar
      DOUBLE PRECISION Ad,dAd,Ap,dAp,a2,da2,b2,db2,c2,dc2,Bd,dBd,Bp,dBp
      INTEGER IPDF
      logical DEBUG/.false./

      integer STAT

      DOUBLE PRECISION einj
      data einj/0.0/
      CHARACTER pdfname*64
        
c$$$        print *,'warning! inconsistency between struc. func.'
c$$$        print *,'definition and PDF'
C        pause

      e=e0
      xbj=xbj0
      qmu2=qmu20
      
      if (IPDF.eq.12400) then         
         pdfname="CJ15nlo"
       
      elseif (IPDF.eq.14400) then         
         pdfname="CT18NLO"
       
      elseif (IPDF.eq.25100) then
         pdfname="MMHT2014nlo68cl"

      elseif (IPDF.eq.303400) then
         pdfname="NNPDF31_nlo_as_0118"
         
      elseif (IPDF.eq.952000) then
         pdfname="JAM22_STAR_W-PDF_proton_nlo"
         pdfname="JAM22-PDF_proton_nlo"
      else
         print *,'GETF2PDF invalid iPDF:',iPDF
      endif

      STAT=0
      
      if (DEBUG) print *,'in GETF2PDF using ',IPDF,pdfname,e
         CALL CRPDF_LHAPDF_SINGLE(pdfname,0,1,0,0,0,
     >        e,xbj,qmu2,uv,duv,dv,ddv,u,du,d,dd,s,ds,c,dc,
     >        ubar,dubar,dbar,ddbar,sbar,dsbar,cbar,dcbar,
     >        Rs,dRs,Rv,dRv,Rc,dRc,Ad,dAd,Ap,dAp,
     >        Bp,dBp,Bd,dBd,STAT)
      
      if (DEBUG) print *,'u,d,s=',u,ubar,d,dbar,s,sbar
      f2p=xbj*((2./3)**2*(u+ubar)+(1./3)**2.*(d+dbar)
     >     +(1./3)**2.*(s+sbar))
      f2n=xbj*((1./3)**2*(u+ubar)+(2./3)**2.*(d+dbar)
     >     +(1./3)**2.*(s+sbar))

      return
      end


C------------------------------------------------------------
c$$$C     subroutine to calculate F1 from PDF fits
c$$$C     IPDF=1: CTEQ/2:MRST
c$$$
c$$$      SUBROUTINE GETF1PDF2(e,xbj,qmu2,f1p,df1p,f1n,df1n,IPDF)
c$$$      DOUBLE PRECISION e,xbj,qmu2,f2p,df2p,f2n,df2n
c$$$      DOUBLE PRECISION f1p,df1p,f1n,df1n,gamma2,AMP
c$$$      DOUBLE PRECISION R,dR,r1998,dr1998
c$$$      INTEGER IPDF
c$$$
c$$$      logical DEBUG
c$$$
c$$$      DEBUG=.false.
c$$$
c$$$c$$$        print *,'warning! inconsistency between struc. func.'
c$$$c$$$        print *,'definition and PDF'
c$$$C        pause
c$$$
c$$$      AMP=0.93827
c$$$
c$$$      if (DEBUG) print *,'starting GETF1PDF2',IPDF
c$$$      call GETF2PDF2(e,xbj,qmu2,f2p,df2p,f2n,df2n,R,dR,IPDF)
c$$$
c$$$      gamma2=4.*AMP**2.*xbj**2./qmu2
c$$$      
c$$$      if (DEBUG) print *,'gamma=',r1998(xbj,qmu2),gamma2,
c$$$     >     (1+r1998(xbj,qmu2))/(1+gamma2),
c$$$     >     (1+gamma2)**2./(1+r1998(xbj,qmu2))**2.
c$$$
c$$$      f1p=f2p*(1+gamma2)/(2*xbj*(1+r1998(xbj,qmu2)))
c$$$      f1n=f2n*(1+gamma2)/(2*xbj*(1+r1998(xbj,qmu2)))
c$$$
c$$$      if (DEBUG)print *,'prev F1p=',f1p/((1+gamma2)/(1+r1998(xbj,qmu2)))
c$$$      if (DEBUG)print *,'prev F1n=',f1n/((1+gamma2)/(1+r1998(xbj,qmu2)))
c$$$ 
c$$$      df1p=sqrt((df2p*(1+gamma2)/(2*xbj*(1+r1998(xbj,qmu2))))**2.
c$$$     >     +(f2p*(1+gamma2)/(2*xbj)/(1+r1998(xbj,qmu2))**2.
c$$$     >     **dr1998(xbj,qmu2))**2.)
c$$$      df1n=sqrt((df2n*(2*xbj*(1+r1998(xbj,qmu2))/(1+gamma2)))**2.
c$$$     >     +(f2n*(1+gamma2)/(2*xbj)/(1+r1998(xbj,qmu2))**2.
c$$$     >     *dr1998(xbj,qmu2))**2.)
c$$$
c$$$      return
c$$$      end
c$$$
c$$$C------------------------------------------------------------
c$$$
c$$$      SUBROUTINE GETF2PDF2(e,xbj,qmu2,f2p,df2p,f2n,df2n,R,dR,IPDF)
c$$$      dimension PDF(11)
c$$$      DOUBLE PRECISION e,xbj,qmu2,f2p,df2p,f2n,df2n,f2p0,f2n0
c$$$      DOUBLE PRECISION uv,dv,u,d,s,c,ubar,dbar
c$$$      DOUBLE PRECISION f2p1,f2n1,f2p2,f2n2
c$$$
c$$$      DOUBLE PRECISION upv,dnv,usea,dsea,str,chm,bot,glu
c$$$      DOUBLE PRECISION upv1,dnv1,usea1,dsea1,str1,chm1,bot1,glu1
c$$$      DOUBLE PRECISION upv2,dnv2,usea2,dsea2,str2,chm2,bot2,glu2
c$$$
c$$$      DOUBLE PRECISION X,Q
c$$$      DOUBLE PRECISION r1998,dr1998
c$$$      DOUBLE PRECISION R,dR,R1,R2 ! R=f2p/f2n
c$$$      INTEGER IPDF,K
c$$$
c$$$      LOGICAL DEBUG/.false./
c$$$
c$$$C
c$$$c$$$        print *,'warning! inconsistency between struc. func.'
c$$$c$$$        print *,'definition and PDF'
c$$$C        pause
c$$$
c$$$      X=xbj
c$$$      Q=sqrt(qmu2)
c$$$
c$$$      if (IPDF.eq.1) then  ! CTEQ
c$$$         call SetCtq6(1)
c$$$C            
c$$$         DO I=1,11
c$$$            PDF(I)=Ctq6Pdf(I-6,X,Q)
c$$$         ENDDO
c$$$C     nominal value:
c$$$         u=PDF(7)
c$$$         d=PDF(8)
c$$$         s=PDF(9)
c$$$         c=PDF(10)
c$$$         ubar=PDF(5)
c$$$         dbar=PDF(4)
c$$$C     
c$$$         uv=u-ubar
c$$$         dv=d-dbar
c$$$         if (DEBUG) print *,'u,d,s,c,ubar,dbar=',u,d,s,c,ubar,dbar
c$$$         f2p=xbj*((2./3)**2*(u+ubar)+(1./3)**2.*(d+dbar)+(1./3)**2.*2*s)
c$$$         f2n=xbj*((1./3)**2*(u+ubar)+(2./3)**2.*(d+dbar)+(1./3)**2.*2*s)
c$$$         f2p0=f2p
c$$$         f2n0=f2n
c$$$         R=f2p/f2n
c$$$         if (DEBUG) print *,'f2p,f2n=',f2p,f2n
c$$$C     
c$$$C     now calculating the uncertainties 
c$$$C     
c$$$         df2p=0
c$$$         df2n=0
c$$$         dR=0
c$$$C
c$$$         do K=1,20
c$$$            KK=100+2*K
c$$$            call SetCtq6(KK)
c$$$C     
c$$$C     the first set for uncertainties
c$$$C     
c$$$            DO I=1,11
c$$$               PDF(I)=Ctq6Pdf(I-6,X,Q)
c$$$            ENDDO
c$$$C     
c$$$C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
c$$$C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
c$$$C   now I =                       11,10,9, 8, 7, 6, 5, ...       1)
c$$$            u=PDF(7)
c$$$            d=PDF(8)
c$$$            s=PDF(9)
c$$$            c=PDF(10)
c$$$            ubar=PDF(5)
c$$$            dbar=PDF(4)
c$$$C     
c$$$            uv=u-ubar
c$$$            dv=d-dbar
c$$$C     
c$$$            f2p1=xbj*((2./3)**2*(u+ubar)
c$$$     >           +(1./3)**2.*(d+dbar)+(1./3)**2.*2*s)
c$$$            f2n1=xbj*((1./3)**2*(u+ubar)
c$$$     >           +(2./3)**2.*(d+dbar)+(1./3)**2.*2*s)
c$$$            R1=f2p1/f2n1
c$$$C     
c$$$C     the second set for uncertainties
c$$$C     
c$$$            KK=100+2*K-1
c$$$            call SetCtq6(KK)
c$$$C     
c$$$            DO I=1,11
c$$$               PDF(I)=Ctq6Pdf(I-6,X,Q)
c$$$            ENDDO
c$$$            u=PDF(7)
c$$$            d=PDF(8)
c$$$            s=PDF(9)
c$$$            c=PDF(10)
c$$$            ubar=PDF(5)
c$$$            dbar=PDF(4)
c$$$C     
c$$$            uv=u-ubar
c$$$            dv=d-dbar
c$$$C     
c$$$            f2p2=xbj*((2./3)**2*(u+ubar)
c$$$     >           +(1./3)**2.*(d+dbar)+(1./3)**2.*2*s)
c$$$            f2n2=xbj*((1./3)**2*(u+ubar)
c$$$     >           +(2./3)**2.*(d+dbar)+(1./3)**2.*2*s)
c$$$            R2=f2p2/f2n2
c$$$C     
c$$$C     the uncertainty is the sum of (difference between 1st set and 2nd set)^2
c$$$C
c$$$            df2p=df2p+(f2p1-f2p2)**2.
c$$$            df2n=df2n+(f2n1-f2n2)**2.
c$$$            dR=dR+(R1-R2)**2.
c$$$C
c$$$         ENDDO ! K
c$$$c     
c$$$         df2p=sqrt(df2p)/2.
c$$$         df2n=sqrt(df2n)/2.
c$$$         dR=sqrt(dR)/2.
c$$$
c$$$      else ! MRST
c$$$
c$$$         call mrst2001E(X,Q,0,upv,dnv,usea,dsea,
c$$$     >        str,chm,bot,glu)
c$$$C
c$$$         f2p=xbj*((2./3)**2*(upv+usea*2.)/X
c$$$     >        + (1./3)**2.*(dnv+dsea*2.)/X + (1./3)**2.*2*str/X)
c$$$         f2n=xbj*((1./3)**2*(upv+usea*2)/X
c$$$     >        + (2./3)**2.*(dnv+dsea*2.)/X + (1./3)**2.*2*str/X)
c$$$C
c$$$
c$$$         R=f2p/f2n
c$$$C     
c$$$C     now calculate uncertainties
c$$$C
c$$$         df2p=0
c$$$         df2n=0
c$$$         dR=0.
c$$$C
c$$$C the error is given by
c$$$C sigma(0) +- 1/2 sqrt[sum_i=1,15 {sigma(2i-1) - sigma(2i)}^2 ] C 
c$$$C
c$$$         do k=1,15
c$$$            xb_r=xb
c$$$            q_r=q
c$$$            call mrst2001E(X,Q,k*2,upv,dnv,usea,dsea,str,
c$$$     >           chm,bot,glu)
c$$$C     
c$$$            f2p1=xbj*((2./3)**2*(upv+usea*2.)/X
c$$$     >           +(1./3)**2.*(dnv+dsea*2.)/X
c$$$     >           +(1./3)**2.*2*str/X)
c$$$            f2n1=xbj*((1./3)**2*(upv+usea*2)/X
c$$$     >           +(2./3)**2.*(dnv+dsea*2.)/X
c$$$     >           +(1./3)**2.*2*str/X)
c$$$            R1=f2p1/f2n1
c$$$
c$$$            call mrst2001E(X,Q,k*2-1,upv,dnv,usea,dsea,str,
c$$$     >           chm,bot,glu)
c$$$C     
c$$$            f2p2=xbj*((2./3)**2*(upv+usea*2.)/X
c$$$     >           +(1./3)**2.*(dnv+dsea*2.)/X
c$$$     >           +(1./3)**2.*2*str/X)
c$$$            f2n2=xbj*((1./3)**2*(upv+usea*2)/X
c$$$     >           +(2./3)**2.*(dnv+dsea*2.)/X
c$$$     >           +(1./3)**2.*2*str/X)
c$$$            R2=f2p2/f2n2
c$$$C     
c$$$            df2p=df2p+(f2p1-f2p2)**2.
c$$$            df2n=df2n+(f2n1-f2n2)**2.
c$$$            dR=dR+(R1-R2)**2.
c$$$         enddo                  ! K
c$$$c     
c$$$         df2p=sqrt(df2p)/2.
c$$$         df2n=sqrt(df2n)/2.
c$$$         dR=sqrt(dR)/2.
c$$$
c$$$      endif
c$$$      return
c$$$      end
c$$$
c$$$C--------------------------------------------------------
c$$$c   subroutine to calculate PDF ratio d/u and its error
c$$$C  
c$$$      SUBROUTINE getdu(xx, du, dumax, dumin)
c$$$      DOUBLE PRECISION xx, du, dumax, dumin,x
c$$$C
c$$$C  conservative estimation from \cite{theory:duratio_2000}
c$$$C  x=xx-0.1 is because i made a mistake in the fit program
c$$$      x=xx-0.1
c$$$      du=1.0-2.2901*x+2.33738*x*x-0.8718*x*x*x
c$$$      dumax=1.0-1.789444*x+0.48548*x*x+1.048846*x*x*x
c$$$      dumin=1.0-2.1791*x+1.1766*x*x
c$$$C
c$$$C  reasonable estimation from \cite{theory:duratio}
c$$$      x=xx
c$$$C     du=1.382085-2.910495*x+2.340232*x*x-0.6266624*x*x*x;
c$$$C     dumin=1.034286-1.448902*x+0.429236*x*x;
c$$$C     dumax=*du+(*du-*dumin)/4;
c$$$C using Wally's number (email 09/15/2004)
c$$$      du=0.9850079-1.907526*x+1.320247*x*x-.2068351*x*x*x
c$$$      dumin=0.7122418-.8040541*x
c$$$      dumax=du+(du-dumin)/4
c$$$      return
c$$$      end

c--------------------------------------------------------
c   function to compute DIS cross sections using 
c   NMC95 DIS structure functions for proton and 
c   deuteron;
c
c   for A>=3, EMC effect has not been included yet
c
c   inputs:
c       A: target atomic number (double precision number!)
c       e: beam energy in GeV
c       th: scattering angle in rad
c       ef: scattered electron momentum in GeV
c
c   output:
c       DIS cross section in nbarn/sr/GeV, per nuclei
c--------------------------------------------------------
        double precision function cross_section(Z,A,e,th,ef)

        implicit none

        double precision Z,A,alpha,e,ef
        double precision mn,mp,mott,nbarn,nu
        double precision s2,t2,theta,th
        double precision xbj,qmu2,wmm2
        double precision f1psfun,f2psfun,f1nsfun,f2nsfun
        double precision f1dsfun,f2dsfun,f1hesfun,f2hesfun
        double precision f1,f2,w1,w2,tmpR

        logical DEBUG/.false./
C        logical DEBUG/.true./
c
        if (th.le.0) then
           theta=-th
        else
           theta=th
        endif
        mp=0.9382727
        mn=0.9395653

        alpha=7.3e-3 !  alpha = 1/137
        nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn
c
        s2=sin(abs(theta)/2.)**2

c       Calculate kinematic variables
        nu = e - ef             !energy transfer
        qmu2 = 4*e*ef*s2        !Q2
        xbj=qmu2/2/nu/mp        !Bjorken x

        wmm2=mp*mp+2*mp*nu-qmu2  !invariant mass

c       Calculate Mott cross section in nbarn/(sr)
        mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn

        call F1F2IN09(Z,A,qmu2, wmm2, f1, f2, tmpR)

        w1 = f1/mp
        w2 = f2/nu
c
        t2 = tan(abs(theta)/2.)**2
c
c       DIS cross section in nbarn/sr-GeV
c
        cross_section = mott * (w2 + 2 * w1 * t2)
        if (DEBUG) write(*,'("in xsec F1F2: Z=",F2.0," A=",F2.0,
     >" Eb,th,Ep=",3F6.2," Mott=",F10.5," F1,F2=",2F8.4," xsec=",
     >F10.5)')Z,A,e,theta*180./3.1416,ef,mott,f1,f2,cross_section
 103    continue
        end

c--------------------------------------------------------
c   function to compute DIS cross sections using 
c   NMC95 DIS structure functions for proton and 
c   deuteron;
c
c   for A>=3, EMC effect has not been included yet
c
c   inputs:
c       A: target atomic number (double precision number!)
c       e: beam energy in GeV
c       th: scattering angle in rad
c       ef: scattered electron momentum in GeV
c
c   output:
c       DIS cross section in nbarn/sr/GeV, per nuclei
c--------------------------------------------------------
        double precision function cross_section_nmc95(Z,A,e,th,ef)

        implicit none

        double precision Z,A,alpha,e,ef
        double precision mn,mp,mott,nbarn,nu
        double precision s2,t2,theta,th
        double precision xbj,qmu2,wmm2
        double precision f1psfun,f2psfun,f1nsfun,f2nsfun
        double precision f1dsfun,f2dsfun,f1hesfun,f2hesfun
        double precision f1,f2,w1,w2

        logical DEBUG/.false./
C        logical DEBUG/.true./
c
        DEBUG=.false.
        if (DEBUG) print *,'in NMC92:',Z,A,e,th,ef
        if (th.le.0) then
           theta=-th
        else
           theta=th
        endif
        mp=0.9382727
        mn=0.9395653

        alpha=7.3e-3 !  alpha = 1/137
        nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn
c
        s2=sin(abs(theta)/2.)**2

c       Calculate kinematic variables
        nu = e - ef             !energy transfer
        qmu2 = 4*e*ef*s2        !Q2
        xbj=qmu2/2/nu/mp        !Bjorken x

        wmm2=mp*mp+2*mp*nu-qmu2  !invariant mass

c       Calculate Mott cross section in nbarn/(sr)
        mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn
        if (A==1..and.Z==1.) then    ! for proton
           f1 = f1psfun(xbj,qmu2)
           f2 = f2psfun(xbj,qmu2)
           if (DEBUG) print *,'for proton',f1,f2,xbj,qmu2
        else if (A==1..and.Z==0.) then ! for neutron
           f1 = f1nsfun(xbj,qmu2)
           f2 = f2nsfun(xbj,qmu2)
           if (DEBUG) print *,'for neutron',f1,f2,xbj,qmu2
        else if (A==2.and.Z==1.) then  ! for deuteron
           f1 = 2.*f1dsfun(xbj,qmu2)
           f2 = 2.*f2dsfun(xbj,qmu2)
           if (DEBUG) print *,'for deuteron',f1,f2,xbj,qmu2
        else if (A==3.and.Z==2.) then     ! for 3He.
           f1 = 3.*f1hesfun(xbj,qmu2)
           f2 = 3.*f2hesfun(xbj,qmu2)
           if (DEBUG) print *,'for 3He',f1,f2,xbj,qmu2
        else                    ! for heavy nuclei
           f1 = Z*f1psfun(xbj,qmu2)+(A-Z)*f1nsfun(xbj,qmu2)
           f2 = Z*f2psfun(xbj,qmu2)+(A-Z)*f2nsfun(xbj,qmu2)
        endif

        w1 = f1/mp
        w2 = f2/nu
c
        t2 = tan(abs(theta)/2.)**2
c
c       DIS cross section in nbarn/sr-GeV
c
        cross_section_nmc95 = mott * (w2 + 2 * w1 * t2)
        if (DEBUG) write(*,'("in xsec NMC: Z=",F2.0," A=",F2.0,
     >" Eb,th,Ep=",3F6.2," Mott=",F10.5," F1,F2=",2F8.4," xsec=",
     >F10.5)')Z,A,e,theta*180./3.1416,ef,mott,f1,f2,cross_section_nmc95
 103    continue
        end

********************** f1psfun ***********************************
* nucleon averaged structure functions

      double precision function f1psfun(aks,t)
      implicit double precision(a-h,o-z)
      f2p=d95f2h8(t,aks)
      amp=0.9382727
      dnu=t/(2*amp*aks)
      f1psfun=f2p*(1+t/dnu**2)/(2.*aks*(1.+r1998(aks,t)))
      end


********************** f2psfun ***********************************
* nucleon averaged structure functions

      double precision function f2psfun(aks,t)
      implicit double precision(a-h,o-z)
      f2psfun=d95f2h8(t,aks)
      end


********************** f1nsfun ***********************************
      double precision function f1nsfun(aks,t)
      implicit double precision(a-h,o-z)
      df1d=f1dsfun(aks,t)
      df1p=f1psfun(aks,t)
      f1nsfun=2*df1d-df1p
      end



********************** f2nsfun ***********************************
      double precision function f2nsfun(aks,t)
      implicit double precision(a-h,o-z)
      df2d=f2dsfun(aks,t)
      df2p=f2psfun(aks,t)
      f2nsfun=2*df2d-df2p
      end


********************** f1dsfun ***********************************
* note: per nucleon definition
*
      double precision function f1dsfun(aks,t)
      implicit double precision(a-h,o-z)
      f2=d95f2d8(t,aks)
      amp=0.9382727
      dnu=t/(2*amp*aks)
      f1dsfun=f2*(1+t/dnu**2)/(2.*aks*(1.+r1998(aks,t)))
C      print *,'f1dsfun:',f1dsfun,f2,r1998(aks,t)
      end

********************** f2dsfun ***********************************
* note: per nucleon definition
*
      double precision function f2dsfun(aks,t)
      implicit double precision(a-h,o-z)
      f2dsfun=d95f2d8(t,aks)
      end


********************** f1hesfun ***********************************
      double precision function f1hesfun(aks,t)
      implicit double precision(a-h,o-z)
      f2=f2hesfun(aks,t)
      amp=0.9382727
      dnu=t/(2*amp*aks)
      f1hesfun=f2*(1+t/dnu**2)/(2.*aks*(1.+r1998(aks,t)))
c      print *,'nmc95:',aks,t,f1hesfun
      end

********************** f2hesfun ***********************************
      double precision function f2hesfun(aks,t)
      implicit double precision(a-h,o-z)
      f2p=d95f2h8(t,aks)
      f2d=d95f2d8(t,aks)
      f2hesfun=(f2p+2d0*f2d)/3d0
c      print *,'nmc95:',f2hesfun
      end

********************** d95f2d8  ***********************************
      double precision function d95f2d8(dq2,dx)
*:=====================================================================:
*:                                                                     :
*:      author:    x.zheng        last update: 11.25.2001              :
*:                                 tested: xxx                         :
*:                                                                     :
*:      arguments: dq2,dx: double prec. input xbj,q2                   :
*:                 d95f2h8* double prec f2  output                     :
*:                                                                     :
*:      called by: mkf2                                                :
*:                                                                     :
*:      action:    calculate f2 structure function of the deuteron     :
*:                 nmc fit of dis-region with 15 parameters fit        :
*:                 kinematics range: 0.5<q2<75 GeV^2, 0.006<x<0.9      :
*:                                                                     :
*:                 parametrized with a1-7,b1-4,c1-4  as                :
*:                                                                     :
*:                 f2_dis(x,q2) ~prop.                                 :
*:                   A(x)*(ln(q2/dl2)/ln(q20/dl2))**(B(x))*(1+C(x)/q2) :
*:                   with x = (q2+m_a)/(2m*nu + m_b**2)                :
*:                        dl2 = (0.250 GeV)^2, (so-called lambda^2)    :
*:                        q20 = 20 GeV^2                               :
*:                        A(x) = x**a1*(1-x)**a2*(a3+a4*(1-x)          :
*:                               +a5*(1-x)**2+a6*(1-x)**3+a7*(1-x)**4) :
*:                        B(x) = b1+b2*x+b3/(x+b4)                     :
*:                        C(x) = c1*x+c2*x**2+c3*x**3+c4*x**4          :
*:                 reference:                                          :
*:                 the new muon collaboration                          :
*:                 Phys.Lett.B364(1995)107~115                         :
*:                 hep-ph/9509406                                      :
*:                                                                     :
*:                 resonance contribution is calculated in the same    :
*:                 way as df2d8(dq2,dx)                                :
*:                                                                     :
*:       comments: The proton and deuteron structure functions F2p and :
*:                 F2d were measured in the kinematic range 0.006<x<0.6:
*:                 and 0.5<q2<75 GeV^2, by inclusive deep inelastic    :
*:                 muon scattering at 90,120,300 and 280 GeV.  The     :
*:                 measurements are in good agreement with earlier high:
*:                 precision results.  The present and earlier results :
*:                 together have been parametrised to give descriptions:
*:                 of the proton and deuteron structure functions F2   :
*:                 and their uncertainties over the range 0.006<x<0.9  :
*:=====================================================================:
c
      implicit double precision (d)
c
c
c *** a1,..a7,b1..b4,c1..c4 = 15 param of nmc, slac, bcdms (95)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c            values: daw=1.512(gev2), dbw=0.351(gev2)
c            ref:  bodek et al., p.r.d20(1979)1427.
c            see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.2**2 = 0.04 (gev2)
c *** q0**2 = 2 gev2 ... (2+0.351)/0.04 = 58.771
c *** fit by y.m.(25-nov-88 19h43m14s)
c
      data a1,a2,a3,a4,a5,a6
     :     ,a7,b1,b2,b3,b4
     :     ,c1,c2,c3,c4
     :     ,d9,d10
     :     ,daw,dbw
c
c     f2 from nmc phys.lett.b364(1995)107
     :     /-0.04858,2.863,0.8367,-2.532,9.145,-12.504
     :     ,5.473,-0.008,-2.227,0.0551,0.0570
     :     ,-1.509,8.553,-31.20,39.98
c     resonance-region
     :     ,.89456,.16452
     :     ,1.512,.351 /
c
c
      d95f2d8=1.d-30
      dl2 = 0.25**2
      dq20 = 20
      d95f2d8 = (dx**a1*(1-dx)**a2*(a3+a4*(1-dx)
     :      +a5*(1-dx)**2+a6*(1-dx)**3+a7*(1-dx)**4))
     :      *(log(dq2/dl2)/log(dq20/dl2))
     :      **(b1+b2*dx+b3/(dx+b4))
     :      *(1+(c1*dx+c2*dx**2+c3*dx**3+c4*dx**4)/dq2)
C
      if (d95f2d8.gt.0d0) return
      d95f2d8 = 1.d-30
      return
      end
      
      double precision function d95f2d8_lo(dq2,dx)
      implicit double precision (d)

      data a1,a2,a3,a4,a5,a6
     :     ,a7,b1,b2,b3,b4
     :     ,c1,c2,c3,c4
C     :     ,d9,d10
C     :     ,daw,dbw
     :     /-0.04715,2.814,0.7286,-2.151,8.662,-12.258, 
     :     5.452, -0.048, -2.114, 0.0672, 0.0677, 
     :     -1.517, 9.515, -34.94, 44.42/

      d95f2d8_lo=1.d-30
      dl2 = 0.25**2
      dq20 = 20
      d95f2d8_lo = (dx**a1*(1-dx)**a2*(a3+a4*(1-dx)
     :      +a5*(1-dx)**2+a6*(1-dx)**3+a7*(1-dx)**4))
     :      *(log(dq2/dl2)/log(dq20/dl2))
     :      **(b1+b2*dx+b3/(dx+b4))
     :      *(1+(c1*dx+c2*dx**2+c3*dx**3+c4*dx**4)/dq2)
C
      if (d95f2d8_lo.gt.0d0) return
      d95f2d8_lo = 1.d-30
      return
      end
      
      double precision function d95f2d8_hi(dq2,dx)
      implicit double precision (d)

      data a1,a2,a3,a4,a5,a6
     :     ,a7,b1,b2,b3,b4
     :     ,c1,c2,c3,c4
C     :     ,d9,d10
C     :     ,daw,dbw
     :     /-0.02732, 2.676, 0.3966, -0.608, 4.946, -7.994, 
     :     3.686, 0.141, -2.464, 0.0299, 0.0396, 
     :     -2.128, 14.378, -47.76, 53.63/
      
      d95f2d8_hi=1.d-30
      dl2 = 0.25**2
      dq20 = 20
      d95f2d8_hi = (dx**a1*(1-dx)**a2*(a3+a4*(1-dx)
     :      +a5*(1-dx)**2+a6*(1-dx)**3+a7*(1-dx)**4))
     :      *(log(dq2/dl2)/log(dq20/dl2))
     :      **(b1+b2*dx+b3/(dx+b4))
     :      *(1+(c1*dx+c2*dx**2+c3*dx**3+c4*dx**4)/dq2)
C
      if (d95f2d8_hi.gt.0d0) return
      d95f2d8_hi = 1.d-30
      return
      end


********************** d95f2h8  ***********************************
      double precision function d95f2h8(dq2,dx)
*:=====================================================================:
*:                                                                     :
*:      author:    x.zheng        last update: 11.25.2001              :
*:                                tested: xxx                          :
*:                                                                     :
*:      arguments: dq2,dx: double prec. input xbj,q2                   :
*:                 d95f2h8* double prec f2  output                   :
*:                                                                     :
*:      called by: mkf2                                                :
*:                                                                     :
*:      action:    calculate f2 structure function of the proton       :
*:                 nmc fit of dis-region with 15 parameters fit        :
*:                 kinematics range: 0.5<q2<75 GeV^2, 0.006<x<0.9      :
*:                                                                     :
*:                 parametrized with a1-7,b1-4,c1-4  as                :
*:                                                                     :
*:                 f2_dis(x,q2) ~prop.                                 :
*:                   A(x)*(ln(q2/dl2)/ln(q20/dl2))**(B(x))*(1+C(x)/q2) :
*:                   with x = (q2+m_a)/(2m*nu + m_b**2)                :
*:                        dl2 = (0.250 GeV)^2, (so-called lambda^2)    :
*:                        q20 = 20 GeV^2                               :
*:                        A(x) = x**a1*(1-x)**a2*(a3+a4*(1-x)          :
*:                               +a5*(1-x)**2+a6*(1-x)**3+a7*(1-x)**4) :
*:                        B(x) = b1+b2*x+b3/(x+b4)                     :
*:                        C(x) = c1*x+c2*x**2+c3*x**3+c4*x**4          :
*:                 reference:                                          :
*:                 the new muon collaboration                          :
*:                 Phys.Lett.B364(1995)107~115                         :
*:                 hep-ph/9509406                                      :
*:                                                                     :
*:                 resonance contribution is calculated in the same    :
*:                 way as df2d8(dq2,dx)                                :
*:                                                                     :
*:       comments: The proton and deuteron structure functions F2p and :
*:                 F2d were measured in the kinematic range 0.006<x<0.6:
*:                 and 0.5<q2<75 GeV^2, by inclusive deep inelastic    :
*:                 muon scattering at 90,120,300 and 280 GeV.  The     :
*:                 measurements are in good agreement with earlier high:
*:                 precision results.  The present and earlier results :
*:                 together have been parametrised to give descriptions:
*:                 of the proton and deuteron structure functions F2   :
*:                 and their uncertainties over the range 0.006<x<0.9  :
*:=====================================================================:
c
      implicit double precision (d)
*
c *** a1,..a7,b1..b4,c1..c4 = 15 param of nmc, slac, bcdms (95)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c            values: daw=1.512(gev2), dbw=0.351(gev2)
c            ref:  bodek et al., p.r.d20(1979)1427.
c            see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.25**2 = 0.0625 (gev2)
c *** q0**2 = 20 gev2
*
      data
     .amm/2.7928456d0/,amn/-1.913148d0/,pi/3.1415926d0/
     .,alfa/.729735d-2/,amh/.938272d0/,ampi/.104/
c
      data a1,a2,a3,a4,a5,a6
     :     ,a7,b1,b2,b3,b4
     :     ,c1,c2,c3,c4
     :     ,d9,d10,d11,d12,d13,d14
     :     ,d15,d16
     :     ,daw,dbw
c
c     f2 from nmc phys.lett.b364(1995)107
     :     /-.02778,2.926,1.0362,-1.840,8.123,-13.074
     :     ,6.215,0.285,-2.694,0.0188,0.0274
     :     ,-1.413,9.366,-37.79,47.10
c     resonance-region:
     :     ,.1179, .044735, .038445, .27921, 8.8228d-5, 6.2099d-5
     :     ,1.421,1.2582
     :     ,1.642, .376/
c
      d95f2h8=1.d-30
      dl2 = 0.25**2
      dq20 = 20.
c
      d95f2h8=1.d-30
      dl2 = 0.25**2
      dq20 = 20.
      d95f2h8 = (dx**a1*(1-dx)**a2*(a3+a4*(1-dx)
     :      +a5*(1-dx)**2+a6*(1-dx)**3+a7*(1-dx)**4))
     :      *(log(dq2/dl2)/log(dq20/dl2))
     :      **(b1+b2*dx+b3/(dx+b4))
     :      *(1+(c1*dx+c2*dx**2+c3*dx**3+c4*dx**4)/dq2)
c
      dres = 0.d0
c
c *** (total) = (qcd part) + (resonance region)
c
      d95f2h8 = d95f2h8 + dres
c
      if(d95f2h8 .gt. 0.d0) return
      d95f2h8=1.d-30
c
      return
      end
c




c--------------------------------------------------------
c   function to compute RES cross sections using F1F2
c     from P. Bosted and E. Christy
c
c   inputs:
c       A: target atomic number (double precision number!)
c       e: beam energy in GeV
c       th: scattering angle in rad
c       ef: scattered electron momentum in GeV
c
c   output:
c       RES cross section in nbarn/sr/GeV, per nuclei
c--------------------------------------------------------
        double precision function cross_section_res(Z,A,e,th,ef)

        implicit none

        double precision Z,A,alpha,e,ef
        double precision mn,mp,mott,nbarn,nu
        double precision s2,t2,theta,th
        double precision xbj,qmu2,wmm
        double precision w1,w2

        real*8 resZ, resA, QSQ, Wsq, F1, F2, Rc

        logical RESDEBUG/.false./
c
        if (th.le.0) then
           theta=-th
        else
           theta=th
        endif
        mp=0.9382727
        mn=0.9395653

        alpha=7.3e-3 !  alpha = 1/137
        nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn
c
        s2=sin(abs(theta)/2.)**2

c       Calculate kinematic variables
        nu = e - ef             !energy transfer
        qmu2 = 4*e*ef*s2        !Q2
        xbj=qmu2/2/nu/mp        !Bjorken x

        wmm=mp*mp+2*mp*nu-qmu2  !invariant mass

c       Calculate Mott cross section in nbarn/(sr)
        mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn

        resZ=Z
        resA=A
        QSQ=qmu2
        Wsq=wmm

        if (wmm.lt.1.12**2.) then
           call F1F2QE09(resZ, resA, QSQ, Wsq, F1, F2, Rc)
        if (RESDEBUG) print *,'in QE: Z=',resZ,' A=',resA,
     >          ' f1=',f1,' f2=',f2
        if (RESDEBUG) print *,'in QE: w,q2,x=',sqrt(wmm),qmu2,xbj
        else
           call F1F2IN09(resZ, resA, QSQ, Wsq, F1, F2, Rc)
        if (RESDEBUG) print *,'in RES: Z=',resZ,' A=',resA,
     >          ' f1=',f1,' f2=',f2
        if (RESDEBUG) print *,'in RES: w,q2,x=',sqrt(wmm),qmu2,xbj
        endif
           
        w1 = f1/mp
        w2 = f2/nu
c
        t2 = tan(abs(theta)/2.)**2
c
c       DIS cross section in nbarn/sr-GeV
c
        cross_section_res = mott * (w2 + 2 * w1 * t2)
        if (RESDEBUG) 
     >       print *,'mott=',mott,' f1=',f1,' f2=',f2,
     >       ' w1=',w1,' w2=',w2,' xsec=', cross_section_res
 103    continue
        end

