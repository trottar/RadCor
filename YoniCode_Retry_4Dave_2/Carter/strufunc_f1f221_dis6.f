
      
********************** f1sfun ***********************************
      double precision function f1sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      anu=t/(2.*amp*aks)
      wsq=amp*amp+2*amp*anu-t
c      print *,'f1sfun:',aks,t,wsq
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
      f1sfun=(f1p*2+f1n)/3. !polrad convention is per nucleon
      end

********************** f2sfun ***********************************
      double precision function f2sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      anu=t/(2.*amp*aks)
      wsq=amp*amp+2*amp*anu-t
c      print *,'f2sfun:',aks,t,wsq
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
      f2sfun=(f2p*2+f2n)/3. !polrad convention is per nucleon
      end

 
      
CDECK  ID>, R1990.  
********************** r1990 ************************************

      double precision function r1990_orig(aks,tc)
C      double precision function r1990(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      t=tc
      if(tc.lt..35d0)t=0.35
      teta=1.+12.*t/(1.+t)*(.125**2/(aks**2+.125**2))
      zn=teta/log(t/.04)
      ra=.0672*zn+.4671/(t**4+1.8979**4)**(.25d0)
      rb=.0635*zn+.5747/t-.3534/(t**2+.09)
      rc=.0599*zn+.5088/sqrt((t-5.*(1.-aks)**5)**2+2.1081**2)
      rrr=(ra+rb+rc)/3.
      r1990=rrr
      return
      end
      
CDECK  ID>, R1990.  
********************** "r1990" NEW ************************************

c C      double precision function r1990_f1f2in(aks,t)
c       double precision function r1990(aks,t)
c       implicit real*8(a-h,o-z)
c       common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
c      .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
c       anu=t/(2.*amp*aks)
c       wsq=amp*amp+2*amp*anu-t
c       Z=1.
c       A=1.
c       call F1F2IN21(Z,A,t,wsq,f1p,f2p)
c       Z=0.
c       A=1.
c       call F1F2IN21(Z,A,t,wsq,f1n,f2n)
c       f1he=(f1p*2+f1n)/3. ! polrad convention is per-nucleon
c       f2he=(f2p*2+f2n)/3. ! polrad convention is per-nucleon
c       gam2=t/anu/anu
c       r1990=(f2he*(1+gam2)/(2.*aks))/f1he-1
c c      write(*,'("F1F221 at x=",F8.4," Q2=",F8.3," W2=",F8.2," is",F8.4)')
c c     >aks,t,wsq,r1990
c       return
c       end

CDECK  ID>, G1SFUN. XZ: modified version with F1F2_21_Nucleon and 6 GeV fits (DIS) for g1/F1
********************** g1sfun ***********************************
C      double precision function g1sfun_f1f221n_6gevfit(aks,t)
      double precision function g1sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20

C      write(*,'("calling g1sfun:",2F8.4)')aks,t
      anu=t/(2.*amp*aks)
      wsq=amp*amp+2*amp*anu-t
      
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
      
C      proton and neutron g1/F1 fits from Zheng thesis:
      asyp=(aks**(0.8126))*(1.2307-0.4128*aks)*(1+0.0303/t)
C      if (aks.gt.0.2) then
      asyn=(-0.0490-0.1618*aks+0.6979*aks*aks)*(1+0.7510/t)

      g1p=asyp*f1p
      g1n=asyn*f1n
      
      g1sfun=(-0.028*g1p*2 + 0.86*g1n)/3. ! polrad convention is per-nucleon
c      f1he=f1sfun(aks,t)
c      write(*,'("g1strf:",12F8.4)')aks,f2p,f2n,f1p,f1n,f1he,asyp,asyn,
c     >     g1p,g1n,g1sfun,g1sfun/f1he
      END FUNCTION g1sfun

     
CDECK  ID>, G2SFUN. 
********************** g2sfun ***********************************
C      double precision function g2sfun_g2ww(aks,t)
      double precision function g2sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/g2ww/tww
      external g1scal
     
c      write(*,'("calling g2ww:",2F8.4)')aks,t
      tww=t
      aksm=t/(1.08**2+t-amp2)
      if(aksm.gt.aks)then
      call dqg32(log(aks),log(aksm),g1scal,g2wwa)
      else
      g2wwa=0d0
      endif
      g2sfun=-g1sfun(aks,t)+g2wwa
      END FUNCTION g2sfun

      double precision function g1scal(aksl)
      implicit real*8(a-h,o-z)
      common/g2ww/tww
      aks=exp(aksl)
      g1scal=g1sfun(aks,tww)
      end

CDECK  ID>, DQG32.
      subroutine dqg32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end
c
