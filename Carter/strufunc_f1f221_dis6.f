c include "F1F2IN21.f"
      

********************** f1sfun ***********************************
      double precision function f1sfun(aks,t)
c      &BIND(C, name = "f1sfun_")
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      anu=t/ap/aks
      wsq=ap*ap+2*ap*anu-t
c      PRINT *, "Test 1"
c      PRINT *, aks, t, wsq, ap
c      print *,'f1sfun:',aks,t,wsq
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
c      PRINT *, f1p
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
c      PRINT *, f1n
      f1sfun=f1p*2+f1n
      end

********************** f2sfun ***********************************
      double precision function f2sfun(aks,t)
c      &BIND(C, name = "f2sfun_")
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      anu=t/ap/aks
      wsq=ap*ap+2*ap*anu-t
c      print *,'f2sfun:',aks,t,wsq
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
      f2sfun=f2p*2+f2n
      end

      
      
CDECK  ID>, R1990.  
********************** r1990 NEW ************************************

      double precision function r1990_f1f2in(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      anu=t/ap/aks
      wsq=ap*ap+2*ap*anu-t
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
      f1he=f1p*2+f1n
      f2he=f2p*2+f2n
      gam2=t/anu/anu
      r1990=(f2he*(1+gam2)/(2.*aks))/f1he-1
      print *,'r1990:',aks,t,wsq,r1990
      return
      end

CDECK  ID>, G1SFUN. XZ: modified version with F1F2_21_Nucleon and 6 GeV fits (DIS) for g1/F1
********************** g1sfun ***********************************
C      double precision function g1sfun_f1f221n_6gevfit(aks,t)
      double precision function g1sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20

C      write(*,'("calling g1sfun:",2F8.4)')aks,t
      anu=t/ap/aks
      wsq=ap*ap+2*ap*anu-t
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
      
      g1sfun=(-0.028*g1p*2 + 0.86*g1n)/3. ! Since g1 is calculated as A1*F1, it should be per-nucleon
      f1he=f1sfun(aks,t)
c      write(*,'("g1strf:",12F8.4)')aks,f2p,f2n,f1p,f1n,f1he,asyp,asyn,
c     >     g1p,g1n,g1sfun,g1sfun/f1he

      end

     
CDECK  ID>, G2SFUN. 
********************** g2sfun ***********************************
c C      double precision function g2sfun_g2ww(aks,t)
c       double precision function g2sfun(aks,t)
c       implicit real*8(a-h,o-z)
c       common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
c      .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
c       common/g2ww/tww
c       external g1scal
      
c C      write(*,'("calling g2ww:",2F8.4)')aks,t
c       tww=t
c       aksm=t/(1.08**2+t-amp2)
c       if(aksm.gt.aks)then
c        call dqg32(log(aks),log(aksm),g1scal,g2wwa)
c        else
c        g2wwa=0d0
c        endif
c        g2sfun=-g1sfun(aks,t)+g2wwa
c       end
c       double precision function g1scal(aksl)
c       implicit real*8(a-h,o-z)
c       common/g2ww/tww
c       aks=exp(aksl)
c       g1scal=g1sfun(aks,tww)
c       end
