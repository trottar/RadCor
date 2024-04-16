
      
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
      f1sfun=(f1p*2+f1n)/3. ! polrad convention is per-nucleon
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
      f2sfun=(f2p*2+f2n)/3. ! polrad convention is per-nucleon
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

C      double precision function r1990_f1f2in(aks,t)
      double precision function r1990(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      anu=t/(2.*amp*aks)
      wsq=amp*amp+2*amp*anu-t
      Z=1.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1p,f2p)
      Z=0.
      A=1.
      call F1F2IN21(Z,A,t,wsq,f1n,f2n)
      f1he=(f1p*2+f1n)/3. ! polrad convention is per-nucleon
      f2he=(f2p*2+f2n)/3. ! polrad convention is per-nucleon
      gam2=t/anu/anu
      r1990=(f2he*(1+gam2)/(2.*aks))/f1he-1
      write(*,'("F1F221 at x=",F8.4," Q2=",F8.3," W2=",F8.2," is",F8.4)')
     >aks,t,wsq,r1990
      return
      end

      
CDECK  ID>, G1SFUN. XZ: modified version with F1F2_21N (nucleon) and CLAS fit for g1/F1
********************** g1sfun ***********************************
C      double precision function g1sfun_clas(aks,t)
      double precision function g1sfun(aks,t)
      implicit real*8(a-h,o-z)
      integer iA,iZ
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
C      common/inkin/bmom,tmom

C      write(*,'("calling g1sfun:",2F8.4)')aks,t
      anu=t/(2.*amp*aks)
      wsq=amp*amp+2*amp*anu-t

      goto 999
      
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
      
      g1sfun=(-0.028*g1p*2 + 0.86*g1n)/3. ! )/3. ! polrad convention is per-nucleon
      f1he=f1sfun(aks,t)
C     --- DIS approach ends here ---
 999  continue
C     --- CLAS model begins ---
      iA=3
      iZ=2
      call sig_pol_dis_calc(iA, iZ, aks, t, g1f1, g2f1)
      
C      write(*,'("g1strf:",F8.4," F2p,n=",2F8.4," F1p,n=",2F8.4," F1he3=",
C     >F8.4," A1p,n=",2F8.4," G1p,n=",2F8.4," G1he3=",F8.4," g1f1he3 DIS=",F8.4,
C     >" CLAS g1f1he3=",F8.4)')aks,f2p,f2n,f1p,f1n,f1he,asyp,asyn,
C     >     g1p,g1n,g1sfun,g1sfun/f1he,g1f1
      
      f1he=f1sfun(aks,t)
      g1sfun=g1f1*f1he

      end
C     
      
CDECK  ID>, G2SFUN. 
********************** g2sfun ***********************************
C      double precision function g2sfun_clas(aks,t)
      double precision function g2sfun(aks,t)
      implicit real*8(a-h,o-z)
      integer iA,iZ
C      write(*,'("calling g2sfun:",2F8.4)')aks,t
C     --- CLAS model begins ---
      iA=3
      iZ=2
      call sig_pol_dis_calc(iA, iZ, aks, t, g1f1, g2f1)
      
C      write(*,'("g1strf:",F8.4," F2p,n=",2F8.4," F1p,n=",2F8.4," F1he3=",
C     >F8.4," A1p,n=",2F8.4," G1p,n=",2F8.4," G1he3=",F8.4," g1f1he3 DIS=",F8.4,
C     >" CLAS g1f1he3=",F8.4)')aks,f2p,f2n,f1p,f1n,f1he,asyp,asyn,
C     >     g1p,g1n,g1sfun,g1sfun/f1he,g1f1
      
      f1he=f1sfun(aks,t)
      g2sfun=g2f1*f1he

      end
