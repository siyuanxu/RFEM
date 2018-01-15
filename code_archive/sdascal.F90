
!****************************************************************************
!
!  PROGRAM: DASAC
!  (DEFORMATION AND STRESS ANALYSIS CONSIDERING CONSOLIDATION.)
!
!****************************************************************************
MODULE MAXDEF

	IMPLICIT NONE
	INTEGER,PARAMETER::NE_MAX=150000
	INTEGER,PARAMETER::NX_MAX=150000
	INTEGER,PARAMETER::NP_MAX=1000
	INTEGER,PARAMETER::NM_MAX=1000
	INTEGER,PARAMETER::NF_MAX=25000
	INTEGER,PARAMETER::NQ_MAX=20000
    INTEGER,PARAMETER::NF1_MAX=50000
	INTEGER,PARAMETER::NP1_MAX=80000
	INTEGER,PARAMETER::NTEM=70000000
END MODULE MAXDEF


MODULE CONDEF
	IMPLICIT NONE
	INTEGER,PARAMETER::NDM=3
	INTEGER,PARAMETER::NDS=8
	INTEGER,PARAMETER::NFS=4
	INTEGER,PARAMETER::NFR=4
	INTEGER,PARAMETER::NGS=2
	INTEGER,PARAMETER::NGF=7
	REAL(KIND=8),PARAMETER::GATA=0.6670
	REAL(KIND=8),PARAMETER::GMAW=9.810
	REAL(KIND=8),PARAMETER::GMAG=9.810
	REAL(KIND=8),PARAMETER::PAR=101.325
	REAL(KIND=8),PARAMETER::PAI=3.142
	REAL(KIND=8),PARAMETER::RAD=57.296
END MODULE CONDEF


MODULE COMM
	USE MAXDEF
	USE CONDEF
	IMPLICIT NONE
!---------------------------------------------------------
	INTEGER,DIMENSION(10)::ICON
	real(kind=8)::hsur
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
	INTEGER::NP,NM,NE,NX,NF,NQ,NF1,NP1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER,DIMENSION(NP_MAX,10)::LQ
	REAL(KIND=8),DIMENSION(NP_MAX,2)::HWU,HWD
	REAL(KIND=8),DIMENSION(NP_MAX)::conh
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(KIND=8),DIMENSION(NP1_MAX,2)::PSJ
	REAL(KIND=8),DIMENSION(NE_MAX,3)::SEEP
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER,DIMENSION(NM_MAX,7)::KGSS
	REAL(KIND=8),DIMENSION(NM_MAX,15)::CS,CK,CW,CR
	INTEGER,DIMENSION(NE_MAX,NDS)::NOD,nod2
	INTEGER,DIMENSION(NE_MAX)::MAT,JAD,JRD,SRD
	REAL(KIND=8),DIMENSION(NX_MAX,NDM)::COR
	INTEGER,DIMENSION(NX_MAX,NFR)::ID,ig
	REAL(KIND=8),DIMENSION(NX_MAX,NFR)::GG
	REAL(KIND=8),DIMENSION(NE_MAX,12)::STR
	INTEGER,DIMENSION(NF_MAX,NFS+3)::KSF
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTEGER,DIMENSION(NF1_MAX,NFS+3)::KSP
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER,DIMENSION(NQ_MAX,NFS+2)::KSQ
	REAL(KIND=8),DIMENSION(NQ_MAX)::DSQ
!---------------------------------------------------------
	REAL(KIND=8),DIMENSION(NGS)::WGH,WXYZ
	REAL(KIND=8),DIMENSION(NGF)::VGH,VXYZ
!---------------------------------------------------------
	INTEGER,DIMENSION(NX_MAX)::JR,JST
	INTEGER,DIMENSION(NE_MAX)::LUS,IST,jplane
	INTEGER,DIMENSION(NX_MAX*NFR)::MA,MB
!	INTEGER,DIMENSION(NTEM)::NDEX
	INTEGER,DIMENSION(NX_MAX,NFR)::JRR
	INTEGER,DIMENSION(NP_MAX)::NDOF,NVOL
!---------------------------------------------------------
	INTEGER,DIMENSION(NE_MAX)::LDG
	REAL(KIND=8),DIMENSION(NE_MAX)::PWS,thata
	REAL(KIND=8),DIMENSION(NE_MAX,2)::YLD,PLS,YLDM
	REAL(KIND=8),DIMENSION(NE_MAX,4)::EST
	REAL(KIND=8),DIMENSION(NE_MAX,6,6)::DMT
	REAL(KIND=8),DIMENSION(NE_MAX,NDS)::RLD
	REAL(KIND=8),DIMENSION(NE_MAX,12)::STS,STN
	REAL(KIND=8),DIMENSION(NE_MAX,6)::DSTS,DSTN,rhstn
	REAL(KIND=8),DIMENSION(NE_MAX,NDM,NDM)::VCS,VCN
!---------------------------------------------------------
	REAL(KIND=8),DIMENSION(NX_MAX)::TWH,DTWH
	REAL(KIND=8),DIMENSION(NX_MAX,NDM)::DIS,DDIS
	REAL(KIND=8),DIMENSION(NX_MAX,NDM)::NOFCE
	REAL(KIND=8),DIMENSION(NX_MAX,NFR)::NODIS
!---------------------------------------------------------
	REAL(KIND=8),DIMENSION(NE_MAX,2)::RHS
	REAL(KIND=8),DIMENSION(NE_MAX,NDM)::EQK
	real(kind=8)::totaltime
	integer::ntlb,ntep
	real(kind=8),dimension(10000,2)::shuiwei
!---------------------------------------------------------
END MODULE COMM


PROGRAM DASAC
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,i,j,mt,ie
	CHARACTER(LEN=30)::MAINFILE
	CALL SEARCHFILE(MAINFILE)
	CALL FILEOPEN(MAINFILE)
	CALL DATAINP
	CALL DATAOUT
	CALL PREPROC
	CALL INITCOND
    totaltime=0.0 
	ntep=0
	call svrst(0)
	DO IP=1,NP
      if(ip.eq.10)then
        do ie=1,ne
		  !if(mat(ie).eq.6.or.mat(ie).eq.16)mat(ie)=18             !平台建设完毕后进行防渗墙浇筑，材料变换
		  if(mat(ie).eq.6.or.mat(ie).eq.16)MAT(IE)=17+SRD(IE)
		end do
	  end if
		CALL MQQ(IP)
		CALL MPP(IP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF(LQ(IP,1)<0)THEN
			CALL PRTIME
			PRINT*,'STEP=',IP
			PRINT*,'...EXCAVATING CONDITION...'
			CALL CALCEX(IP)
			IF(LQ(IP,8)/=0) CALL CALCRH(IP)
			totaltime=totaltime+lq(ip,3)
		ELSE IF(LQ(IP,1)==1)THEN
			CALL PRTIME
			PRINT*,'STEP=',IP
			PRINT*,'...DAM FILLING CONDITION...'
			CALL CALCFL(IP)
			IF(LQ(IP,8)/=0) CALL CALCRH(IP)
			totaltime=totaltime+lq(ip,3)
		ELSE IF(LQ(IP,1)==2)THEN
			CALL PRTIME
			PRINT*,'STEP=',IP
			PRINT*,'...WATER FILLING CONDITION...'
			CALL CALCWF(IP)
			IF(LQ(IP,8)/=0) CALL CALCRH(IP)
			totaltime=totaltime+lq(ip,3)
		ELSE IF(LQ(IP,1)==3)THEN
			CALL PRTIME
			PRINT*,'STEP=',IP
			PRINT*,'...EARTHQUAKING CONDITION...'
			CALL CALCQK(IP)
		ELSE IF(LQ(IP,1)==4)THEN
			CALL PRTIME
			PRINT*,'STEP=',IP
			PRINT*,'...RHELOGICAL DEFORMATION CONDITION...'
			CALL CALCRH(IP)
		ELSE IF(LQ(IP,1)==5)THEN
			CALL PRTIME
			PRINT*,'STEP=',IP
			PRINT*,'...CONSOLIDATION CONDITION...'
			CALL CALCSD(IP)
		END IF
		IF(LQ(IP,10)/=0)THEN
			CALL SVRST(IP)
		END IF
	END DO
	CALL FILECLOSE
	PRINT*,'...CALCULATION SUCCESSFULLY FINISHED...'
	CALL PRTIME
END PROGRAM DASAC


SUBROUTINE CALCEX(IP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::DT,RT
	INTEGER::IP,IE,IX,J,K,L,KL,ME,LE,NS,ND,NJ,NK,IM,IL,NN,NV,NSK,LC
	REAL(KIND=8),DIMENSION(NDS)::QV
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	REAL(KIND=8),ALLOCATABLE::SK(:),RD(:)
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	INTEGER,ALLOCATABLE::NEX(:)
	ME=LQ(IP,2)
	ND=LQ(IP,3)
	NJ=LQ(IP,4)
	NK=LQ(IP,5)
	NS=LQ(IP,9)
	IM=ICON(10)
	IF(NS<1) NS=1
	DT=REAL(ND)/REAL(NS)
	RT=1.0/REAL(NS)
	NSK=1
	IF(ICON(1)/=0) NSK=2
	DO IX=1,NX
		DO J=1,NDM
			NOFCE(IX,J)=0.0
		END DO
		DO K=1,NFR
			NODIS(IX,K)=0.0
		END DO
	END DO
	IF(NJ/=0) CALL NOFCEIN(NJ)
	IF(NK/=0) CALL NODISIN(NK)
	DO IX=1,NX
 		DO J=1,NDM
			NOFCE(IX,J)=NOFCE(IX,J)*RT
		END DO
		DO K=1,NFR
			NODIS(IX,K)=NODIS(IX,K)*RT
		END DO
	END DO
	DO IE=1,ME
		LE=LUS(IE)
		IF(IST(LE)/=-1) CYCLE
		FORALL(L=1:NDS*NDM) FF(L)=0.0
		CALL UNBALFRC(LE,FF)
		KL=0
		DO J=1,NDS
			K=NOD(LE,J)
			DO L=1,NDM
				KL=KL+1
				NOFCE(K,L)=NOFCE(K,L)+FF(KL)*RT
			END DO
		END DO
		DO J=1,12
			STN(LE,J)=0.0
			STS(LE,J)=0.0
		END DO
		PWS(LE)=0.0
	END DO
!	ALLOCATE(NEX(NTEM))
	IF(IM==0)THEN
		CALL MRR(IP)
	ELSE IF(IM==1)THEN
		CALL MSS(IP)
	ELSE IF(IM==2)THEN
	    CALL MMM(IP)
	END IF
	NN=NDOF(IP)
	NV=NVOL(IP)
	ALLOCATE(SK(NV),RD(NN),NEX(NV))
	CALL MMKK(IP,NN,NV,NEX)
	DO LC=1,NS
		PRINT('(A10,2X,A8,2X,I2)'),'*****','SUBSTEP:',LC
		DO IE=1,NE
			DO J=1,NDS
				RLD(IE,J)=0.0
			END DO
		END DO
		IF(ICON(2)/=0)THEN
			DO J=1,NQ
				IE=KSQ(J,1)
				IF(IST(IE)==0) CYCLE
				FORALL(K=1:NDS) QV(K)=0.0
				CALL SRFQ(IP,J,QV)
				DO K=1,NDS
					RLD(IE,K)=RLD(IE,K)+QV(K)*DT
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				FORALL(L=1:NDS) QV(L)=0.0
				CALL VOLQ(LE,QV)
				DO L=1,NDS
					RLD(LE,L)=RLD(LE,L)+QV(L)*DT
				END DO
			END DO
		END IF
		DO IL=1,NSK
			DO J=1,NN
				RD(J)=0.0
			END DO
			DO K=1,NV
				SK(K)=0.0
			END DO
			DO IX=1,NX
				DO J=1,NDM
					DDIS(IX,J)=0.0
				END DO
				DTWH(IX)=0.0
			END DO
			DO IE=1,NE
				DO J=1,6
					DSTN(IE,J)=0.0
					DSTS(IE,J)=0.0
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				CALL EVSNCALC(LE)
				CALL STFNS(LE,STF,DT)
				FORALL(K=1:NDS*NFR) GQ(K)=0.0
				CALL BNDAPP(LE,STF,GQ)
				DO K=1,NDS*NFR
					GQ(K)=-GQ(K)
				END DO
				DO K=1,NDS
					L=NFR*K
					GQ(L)=GQ(L)+RLD(LE,K)
				END DO
				CALL SUMLD(LE,GQ,NN,RD)
				IF(IM==0) CALL SUMRR(LE,STF,NV,SK)
				IF(IM==1) CALL SUMSS(LE,STF,NV,SK)
				IF(IM==2) CALL SUMKK(LE,STF,NV,SK,NEX)
			END DO
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				DO J=1,NDM
					K=JRR(IX,J)
					IF(K>0)THEN
						RD(K)=RD(K)+NOFCE(IX,J)
					END IF
				END DO
			END DO
			IF(IM==0)THEN
				CALL DECOMP(NN,NV,SK)
				CALL FORBACK(NV,SK,NN,RD)
			ELSE IF(IM==1)THEN
				CALL SOLVE(NV,SK,NN,RD)
			ELSE IF(IM==2)THEN
			    CALL SOL(SK,NV,NN,RD,NEX)
			END IF
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				CALL DDISCALC(IX,NN,RD)
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				CALL DSTNCALC(LE)
				CALL DSTSCALC(LE)
			END DO
			IF(IL<NSK.AND.ICON(1)/=0)THEN
				DO IE=1,ME
					LE=LUS(IE)
					IF(IST(LE)==-1) CYCLE
					CALL LDING(LE)
				END DO
			END IF
		END DO
		CALL ADRSLT(IP)
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				CALL YLDMCALC(LE)
			END DO
		END IF
	END DO
	DEALLOCATE(SK,RD,NEX)
	RETURN
END SUBROUTINE CALCEX


SUBROUTINE CALCFL(IP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::DT,RT,zd,hm
	INTEGER::IP,IE,IX,J,K,L,ME,LE,NS,ND,NJ,NK,IM,IL,NN,NV,NSK,LC,kk,ll
	REAL(KIND=8),DIMENSION(NDS)::GV,QV
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS)::KC
	REAL(KIND=8),DIMENSION(NDS*NDM)::GL,sf
	REAL(KIND=8),ALLOCATABLE::SK(:),RD(:)
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	INTEGER,ALLOCATABLE::NEX(:)
	ME=LQ(IP,2)
	ND=LQ(IP,3)
	NJ=LQ(IP,4)
	NK=LQ(IP,5)
	NS=LQ(IP,9)
	IM=ICON(10)
	IF(NS<1) NS=1
	DT=REAL(ND)/REAL(NS)
	RT=1.0/REAL(NS)
	NSK=1
	IF(ICON(1)/=0) NSK=2
	DO IX=1,NX
		DO J=1,NDM
			NOFCE(IX,J)=0.0
		END DO
		DO K=1,NFR
			NODIS(IX,K)=0.0
		END DO
	END DO
	IF(NJ/=0) CALL NOFCEIN(NJ)
	IF(NK/=0) CALL NODISIN(NK)
	DO IX=1,NX
 		DO J=1,NDM
			NOFCE(IX,J)=NOFCE(IX,J)*RT
		END DO
		DO K=1,NFR
			NODIS(IX,K)=NODIS(IX,K)*RT
		END DO
	END DO
	DO IE=1,ME
	    le=lus(ie)
		IF(IST(le).ne.2) CYCLE
		FORALL(L=1:NDS) GV(L)=0.0
		CALL GRVT(LE,GV)
		DO J=1,NDS
			K=NOD(LE,J)
			NOFCE(K,NDM)=NOFCE(K,NDM)+GV(J)*RT
		END DO
	END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   施加水压力
	DO J=1,NF
		IE=KSF(J,1)
		IF(IST(IE)==0.or.kgss(mat(ie),1)<0) CYCLE
		FORALL(L=1:NDS*NDM) SF(L)=0.0
		L=KSF(J,2)
		ZD=COR(L,NDM)
		DO K=3,5
			L=KSF(J,K)
			IF(COR(L,NDM)<ZD) ZD=COR(L,NDM)
		END DO
		IF(KSF(J,7)<0)THEN
			HM=MAX(HWU(IP,1),HWU(IP,2))
		ELSE IF(KSF(J,7)>0)THEN
			HM=MAX(HWD(IP,1),HWD(IP,2))
		END IF
		IF(ZD>=HM) CYCLE
		CALL SRFP(IP,J,SF)
		DO K=1,NDS
			KK=NOD(IE,K)
			LL=(K-1)*NDM
			DO L=1,NDM
				NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)*RT
			END DO
		END DO
	END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	ALLOCATE(NEX(NTEM))
	IF(IM==0)THEN
		CALL MRR(IP)
	ELSE IF(IM==1)THEN
		CALL MSS(IP)
	ELSE IF(IM==2)THEN
	    CALL MMM(IP)
	END IF
	NN=NDOF(IP)
	NV=NVOL(IP)
	ALLOCATE(SK(NV),RD(NN),NEX(NV))
	CALL MMKK(IP,NN,NV,NEX)
	DO LC=1,NS
		PRINT('(A10,2X,A8,2X,I2)'),'*****','SUBSTEP:',LC
		DO IE=1,NE
			DO J=1,NDS
				RLD(IE,J)=0.0
			END DO
		END DO
		IF(ICON(2)/=0)THEN
			DO J=1,NQ
				IE=KSQ(J,1)
				IF(IST(IE)==0) CYCLE
				FORALL(K=1:NDS) QV(K)=0.0
				CALL SRFQ(IP,J,QV)
				DO K=1,NDS
					RLD(IE,K)=RLD(IE,K)+QV(K)*DT
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				FORALL(L=1:NDS) QV(L)=0.0
				CALL VOLQ(LE,QV)
				DO L=1,NDS
					RLD(LE,L)=RLD(LE,L)+QV(L)*DT
				END DO
			END DO
		END IF
		DO IL=1,NSK
			DO J=1,NN
				RD(J)=0.0
			END DO
			DO K=1,NV
				SK(K)=0.0
			END DO
			DO IX=1,NX
				DO J=1,NDM
					DDIS(IX,J)=0.0
				END DO
				DTWH(IX)=0.0
			END DO
			DO IE=1,NE
				DO J=1,6
					DSTN(IE,J)=0.0
					DSTS(IE,J)=0.0
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				CALL EVSNCALC(LE)
				CALL STFNS(LE,STF,DT)
				FORALL(K=1:NDS*NFR) GQ(K)=0.0
				CALL BNDAPP(LE,STF,GQ)
				DO K=1,NDS*NFR
					GQ(K)=-GQ(K)
				END DO
				DO K=1,NDS
					L=NFR*K
					GQ(L)=GQ(L)+RLD(LE,K)
				END DO
				CALL SUMLD(LE,GQ,NN,RD)
				IF(IM==0) CALL SUMRR(LE,STF,NV,SK)
				IF(IM==1) CALL SUMSS(LE,STF,NV,SK)
				IF(IM==2) CALL SUMKK(LE,STF,NV,SK,NEX)
			END DO
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				DO J=1,NDM
					K=JRR(IX,J)
					IF(K>0)THEN
						RD(K)=RD(K)+NOFCE(IX,J)
					END IF
				END DO
			END DO
			do ix=1,nn
             rd(ix)=rd(ix)+1.0e-20
			end do
			print*,'ok'
			IF(IM==0)THEN
				CALL DECOMP(NN,NV,SK)
				CALL FORBACK(NV,SK,NN,RD)
			ELSE IF(IM==1)THEN
				CALL SOLVE(NV,SK,NN,RD)
			ELSE IF(IM==2)THEN
			    CALL SOL(SK,NV,NN,RD,NEX)
			END IF
			print*,'also ok!'
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				CALL DDISCALC(IX,NN,RD)
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				CALL DSTNCALC(LE)
				CALL DSTSCALC(LE)
			END DO
			IF(IL<NSK.AND.ICON(1)/=0)THEN
				DO IE=1,ME
					LE=LUS(IE)
					CALL LDING(LE)
				END DO
			END IF
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			IF(IST(LE)==2)THEN
				DO J=1,6
					DSTN(LE,J)=0.0
				END DO
				DO K=1,6
					DSTS(LE,K)=0.0
				END DO
			END IF
		END DO
		CALL ADRSLT(IP)
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				CALL YLDMCALC(LE)
			END DO
		END IF
	END DO
	DEALLOCATE(SK,RD,NEX)
	RETURN
END SUBROUTINE CALCFL

SUBROUTINE CALCWF(IP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::ZD,HM,HX,HN,VR,DT,RT
	INTEGER::IP,IE,IX,J,K,L,KK,LL,ME,LE,ND,NS,NJ,NK,NW,MW,	&
					IM,IL,NN,NV,NSK,LM,LW,LC
	REAL(KIND=8),DIMENSION(NDS)::GV,QV
	REAL(KIND=8),DIMENSION(NDS*NDM)::SF
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	REAL(KIND=8),ALLOCATABLE::SK(:),RD(:)
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	INTEGER,ALLOCATABLE::NEX(:)
	ME=LQ(IP,2)
	ND=LQ(IP,3)
	NJ=LQ(IP,4)
	NK=LQ(IP,5)
	NW=LQ(IP,6)  !/=1不考虑湿化；=1考虑湿化
	MW=LQ(IP,7)  !/=1不考虑浮力；=1考虑浮力
	NS=LQ(IP,9)
	IM=ICON(10)
	IF(NS<1) NS=1
	DT=REAL(ND)/REAL(NS)
	RT=1.0/REAL(NS)
	NSK=1
	IF(ICON(1)/=0) NSK=2
	DO IX=1,NX
		DO J=1,NDM
			NOFCE(IX,J)=0.0
		END DO
		DO K=1,NFR
			NODIS(IX,K)=0.0
		END DO
	END DO
	IF(NJ/=0) CALL NOFCEIN(NJ)
	IF(NK/=0) CALL NODISIN(NK)
	DO IX=1,NX
		DO J=1,NDM
			NOFCE(IX,J)=NOFCE(IX,J)*RT
		END DO
		DO K=1,NFR
			NODIS(IX,K)=NODIS(IX,K)*RT
		END DO
	END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!增加节点水头的计算程序
	IF(ICON(2)/=0)THEN
		DO J=1,NP1
			IE=PSJ(J,1)
		!	WRITE(40,*)IE
			IF(JST(IE)==0) CYCLE
			IF(PSJ(J,2)<0)THEN
				HM=SIGN(1.0,(HWU(IP,2)-HWU(IP,1)))
				HX=MAX(HWU(IP,2),HWU(IP,1))
				HN=MIN(HWU(IP,2),HWU(IP,1))
			ELSE IF(PSJ(J,2)>0)THEN
				HM=SIGN(1.0,(HWD(IP,2)-HWD(IP,1)))
				HX=MAX(HWD(IP,2),HWD(IP,1))
				HN=MIN(HWD(IP,2),HWD(IP,1))
			END IF
	!		DO K=2,5
	!			L=KSF(J,K)
				ZD=COR(IE,NDM)
				IF(ZD<HN)THEN
				    NODIS(IE,NFR)=(HX-HN)*HM*RT
				ELSE IF(ZD<HX)THEN
					NODIS(IE,NFR)=(HX-ZD)*HM*RT
				ELSE
					NODIS(IE,NFR)=0.0
				END IF
    !		END DO
		END DO
	END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   施加水压力
	DO J=1,NF
		IE=KSF(J,1)
		IF(IST(IE)==0.or.kgss(mat(ie),1)<0) CYCLE
		FORALL(L=1:NDS*NDM) SF(L)=0.0
		L=KSF(J,2)
		ZD=COR(L,NDM)
		DO K=3,5
			L=KSF(J,K)
			IF(COR(L,NDM)<ZD) ZD=COR(L,NDM)
		END DO
		IF(KSF(J,7)<0)THEN
			HM=MAX(HWU(IP,1),HWU(IP,2))
		ELSE IF(KSF(J,7)>0)THEN
			HM=MAX(HWD(IP,1),HWD(IP,2))
		END IF
		IF(ZD>=HM) CYCLE
		CALL SRFP(IP,J,SF)
!write(40,'(32f18.9)')(sf(l),l=1,32)
		DO K=1,NDS
			KK=NOD(IE,K)
			LL=(K-1)*NDM
			DO L=1,NDM
			
			if(l.eq.2)then
				NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)*RT
             else
                NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)*RT
			end if

		    END DO
		END DO
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!考虑湿化
	IF(NW/=0)THEN
		DO IE=1,ME
			LE=LUS(IE)
			LM=MAT(LE)
			LW=KGSS(LM,3)
			IF(LW==0) CYCLE
			ZD=0.0
			DO J=1,NDS
				K=NOD(LE,J)
				ZD=ZD+COR(K,NDM)
			END DO
			ZD=ZD/REAL(NDS)
			IF(KGSS(LM,5)==-1) HM=(ZD-HWU(IP,1))*(ZD-HWU(IP,2))
			IF(KGSS(LM,5)==+1) HM=(ZD-HWD(IP,1))*(ZD-HWD(IP,2))
			IF(HM>0.0) CYCLE
			FORALL(K=1:NDS*NDM) SF(K)=0.0
			CALL WETDFS(LE,LM,SF)                !调用湿化变形子程序
			DO K=1,NDS
				KK=NOD(LE,K)
				LL=NDM*(K-1)
				DO L=1,NDM
					NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)*RT
				END DO
			END DO
		END DO
	END IF
	IF(MW/=0)THEN                           !!!!!!!!!!!!!!!!!!!!考虑浮力
		DO IE=1,ME
			LE=LUS(IE)
			LM=MAT(LE)
			LW=KGSS(LM,4)
			!VR=CW(LM,15)
			vr=0.3
			VR=1.0-VR
			IF(LW==0) CYCLE
			!pause
			ZD=0.0
			DO J=1,NDS
				K=NOD(LE,J)
				ZD=ZD+COR(K,NDM)
			END DO
			ZD=ZD/REAL(NDS)
			IF(KGSS(LM,5)==-1) HM=(ZD-HWU(IP,1))*(ZD-HWU(IP,2))
			IF(KGSS(LM,5)==+1) HM=(ZD-HWD(IP,1))*(ZD-HWD(IP,2))
			IF(HM>0.0) CYCLE
			FORALL(K=1:NDS) GV(K)=0.0
			CALL GRVT(LE,GV)
			HM=GMAW/CS(LM,15)
			FORALL(L=1:NDS) GV(L)=-GV(L)*HM
			IF(KGSS(LM,5)==-1)THEN
				HM=SIGN(1.0,(HWU(IP,2)-HWU(IP,1)))
			ELSE IF(KGSS(LM,5)==+1)THEN
				HM=SIGN(1.0,(HWD(IP,2)-HWD(IP,1)))
			END IF
			DO J=1,NDS
				L=NOD(LE,J)
			!	NOFCE(L,NDM)=NOFCE(L,NDM)+HM*GV(J)*VR*RT
			END DO
		END DO
	END IF
!	ALLOCATE(NEX(NTEM))
	IF(IM==0)THEN
		CALL MRR(IP)
	ELSE IF(IM==1)THEN
		CALL MSS(IP)
	ELSE IF(IM==2)THEN
	    CALL MMM(IP)
	END IF
	NN=NDOF(IP)
	NV=NVOL(IP)
	ALLOCATE(SK(NV),RD(NN),NEX(NV))
	CALL MMKK(IP,NN,NV,NEX)
	DO LC=1,NS
		PRINT('(A10,2X,A8,2X,I2)'),'*****','SUBSTEP:',LC
		DO IE=1,NE
			DO J=1,NDS
				RLD(IE,J)=0.0
			END DO
		END DO
		IF(ICON(2)/=0)THEN
			DO J=1,NQ
				IE=KSQ(J,1)
				IF(IST(IE)==0) CYCLE
				FORALL(K=1:NDS) QV(K)=0.0
				CALL SRFQ(IP,J,QV)
				DO K=1,NDS
					RLD(IE,K)=RLD(IE,K)+QV(K)*DT
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				FORALL(L=1:NDS) QV(L)=0.0
				CALL VOLQ(LE,QV)
				DO L=1,NDS
					RLD(LE,L)=RLD(LE,L)+QV(L)*DT
				END DO
			END DO
		END IF
		DO IL=1,NSK
			DO J=1,NN
				RD(J)=0.0
			END DO
			DO K=1,NV
				SK(K)=0.0
			END DO
			DO IX=1,NX
				DO J=1,NDM
					DDIS(IX,J)=0.0
				END DO
				DTWH(IX)=0.0
			END DO
			DO IE=1,NE
				DO J=1,6
					DSTN(IE,J)=0.0
					DSTS(IE,J)=0.0
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				CALL EVSNCALC(LE)
				CALL STFNS(LE,STF,DT)
				FORALL(K=1:NDS*NFR) GQ(K)=0.0
				CALL BNDAPP(LE,STF,GQ)
				DO K=1,NDS*NFR
					GQ(K)=-GQ(K)
				END DO
				DO K=1,NDS
					L=NFR*K
					GQ(L)=GQ(L)+RLD(LE,K)
				END DO
				CALL SUMLD(LE,GQ,NN,RD)
				IF(IM==0) CALL SUMRR(LE,STF,NV,SK)
				IF(IM==1) CALL SUMSS(LE,STF,NV,SK)
				IF(IM==2) CALL SUMKK(LE,STF,NV,SK,NEX)
			END DO
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				DO J=1,NDM
					K=JRR(IX,J)
					IF(K>0)THEN
						RD(K)=RD(K)+NOFCE(IX,J)
					END IF
				END DO
			END DO
			do ix=1,nn
             rd(ix)=rd(ix)+1.0e-20
			end do
			print*,'ok'
			IF(IM==0)THEN
				CALL DECOMP(NN,NV,SK)
				CALL FORBACK(NV,SK,NN,RD)
			ELSE IF(IM==1)THEN
				CALL SOLVE(NV,SK,NN,RD)
			ELSE IF(IM==2)THEN
			    CALL SOL(SK,NV,NN,RD,NEX)
			END IF
			print*,'also ok!'
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				CALL DDISCALC(IX,NN,RD)
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				CALL DSTNCALC(LE)
				CALL DSTSCALC(LE)
			END DO
			IF(IL<NSK.AND.ICON(1)/=0)THEN
				DO IE=1,ME
					LE=LUS(IE)
					CALL LDING(LE)
				END DO
			END IF
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			IF(IST(LE)==2)THEN
				DO J=1,6
					DSTN(LE,J)=0.0
				END DO
				DO K=1,6
					DSTS(LE,K)=0.0
				END DO
			END IF
		END DO
		CALL ADRSLT(IP)
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				CALL YLDMCALC(LE)
			END DO
		END IF
	END DO
	DEALLOCATE(SK,RD,NEX)
	RETURN
END SUBROUTINE CALCWF


SUBROUTINE CALCQK(IP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::DT,RT,ACX,ACY,ACZ,AM
	INTEGER::IP,IE,IX,J,K,L,ME,LE,NS,ND,NJ,NK,IM,IL,NN,NV,NSK,LC
	REAL(KIND=8),DIMENSION(NDS)::GV,QV
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	REAL(KIND=8),ALLOCATABLE::SK(:),RD(:)
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	INTEGER,ALLOCATABLE::NEX(:)
	ME=LQ(IP,2)
	ND=LQ(IP,3)
	NJ=LQ(IP,4)
	NK=LQ(IP,5)
	NS=LQ(IP,9)
	IM=ICON(10)
	IF(NS<1) NS=1
	DT=REAL(ND)/REAL(NS)
	RT=1.0/REAL(NS)
	NSK=1
	IF(ICON(1)/=0) NSK=2
	DO IX=1,NX
		DO J=1,NDM
			NOFCE(IX,J)=0.0
		END DO
		DO K=1,NFR
			NODIS(IX,K)=0.0
		END DO
	END DO
	IF(NJ/=0) CALL NOFCEIN(NJ)
	IF(NK/=0) CALL NODISIN(NK)
	DO IX=1,NX
 		DO J=1,NDM
			NOFCE(IX,J)=NOFCE(IX,J)*RT
		END DO
		DO K=1,NFR
			NODIS(IX,K)=NODIS(IX,K)*RT
		END DO
	END DO
	READ(20,*) ACX,ACY,ACZ,AM
	CALL ACCDIS(IP,ACX,ACY,ACZ,AM)
	DO IE=1,ME
		LE=LUS(IE)
		FORALL(L=1:NDS) GV(L)=0.0
		CALL GRVT(LE,GV)
		FORALL(L=1:NDS) GV(L)=-GV(L)
		DO J=1,NDS
			K=NOD(LE,J)
			DO L=1,NDM
				NOFCE(K,L)=NOFCE(K,L)-GV(J)*RT*EQK(LE,L)
			END DO
		END DO
	END DO
!	ALLOCATE(NEX(NTEM))
	IF(IM==0)THEN
		CALL MRR(IP)
	ELSE IF(IM==1)THEN
		CALL MSS(IP)
	ELSE IF(IM==2)THEN
	    CALL MMM(IP)
	END IF
	NN=NDOF(IP)
	NV=NVOL(IP)
	ALLOCATE(SK(NV),RD(NN),NEX(NV))
	CALL MMKK(IP,NN,NV,NEX)
	DO LC=1,NS
		PRINT('(A10,2X,A8,2X,I2)'),'*****','SUBSTEP:',LC
		DO IE=1,NE
			DO J=1,NDS
				RLD(IE,J)=0.0
			END DO
		END DO
		IF(ICON(2)/=0)THEN
			DO J=1,NQ
				IE=KSQ(J,1)
				IF(IST(IE)==0) CYCLE
				FORALL(K=1:NDS) QV(K)=0.0
				CALL SRFQ(IP,J,QV)
				DO K=1,NDS
					RLD(IE,K)=RLD(IE,K)+QV(K)*DT
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				FORALL(L=1:NDS) QV(L)=0.0
				CALL VOLQ(LE,QV)
				DO L=1,NDS
					RLD(LE,L)=RLD(LE,L)+QV(L)*DT
				END DO
			END DO
		END IF
		DO IL=1,NSK
			DO J=1,NN
				RD(J)=0.0
			END DO
			DO K=1,NV
				SK(K)=0.0
			END DO
			DO IX=1,NX
				DO J=1,NDM
					DDIS(IX,J)=0.0
				END DO
				DTWH(IX)=0.0
			END DO
			DO IE=1,NE
				DO J=1,6
					DSTN(IE,J)=0.0
					DSTS(IE,J)=0.0
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				CALL EVSNCALC(LE)
				CALL STFNS(LE,STF,DT)
				FORALL(K=1:NDS*NFR) GQ(K)=0.0
				CALL BNDAPP(LE,STF,GQ)
				DO K=1,NDS*NFR
					GQ(K)=-GQ(K)
				END DO
				DO K=1,NDS
					L=NFR*K
					GQ(L)=GQ(L)+RLD(LE,K)
				END DO
				CALL SUMLD(LE,GQ,NN,RD)
				IF(IM==0) CALL SUMRR(LE,STF,NV,SK)
				IF(IM==1) CALL SUMSS(LE,STF,NV,SK)
				IF(IM==2) CALL SUMKK(LE,STF,NV,SK,NEX)
			END DO
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				DO J=1,NDM
					K=JRR(IX,J)
					IF(K>0)THEN
						RD(K)=RD(K)+NOFCE(IX,J)
					END IF
				END DO
			END DO
			IF(IM==0)THEN
				CALL DECOMP(NN,NV,SK)
				CALL FORBACK(NV,SK,NN,RD)
			ELSE IF(IM==1)THEN
				CALL SOLVE(NV,SK,NN,RD)
		    ELSE IF(IM==2)THEN
			    CALL SOL(SK,NV,NN,RD,NEX)
			END IF
			DO IX=1,NX
				IF(JR(IX)==0) CYCLE
				CALL DDISCALC(IX,NN,RD)
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				CALL DSTNCALC(LE)
				CALL DSTSCALC(LE)
			END DO
			IF(IL<NSK.AND.ICON(1)/=0)THEN
				DO IE=1,ME
					LE=LUS(IE)
					CALL LDING(LE)
				END DO
			END IF
		END DO
		CALL ADRSLT(IP)
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				CALL YLDMCALC(LE)
			END DO
		END IF
	END DO
	DEALLOCATE(SK,RD,NEX)
	RETURN
END SUBROUTINE CALCQK


SUBROUTINE CALCRH(IP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::DT,twhmax,tford,tback,flag,hx,hn,zd,vr,hm
	INTEGER::IP,IE,ME,LE,ND,NS,IM,IX,J,K,L,KK,LL,NN,NV,LM,LR,LC,i,nw,lw,mw,il
	REAL(KIND=8),DIMENSION(NDS)::QV,gv
	REAL(KIND=8),DIMENSION(NDS*NDM)::SF
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	REAL(KIND=8),ALLOCATABLE::SK(:),RD(:)
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	INTEGER,ALLOCATABLE::NEX(:)
	ME=LQ(IP,2)
	ND=LQ(IP,3)
	NS=LQ(IP,9)
	nw=lq(ip,6)
	mw=lq(ip,7)
	DT=REAL(ND)/REAL(NS)
	IM=ICON(10)
!	ALLOCATE(NEX(NTEM))
	IF(IM==0)THEN
		CALL MRR(IP)
	ELSE IF(IM==1)THEN
		CALL MSS(IP)
	ELSE IF(IM==2)THEN
	    CALL MMM(IP)
	END IF
	NN=NDOF(IP)
	NV=NVOL(IP)
	ALLOCATE(SK(NV),RD(NN),NEX(NV))
	CALL MMKK(IP,NN,NV,NEX)

	twhmax=234     !竣工时最高水位
    hwu(ip,2)=234  !竣工时初始水位
	DO LC=1,NS
		PRINT('(A10,2X,A8,2X,I2)'),'*****','SUBSTEP:',LC
		DO IX=1,NX
			DO J=1,NDM
				NOFCE(IX,J)=0.0
			END DO
			DO K=1,NFR
				NODIS(IX,K)=0.0
			END DO
		END DO
		DO IE=1,NE
			DO J=1,NDS
				RLD(IE,J)=0.0
			END DO
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			IF(IST(LE)==-1) CYCLE
			LM=MAT(LE)
			LR=KGSS(LM,6)
			IF(LR==0) CYCLE
			FORALL(J=1:NDS*NDM) SF(J)=0.0
			CALL RLGDF(LE,LM,DT,SF)         !调用流变变形子程序
			DO K=1,NDS
				KK=NOD(LE,K)
				LL=NDM*(K-1)
				DO L=1,NDM
					NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)
				END DO
			END DO
		END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!在此考虑运行过程中的水位变化情况
! ntlb, shuiwei(i,1),时间； shuiwei(i,2),库水位
totaltime=totaltime+dt      !当前时间
do i=1,ntlb
tford=shuiwei(i,1)
tback=shuiwei(i+1,1)
flag=(totaltime-tford)*(totaltime-tback)
if(flag.lt.0)then
hwu(ip,1)=hwu(ip,2)
hwu(ip,2)=(totaltime-tford)/(tback-tford)*(shuiwei(i+1,2)-shuiwei(i,2))+shuiwei(i,2) !根据插值计算本级水位
goto 111
end if
end do
111 continue
print*,hwu(ip,1),hwu(ip,2),totaltime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!增加节点水头的计算程序
	IF(ICON(2)/=0)THEN
		DO J=1,NP1
			IE=PSJ(J,1)
			IF(JST(IE)==0) CYCLE
			IF(PSJ(J,2)<0)THEN
				HM=SIGN(1.0,(HWU(IP,2)-HWU(IP,1)))
				HX=MAX(HWU(IP,2),HWU(IP,1))
				HN=MIN(HWU(IP,2),HWU(IP,1))
			ELSE IF(PSJ(J,2)>0)THEN
				HM=SIGN(1.0,(HWD(IP,2)-HWD(IP,1)))
				HX=MAX(HWD(IP,2),HWD(IP,1))
				HN=MIN(HWD(IP,2),HWD(IP,1))
			END IF
				ZD=COR(IE,NDM)
				IF(ZD<HN)THEN
				 NODIS(IE,NFR)=(HX-HN)*HM
				elseIF(ZD<HX)THEN
                 NODIS(IE,NFR)=(HX-ZD)*HM
				ELSE
				 NODIS(IE,NFR)=0.0
				END IF
		END DO
	END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   施加水压力
	DO J=1,NF
		IE=KSF(J,1)
		IF(IST(IE)==0.or.kgss(mat(ie),1)<0) CYCLE
        if(ksf(j,7).eq.-2)cycle                     !!!!!!!水位超过围堰高程后该水压力自动消失
		FORALL(L=1:NDS*NDM) SF(L)=0.0
		L=KSF(J,2)
		ZD=COR(L,NDM)
		DO K=3,5
			L=KSF(J,K)
			IF(COR(L,NDM)<ZD) ZD=COR(L,NDM)
		END DO
		IF(KSF(J,7)<0)THEN
			HM=MAX(HWU(IP,1),HWU(IP,2))
		ELSE IF(KSF(J,7)>0)THEN
			HM=MAX(HWD(IP,1),HWD(IP,2))
		END IF
		IF(ZD>=HM) CYCLE
		CALL SRFP(IP,J,SF)
		DO K=1,NDS
			KK=NOD(IE,K)
			LL=(K-1)*NDM
			DO L=1,NDM
				NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)
			END DO
		END DO
	END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     考虑湿化
	IF(NW/=0)THEN
		DO IE=1,ME
			LE=LUS(IE)
			LM=MAT(LE)
			LW=KGSS(LM,3)
			IF(LW==0) CYCLE
			ZD=0.0
			DO J=1,NDS
				K=NOD(LE,J)
				ZD=ZD+COR(K,NDM)
			END DO
			ZD=ZD/REAL(NDS)
			IF(KGSS(LM,5)==-1) HM=(ZD-HWU(IP,1))*(ZD-HWU(IP,2))
			IF(KGSS(LM,5)==+1) HM=(ZD-HWD(IP,1))*(ZD-HWD(IP,2))
			IF(HM>0.0.or.zd.le.twhmax) CYCLE     !不考虑反复蓄水湿化变形
			FORALL(K=1:NDS*NDM) SF(K)=0.0
			CALL WETDFS(LE,LM,SF)                !调用湿化变形子程序
			DO K=1,NDS
				KK=NOD(LE,K)
				LL=NDM*(K-1)
				DO L=1,NDM
					NOFCE(KK,L)=NOFCE(KK,L)+SF(LL+L)
				END DO
			END DO
		END DO
	END IF
	IF(MW/=0)THEN                  !!!!!!!!!!!!!!!!!!!!考虑浮力
		DO IE=1,ME
			LE=LUS(IE)
			LM=MAT(LE)
			LW=KGSS(LM,4)
			VR=CW(LM,15)
			VR=1.0-VR
			IF(LW==0) CYCLE
			ZD=0.0
			DO J=1,NDS
				K=NOD(LE,J)
				ZD=ZD+COR(K,NDM)
			END DO
			ZD=ZD/REAL(NDS)
			IF(KGSS(LM,5)==-1) HM=(ZD-HWU(IP,1))*(ZD-HWU(IP,2))
			IF(KGSS(LM,5)==+1) HM=(ZD-HWD(IP,1))*(ZD-HWD(IP,2))
			IF(HM>0.0) CYCLE
			FORALL(K=1:NDS) GV(K)=0.0
			CALL GRVT(LE,GV)
			HM=GMAW/CS(LM,15)
			FORALL(L=1:NDS) GV(L)=-GV(L)*HM
			IF(KGSS(LM,5)==-1)THEN
				HM=SIGN(1.0,(HWU(IP,2)-HWU(IP,1)))
			ELSE IF(KGSS(LM,5)==+1)THEN
				HM=SIGN(1.0,(HWD(IP,2)-HWD(IP,1)))
			END IF
			DO J=1,NDS
				L=NOD(LE,J)
				NOFCE(L,NDM)=NOFCE(L,NDM)+HM*GV(J)*VR
			END DO
		END DO
	END IF
!!!!!!!!计算历史最高蓄水位
   if(hwu(ip,2).ge.twhmax)twhmax=hwu(ip,2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF(ICON(2)/=0)THEN
			DO J=1,NQ
				IE=KSQ(J,1)
				IF(IST(IE)==0) CYCLE
				FORALL(K=1:NDS) QV(K)=0.0
				CALL SRFQ(IP,J,QV)
				DO K=1,NDS
					RLD(IE,K)=RLD(IE,K)+QV(K)*DT
				END DO
			END DO
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				FORALL(L=1:NDS) QV(L)=0.0
				CALL VOLQ(LE,QV)
				DO L=1,NDS
					RLD(LE,L)=RLD(LE,L)+QV(L)*DT
				END DO
			END DO
		END IF
	do il=1,2
		DO J=1,NN
			RD(J)=0.0
		END DO
		DO K=1,NV
			SK(K)=0.0
		END DO
		DO IX=1,NX
			DO J=1,NDM
				DDIS(IX,J)=0.0
			END DO
			DTWH(IX)=0.0
		END DO
		DO IE=1,NE
			DO J=1,6
				DSTN(IE,J)=0.0
				DSTS(IE,J)=0.0
			END DO
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			IF(IST(LE)==-1) CYCLE
			CALL EVSNCALC(LE)
			CALL STFNS(LE,STF,DT)
			FORALL(K=1:NDS*NFR) GQ(K)=0.0
			CALL BNDAPP(LE,STF,GQ)
			DO K=1,NDS*NFR
				GQ(K)=-GQ(K)
			END DO
			DO K=1,NDS
				L=NFR*K
				GQ(L)=GQ(L)+RLD(LE,K)
			END DO
			CALL SUMLD(LE,GQ,NN,RD)
			IF(IM==0) CALL SUMRR(LE,STF,NV,SK)
			IF(IM==1) CALL SUMSS(LE,STF,NV,SK)
			IF(IM==2) CALL SUMKK(LE,STF,NV,SK,NEX)
		END DO
		DO IX=1,NX
			IF(JR(IX)==0) CYCLE
			DO J=1,NDM
				K=JRR(IX,J)
				IF(K>0)THEN
					RD(K)=RD(K)+NOFCE(IX,J)
				END IF
			END DO
		END DO
		IF(IM==0)THEN
			CALL DECOMP(NN,NV,SK)
			CALL FORBACK(NV,SK,NN,RD)
		ELSE IF(IM==1)THEN
			CALL SOLVE(NV,SK,NN,RD)
		ELSE IF(IM==2)THEN
			    CALL SOL(SK,NV,NN,RD,NEX)
		END IF
		DO IX=1,NX
			IF(JR(IX)==0) CYCLE
			CALL DDISCALC(IX,NN,RD)
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			IF(IST(LE)==-1) CYCLE
			CALL DSTNCALC(LE)
			CALL DSTSCALC(LE)
		END DO
		IF(il.lt.2.and.ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				CALL LDING(LE)
			END DO
		END IF
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO IE=1,NE
			DO J=1,6
				RHSTN(IE,J)=0.0                  !流变产生的应变归零
			END DO
		END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CALL ADRSLT(IP)   !计算应力
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				IF(IST(LE)==-1) CYCLE
				CALL YLDMCALC(LE)
			END DO
		END IF
		if(lc/4*4.eq.lc) call svrst(ip)               !每级输出结果
	END DO
	DEALLOCATE(SK,RD,NEX)
	RETURN
END SUBROUTINE CALCRH


SUBROUTINE CALCSD(IP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::DT
	INTEGER::IP,IE,ME,LE,ND,NS,IM,IX,J,K,L,NN,NV,LC
	REAL(KIND=8),DIMENSION(NDS)::QV
	REAL(KIND=8),DIMENSION(NDS*NDM)::SF
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	REAL(KIND=8),ALLOCATABLE::SK(:),RD(:)
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	INTEGER,ALLOCATABLE::NEX(:)
	ME=LQ(IP,2)
	ND=LQ(IP,3)
	NS=LQ(IP,9)
	DT=REAL(ND)/REAL(NS)
	IM=ICON(10)
!	ALLOCATE(NEX(NTEM))
	IF(IM==0)THEN
		CALL MRR(IP)
	ELSE IF(IM==1)THEN
		CALL MSS(IP)
	ELSE IF(IM==2)THEN
	    CALL MMM(IP)
	END IF
	NN=NDOF(IP)
	NV=NVOL(IP)
	ALLOCATE(SK(NV),RD(NN),NEX(NV))
	CALL MMKK(IP,NN,NV,NEX)
	DO LC=1,NS
		PRINT('(A10,2X,A8,2X,I2)'),'*****','SUBSTEP:',LC
		DO IX=1,NX
			DO J=1,NDM
				NOFCE(IX,J)=0.0
			END DO
			DO K=1,NFR
				NODIS(IX,K)=0.0
			END DO
		END DO
		DO IE=1,NE
			DO J=1,NDS
				RLD(IE,J)=0.0
			END DO
		END DO
		DO J=1,NQ
			IE=KSQ(J,1)
			IF(IST(IE)==0) CYCLE
			FORALL(K=1:NDS) QV(K)=0.0
			CALL SRFQ(IP,J,QV)
			DO K=1,NDS
				RLD(IE,K)=RLD(IE,K)+QV(K)*DT
			END DO
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			FORALL(L=1:NDS) QV(L)=0.0
			CALL VOLQ(LE,QV)
			DO L=1,NDS
				RLD(LE,L)=RLD(LE,L)+QV(L)*DT
			END DO
		END DO
		DO J=1,NN
			RD(J)=0.0
		END DO
		DO K=1,NV
			SK(K)=0.0
		END DO
		DO IX=1,NX
			DO J=1,NDM
				DDIS(IX,J)=0.0
			END DO
			DTWH(IX)=0.0
		END DO
		DO IE=1,NE
			DO J=1,6
				DSTN(IE,J)=0.0
				DSTS(IE,J)=0.0
			END DO
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			CALL EVSNCALC(LE)
			CALL STFNS(LE,STF,DT)
			FORALL(K=1:NDS*NFR) GQ(K)=0.0
			CALL BNDAPP(LE,STF,GQ)
			DO K=1,NDS*NFR
				GQ(K)=-GQ(K)
			END DO
			DO K=1,NDS
				L=NFR*K
				GQ(L)=GQ(L)+RLD(LE,K)
			END DO
			CALL SUMLD(LE,GQ,NN,RD)
			IF(IM==0) CALL SUMRR(LE,STF,NV,SK)
			IF(IM==1) CALL SUMSS(LE,STF,NV,SK)
			IF(IM==2) CALL SUMKK(LE,STF,NV,SK,NEX)
		END DO
		DO IX=1,NX
			IF(JR(IX)==0) CYCLE
			DO J=1,NDM
				K=JRR(IX,J)
				IF(K>0)THEN
					RD(K)=RD(K)+NOFCE(IX,J)
				END IF
			END DO
		END DO
		IF(IM==0)THEN
			CALL DECOMP(NN,NV,SK)
			CALL FORBACK(NV,SK,NN,RD)
		ELSE IF(IM==1)THEN
			CALL SOLVE(NV,SK,NN,RD)
		ELSE IF(IM==2)THEN
			    CALL SOL(SK,NV,NN,RD,NEX)
		END IF
		DO IX=1,NX
			IF(JR(IX)==0) CYCLE
			CALL DDISCALC(IX,NN,RD)
		END DO
		DO IE=1,ME
			LE=LUS(IE)
			CALL DSTNCALC(LE)
			CALL DSTSCALC(LE)
		END DO
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				CALL LDING(LE)
			END DO
		END IF
		CALL ADRSLT(IP)
		IF(ICON(1)/=0)THEN
			DO IE=1,ME
				LE=LUS(IE)
				CALL YLDMCALC(LE)
			END DO
		END IF
		CALL SVRST(IP)
	END DO
	DEALLOCATE(SK,RD,NEX)
	RETURN
END SUBROUTINE CALCSD


SUBROUTINE PREPROC
	USE COMM
	IMPLICIT NONE
	INTEGER::I,J,NR,IR,KT,NN,NV,IM,KG,ie,mt
	REAL(KIND=8),DIMENSION(NDS)::GV
	REAL(KIND=8),DIMENSION(NGS)::HA,XB
	REAL(KIND=8),DIMENSION(NGF)::HC,XD
OPEN(233,FILE='jswg.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!	integer,dimension(ntem)::ndex
!	integer,dimension(ntem)::ndex
!!!!!!  adjust interface elememt
write(233,*)ne,nx
    do ie=1,ne 
     mt=mat(ie)
	 if(kgss(mt,1).lt.0)call ajtele(ie)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!
     write(233,'(12i8)')ie,(nod(ie,j),j=1,nds),mat(ie),jad(ie),jrd(ie)
	end do

	do i=1,nx
      write(233,'(i8,3f12.4)')i,(cor(i,j),j=1,ndm)
	end do

	IM=ICON(10)
	NR=0

	CALL BNDCHEK(NR)
	IF(NR/=0)THEN
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE PREPROC:'
		WRITE(40,*) 'THE BOUNDARY CONDITION CONFLICTS!'
		STOP
	END IF
	DO I=1,NM
		KG=KGSS(I,1)
		IF(KG==13.OR.KG==14)THEN
			IF(IM==0)THEN
				WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE PREPROC:'
				WRITE(40,*) 'THE STIFFNESS MATRIX IS NOT SYSMETRIC'
				WRITE(40,*) 'PLEASE RESET THE METHOD TO SOLVING LINEAR EQUATIONS!'
				STOP
			END IF
		END IF
	END DO
	CALL CONFIG(NGS,HA,XB)
	DO I=1,NGS
		WGH(I)=HA(I)
		WXYZ(I)=XB(I)
	END DO
	CALL CONFIG(NGF,HC,XD)
	DO I=1,NGF
		VGH(I)=HC(I)
		VXYZ(I)=XD(I)
	END DO
	NR=0
	DO I=1,NE
		KT=MAT(I)
		KT=KGSS(KT,1)
		IF(KT>=0) CALL CHECI(I,IR)
		IF(KT==-1) CALL CHECJ(I,IR)
		IF(KT==-2) CALL CHECJ(I,IR)
		IF(KT==-3) CALL CHECJ(I,IR)
		IF(KT==-4) CALL CHECK(I,IR)
		NR=NR+IR
	END DO

	IF(NR>0)THEN
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE PREPROC'
		WRITE(40,*) 'FOR SHAPE OF ELEMENT IS BAD OR SEQUENCE OF NODE IS WRONG'
		STOP
	END IF
	WRITE(30,10)
	DO I=1,NE
		FORALL(J=1:NDS) GV(J)=0.0
		CALL GRVT(I,GV)
		WRITE(30,20) I,(GV(J),J=1,NDS)
	END DO
	WRITE(30,30)
	DO I=1,NP
		CALL MQQ(I)
		IF(IM==0) CALL MRR(I)
		IF(IM==1) CALL MSS(I)
		if(im==2) call mmm(i)
		NN=NDOF(I)
		NV=NVOL(I)
		WRITE(30,40) I,NN,NV
	END DO
10	FORMAT(/,'GRAVITY OF ELEMENTS GV(NDS):')
20	FORMAT(I5,8F12.3)
30	FORMAT(/,'TOTAL EQUATIONS INFORMATION:')
40	FORMAT('STEP=',I4,5X,'NN=',I5,5X,'NV=',I12)
	RETURN
END SUBROUTINE PREPROC

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ajtele(ie)
use comm
implicit none
integer::ie,temp,l,kk,kl,k1,k2,k3,k4,k5,k6,k7,k8,j
real(kind=8)::a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
real(kind=8)::s1,s2,s3,mins
k1=nod(ie,1)
k2=nod(ie,2)
k3=nod(ie,3)
k4=nod(ie,4)
k5=nod(ie,5)
k6=nod(ie,6)
k7=nod(ie,7)
k8=nod(ie,8)
a1=(cor(k1,1)+cor(k4,1)+cor(k3,1)+cor(k2,1))/4-(cor(k5,1)+cor(k8,1)+cor(k7,1)+cor(k6,1))/4
a2=(cor(k1,2)+cor(k4,2)+cor(k3,2)+cor(k2,2))/4-(cor(k5,2)+cor(k8,2)+cor(k7,2)+cor(k6,2))/4
a3=(cor(k1,3)+cor(k4,3)+cor(k3,3)+cor(k2,3))/4-(cor(k5,3)+cor(k8,3)+cor(k7,3)+cor(k6,3))/4
s1=sqrt(a1**2+a2**2+a3**2)
b1=(cor(k1,1)+cor(k4,1)+cor(k8,1)+cor(k5,1))/4-(cor(k2,1)+cor(k3,1)+cor(k7,1)+cor(k6,1))/4
b2=(cor(k1,2)+cor(k4,2)+cor(k8,2)+cor(k5,2))/4-(cor(k2,2)+cor(k3,2)+cor(k7,2)+cor(k6,2))/4
b3=(cor(k1,3)+cor(k4,3)+cor(k8,3)+cor(k5,3))/4-(cor(k2,3)+cor(k3,3)+cor(k7,3)+cor(k6,3))/4
s2=sqrt(b1**2+b2**2+b3**2)
c1=(cor(k1,1)+cor(k2,1)+cor(k6,1)+cor(k5,1))/4-(cor(k3,1)+cor(k4,1)+cor(k8,1)+cor(k7,1))/4
c2=(cor(k1,2)+cor(k2,2)+cor(k6,2)+cor(k5,2))/4-(cor(k3,2)+cor(k4,2)+cor(k8,2)+cor(k7,2))/4
c3=(cor(k1,3)+cor(k2,3)+cor(k6,3)+cor(k5,3))/4-(cor(k3,3)+cor(k4,3)+cor(k8,3)+cor(k7,3))/4
s3=sqrt(c1**2+c2**2+c3**2)


mins=min(s1,s2)
mins=min(mins,s3)
if(mins.eq.s1)then
jplane(ie)=1
nod(ie,1)=nod2(ie,1)
nod(ie,2)=nod2(ie,2)
nod(ie,3)=nod2(ie,3)
nod(ie,4)=nod2(ie,4)

nod(ie,5)=nod2(ie,5)
nod(ie,6)=nod2(ie,6)
nod(ie,7)=nod2(ie,7)
nod(ie,8)=nod2(ie,8)

end if
if(mins.eq.s2)then
jplane(ie)=2

nod(ie,1)=nod2(ie,2)
nod(ie,2)=nod2(ie,6)
nod(ie,3)=nod2(ie,7)
nod(ie,4)=nod2(ie,3)

nod(ie,5)=nod2(ie,1)
nod(ie,6)=nod2(ie,5)
nod(ie,7)=nod2(ie,8)
nod(ie,8)=nod2(ie,4)


elseif(mins.eq.s3)then

jplane(ie)=3
nod(ie,1)=nod2(ie,6)
nod(ie,2)=nod2(ie,2)
nod(ie,3)=nod2(ie,1)
nod(ie,4)=nod2(ie,5)

nod(ie,5)=nod2(ie,7)
nod(ie,6)=nod2(ie,3)
nod(ie,7)=nod2(ie,4)
nod(ie,8)=nod2(ie,8)

end if

return
end subroutine ajtele

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SEARCHFILE(MAINFILE)
	IMPLICIT NONE
	CHARACTER(LEN=30)::MAINFILE
	LOGICAL::ALIVE
	DO
		PRINT*,'-----PLEASE INPUT THE CONTROL FILE NAME-----'
		READ*,MAINFILE
		INQUIRE(FILE=MAINFILE,EXIST=ALIVE)
		IF(ALIVE) EXIT
		PRINT*,'------CAN NOT FIND THE SPECIFIED FILE------'
	END DO
	RETURN
END SUBROUTINE SEARCHFILE


SUBROUTINE PRTIME
	IMPLICIT NONE
	CHARACTER(LEN=8)::CHAR_TIME
	CALL TIME(CHAR_TIME)
	PRINT*,CHAR_TIME
	RETURN
END SUBROUTINE PRTIME


SUBROUTINE FILEOPEN(MAINFILE)
	USE COMM
	IMPLICIT NONE
	CHARACTER(LEN=30)::MAINFILE
	OPEN(10,FILE=MAINFILE,FORM='FORMATTED',STATUS='OLD')
	OPEN(20,FILE='elenod.dat',FORM='FORMATTED',STATUS='OLD')
	OPEN(30,FILE='INFORM.CHK',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(40,FILE='WARNING.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(50,FILE='DASAC1.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(60,FILE='DASAC2.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(70,FILE='DASAC3.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(80,FILE='DASAC4.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(90,FILE='DASAC5.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(100,FILE='DASAC6.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN(217,FILE='COFFICESEEP.DAT',FORM='FORMATTED',STATUS='UNKNOWN')
!	OPEN(218,FILE='scqx.dat',FORM='FORMATTED',STATUS='UNKNOWN')
	RETURN
END SUBROUTINE FILEOPEN


SUBROUTINE FILECLOSE
	CLOSE(10)
	CLOSE(20)
	CLOSE(30)
	CLOSE(40)
	CLOSE(50)
	CLOSE(60)
	CLOSE(70)
	CLOSE(80)
	CLOSE(90)
	CLOSE(100)
	CLOSE(217)
!	CLOSE(218)
	RETURN
END SUBROUTINE FILECLOSE


SUBROUTINE DATAINP
	USE COMM
	IMPLICIT NONE
	INTEGER::I,J,K,L,IW,mt
	real(kind=8)::hei
!	OPEN(222,FILE='c2.DAT',FORM='FORMATTED',STATUS='OLD')
	READ(10,*) (ICON(I),I=1,10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111!!
	READ(10,*) NP,NM,NE,NX,NF,NQ,NF1,NP1 !,hsur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO I=1,NP
		READ(10,*) IW,(LQ(I,J),J=1,10),(HWU(I,K),K=1,2),(HWD(I,L),L=1,2)           !,conh(i)
	END DO
	DO I=1,NM
		READ(10,*) IW,(KGSS(I,J),J=1,7)
				print*,iw
	END DO	
	DO I=1,NM
		READ(10,*) IW,(CS(I,J),J=1,15)
		print*,iw,cs(i,1)
	END DO

	IF(ICON(2)/=0)THEN
		DO I=1,NM
		forall(j=1:15)ck(i,j)=0.0
			READ(10,*) IW,ck(i,1)    !,ck(i,15)          !(CK(I,J),J=1,15)
			ck(i,5)=ck(i,1)
			ck(i,9)=ck(i,1)
			if(kgss(i,1).lt.0)ck(i,9)=ck(i,1)/100
		END DO
	END IF
	IF(ICON(3)/=0)THEN   !湿化参数
		DO I=1,NM
			READ(10,*) IW,(CW(I,J),J=1,15)
		END DO
	END IF
	IF(ICON(4)/=0)THEN   !流变参数
		DO I=1,NM
			READ(10,*) IW,(CR(I,J),J=1,15)
		END DO
	END IF
		print*,'coefficient input finished!'

	read(20,*)ne
	k=0
	DO I=1,NE
		READ(20,*) IW,(NOD(I,J),J=1,NDS),MAT(I),JAD(I),jrd(i),SRD(I)
		forall(j=1:nds)nod2(i,j)=nod(i,j)
        jplane(i)=0
	END DO

	read(20,*)nx
	DO I=1,NX
		READ(20,*) IW,(COR(I,J),J=1,NDM),(ID(I,K),K=1,NFR)
		do k=1,nfr
        ig(i,k)=id(i,k)
		end do
	END DO

	read(20,*)nf
	DO I=1,NF
		READ(20,*) IW,(KSF(I,J),J=1,NFS+3)
	END DO
    !!!!!!!!!!!!输入排水节点与第二水压力的数据
    DO I=1,NP1
        READ(20,*) IW,(PSJ(I,J),J=1,2)
	END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO I=1,NQ
		READ(20,*) IW,(KSQ(I,J),J=1,NFS+2),DSQ(I)
	END DO

	DO I=1,NE
		mt=mat(i)
		mt=kgss(mt,1)
	    forall(j=1:6)str(i,j)=0
        if(mt.ge.0)then
         str(i,1)=100
		 str(i,2)=100
		 str(i,3)=200
		 else
		 str(i,3)=-200
		end if
	END DO

    DO I=1,NX
        gg(i,nfr)=cor(i,ndm)
    END DO

	RETURN
END SUBROUTINE DATAINP


SUBROUTINE NOFCEIN(MX)
	USE COMM
	IMPLICIT NONE
	INTEGER::MX,IX,JX,IW,LL
	DO IX=1,MX
		READ(20,*) IW,LL,((NOFCE(LL,JX)),JX=1,NDM)
	END DO
	RETURN
END SUBROUTINE NOFCEIN


SUBROUTINE NODISIN(MX)
	USE COMM
	IMPLICIT NONE
	INTEGER::MX,IX,JX,IW,LL
	DO IX=1,MX
		READ(20,*) IW,LL,(NODIS(LL,JX),JX=1,NFR)
	END DO
	RETURN
END SUBROUTINE NODISIN


SUBROUTINE BNDCHEK(NR)
	USE COMM
	IMPLICIT NONE
	CHARACTER(LEN=10)::CH
	INTEGER::NB,IE,J,K,L,IC,NR,jk,lk
	INTEGER,DIMENSION(NFS)::IA,IB
	INTEGER,ALLOCATABLE::KAB(:,:)
	NB=NF+NQ
	IF(NB==0) RETURN
	ALLOCATE(KAB(NB,6))
	L=0
	DO J=1,NF
		L=L+1

		DO K=1,6
			KAB(L,K)=KSF(J,K)
		END DO
	END DO
	DO J=1,NQ
		L=L+1
		DO K=1,6
			KAB(L,K)=KSQ(J,K)
		END DO
	END DO

	np1=0        !水头结点

	DO L=1,NB
		IE=KAB(L,1)
		if(kgss(mat(ie),1)<0)cycle
		IF(L<=NF)THEN
			CH='STRESS BND'
		ELSE
			CH='FLUX BND'
		END IF
		DO J=2,5
			K=J-1
			IA(K)=KAB(L,J)
		END DO
		IF(KAB(L,6)==-1)THEN
			IB(1)=NOD(IE,1)
			IB(2)=NOD(IE,2)
			IB(3)=NOD(IE,3)
			IB(4)=NOD(IE,4)
		ELSE IF(KAB(L,6)==1)THEN
			IB(1)=NOD(IE,5)
			IB(2)=NOD(IE,6)
			IB(3)=NOD(IE,7)
			IB(4)=NOD(IE,8)
		ELSE IF(KAB(L,6)==-2)THEN
			IB(1)=NOD(IE,2)
			IB(2)=NOD(IE,3)
			IB(3)=NOD(IE,7)
			IB(4)=NOD(IE,6)
		ELSE IF(KAB(L,6)==2)THEN
			IB(1)=NOD(IE,1)
			IB(2)=NOD(IE,4)
			IB(3)=NOD(IE,8)
			IB(4)=NOD(IE,5)
		ELSE IF(KAB(L,6)==-3)THEN
			IB(1)=NOD(IE,1)
			IB(2)=NOD(IE,2)
			IB(3)=NOD(IE,6)
			IB(4)=NOD(IE,5)
		ELSE IF(KAB(L,6)==3)THEN
			IB(1)=NOD(IE,4)
			IB(2)=NOD(IE,3)
			IB(3)=NOD(IE,7)
			IB(4)=NOD(IE,8)
		END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   修改水压力面上结点号
        forall(jk=1:4)ia(jk)=ib(jk)
        forall(jk=1:4)ksf(l,jk+1)=ib(jk)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CALL IAIB(NFS,IA,IB,IC)
		IF(IC==0)THEN
			NR=NR+1
			WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE BNDCHEK:',IE
			WRITE(40,*) 'FOR NODES ARE DIFFERENT.',CH
		END IF
	END DO
	DEALLOCATE(KAB)
	RETURN
END SUBROUTINE BNDCHEK


SUBROUTINE IAIB(N,IA,IB,IC)
	IMPLICIT NONE
	INTEGER::I,J,N,IC,II,JJ,NA,NB
	INTEGER,DIMENSION(N)::IA,IB
	IC=0
	DO I=1,N-1
		DO J=I+1,N
			IF(IA(I)==IA(J)) IA(I)=0
			IF(IB(I)==IB(J)) IB(I)=0
		END DO
	END DO
	NA=0
	NB=0
	DO I=1,N
		IF(IA(I)/=0)THEN
			NA=NA+1
			IA(NA)=IA(I)
		END IF
		IF(IB(I)/=0)THEN
			NB=NB+1
			IB(NB)=IB(I)
		END IF
	END DO
	IF(NA/=NB) RETURN
	DO I=1,NA
		II=IA(I)
		DO J=1,NB
			JJ=IB(J)
			IF(JJ/=II) CYCLE
			IC=IC+1
			EXIT
		END DO
	END DO
	IF(IC/=NA) IC=0
	RETURN
END SUBROUTINE IAIB


SUBROUTINE DATAOUT
	USE COMM
	IMPLICIT NONE
	INTEGER::I,J,K,L
	WRITE(30,5)
	WRITE(30,10) (ICON(I),I=1,10)
	WRITE(30,15)
	WRITE(30,20) NP,NM,NE,NX,NF,NQ
	WRITE(30,25)
	DO I=1,NP
		WRITE(30,30) I,(LQ(I,J),J=1,10),(HWU(I,K),K=1,2),(HWD(I,L),L=1,2)
	END DO
	WRITE(30,35)
	DO I=1,NM
		WRITE(30,40) I,(KGSS(I,J),J=1,7)
	END DO
	WRITE(30,45)
	DO I=1,NM
		WRITE(30,50) I,(CS(I,J),J=1,15)
	END DO
	IF(ICON(2)/=0)THEN
		WRITE(30,55)
		DO I=1,NM
			WRITE(30,60) I,(CK(I,J),J=1,15)
		END DO
	END IF
	IF(ICON(3)/=0)THEN
		WRITE(30,65)
		DO I=1,NM
			WRITE(30,70) I,(CW(I,J),J=1,15)
		END DO
	END IF
	IF(ICON(4)/=0)THEN
		WRITE(30,75)
		DO I=1,NM
			WRITE(30,80) I,(CR(I,J),J=1,15)
		END DO
	END IF
	WRITE(30,85)
	DO I=1,NE
		WRITE(30,90) I,(NOD(I,J),J=1,NDS),MAT(I),JAD(I),JRD(I)
	END DO
	WRITE(30,95)
	DO I=1,NX
		WRITE(30,100) I,(COR(I,J),J=1,NDM),(ID(I,K),K=1,NFR)
	END DO
	IF(NF>0)THEN
		WRITE(30,105)
		DO I=1,NF
			WRITE(30,110) I,(KSF(I,J),J=1,NFS+3)
		END DO
	END IF
	IF(NQ>0)THEN
		WRITE(30,115)
		DO I=1,NQ
			WRITE(30,120) I,(KSQ(I,J),J=1,NFS+2),DSQ(I)
		END DO
	END IF
	WRITE(30,125)
	DO I=1,NE
		WRITE(30,130) I,(STR(I,K),K=1,6)
	END DO
	WRITE(30,135)
	DO I=1,NX
		WRITE(30,140) I,(GG(I,L),L=1,NFR)
	END DO
5	FORMAT('CONTROL INFORMATION ICON(10):')
10	FORMAT(10I5)
15	FORMAT(/,'NP=',5X,'NM=',5X,'NE=',5X,'NX=',5X,'NF=',5X,'NQ=')
20	FORMAT(2I5,2I10,2I8)
25	FORMAT(/,'STAGE CONTROL INFORMATION LQ(10):')
30	FORMAT(11I5,2F10.3,4X,2F10.3)
35	FORMAT(/,'MATERIAL CONTROL INFORMATION KGSS(7):')
40	FORMAT(I2,4X,7I5)
45	FORMAT(/,'STATIC MATERIAL PARAMETERS CS(15):')
50	FORMAT(I5,15F17.7)
55	FORMAT(/,'PERMEABILITY PARAMETERS CK(15):')
60	FORMAT(I5,15F17.7)
65	FORMAT(/,'WETTING DEFORMATION PARAMETERS CW(15):')
70	FORMAT(I5,15F17.7)
75	FORMAT(/,'RHEOLOGICAL DEFORMATION PARAMETERS CR(15):')
80	FORMAT(I5,15F17.7)
85	FORMAT(/,'ELEMENTS INFORMATION NOD(8),MAT(1):')
90	FORMAT(I5,8I10,3I5)
95	FORMAT(/,'NODS INFORMATION COR(3),ID(4),GG(4):')
100	FORMAT(I5,3F12.3,4I5)
105	FORMAT(/,'SURFACE LOAD BOUNDARY KSF(NSF+3):')
110	FORMAT(I5,4X,I5,4I10,2I4)
115	FORMAT(/,'SURFACE FLUX BOUNDARY KSQ(NSF+3):')
120	FORMAT(I5,4X,I5,4I10,I4,4X,F12.3)
125	FORMAT(/,'INITIAL STRESS STATE STR(6):')
130	FORMAT(I5,6F15.4)
135	FORMAT(/,'INITIAL NODAL INFORMATION GG(4):')
140	FORMAT(I5,4F15.4)
	RETURN
END SUBROUTINE DATAOUT


SUBROUTINE SVRST(IP)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,I,J,K
	ntep=ntep+1
	WRITE(40,'(2i8)')  ntep,IP
	WRITE(50,'(2i8)')  ntep,IP
	WRITE(60,'(2i8)')  ntep,IP
	WRITE(70,'(2i8)')  ntep,IP
	WRITE(80,'(2i8)')  ntep,IP
	WRITE(90,'(2i8)')  ntep,IP
	WRITE(100,'(2i8)') ntep,IP
	WRITE(217,'(2i8)') ntep,IP
!	WRITE(218,'(2i8,2f18.9)') ntep,IP,totaltime,hwu(ip,2)
	DO I=1,NE
        WRITE(217,'(2(I5,4X),F19.9)')I,MAT(I),SEEP(I,1)
	END DO
	DO I=1,NX
		WRITE(50,5) I,(DIS(I,J),J=1,NDM),TWH(I)
	END DO
	DO I=1,NE
		WRITE(60,10) I,(STS(I,K),K=1,12)
	END DO
	DO I=1,NE
		WRITE(70,15) I,(STN(I,J),J=1,12)
	END DO
	DO I=1,NE
		WRITE(80,20) I,IST(I),LDG(I),(EST(I,K),K=1,4),(PLS(I,K),K=1,2)
	END DO
	DO I=1,NE
		WRITE(90,25) I,(STS(I,J),J=7,9),PWS(I),((VCS(I,J,K),J=1,NDM),K=1,NDM)
	END DO
	DO I=1,NE
		WRITE(100,30) I,(STN(I,J),J=7,9),((VCN(I,J,K),J=1,NDM),K=1,NDM)
	END DO
5	FORMAT(I8,4X,3F15.4,F18.4)
10	FORMAT(I8,4X,12F15.4)
15	FORMAT(I8,4X,12F15.4)
20	FORMAT(I8,4X,I4,4X,I4,4X,4F17.6,2F12.3)
25	FORMAT(I8,4X,4F17.6,4X,9F15.4)
30	FORMAT(I8,4X,3F17.6,4X,9F15.4)
	RETURN
END SUBROUTINE SVRST


SUBROUTINE CONFIG(N,A,B)
	IMPLICIT NONE
	INTEGER::N
	REAL(KIND=8),DIMENSION(N)::A,B
	IF(N==2)THEN
		A(1)=1.000000000000000
		A(2)=1.000000000000000
		B(1)=-0.577350269189626
		B(2)=+0.577350269189626
	ELSE IF(N==3)THEN
		A(1)=0.555555555555556
		A(2)=0.888888888888889
		A(3)=0.555555555555556
		B(1)=-0.774596669241483
		B(2)=+0.000000000000000
		B(3)=+0.774596669241483
	ELSE IF(N==4)THEN
		A(1)=0.347854845137454
		A(2)=0.652145154862546
		A(3)=0.652145154862546
		A(4)=0.347854845137454
		B(1)=-0.861136311594053
		B(2)=-0.339981043584856
		B(3)=+0.339981043584856
		B(4)=+0.861136311594053
	ELSE IF(N==5)THEN
		A(1)=0.236926885056189
		A(2)=0.478628670499366
		A(3)=0.568888888888889
		A(4)=0.478628670499366
		A(5)=0.236926885056189
		B(1)=-0.906179845938664
		B(2)=-0.538469310105683
		B(3)=+0.000000000000000
		B(4)=+0.538469310105683
		B(5)=+0.906179845938664
	ELSE IF(N==6)THEN
		A(1)=0.171324492379170
		A(2)=0.360761573048139
		A(3)=0.467913934572691
		A(4)=0.467913934572691
		A(5)=0.360761573048139
		A(6)=0.171324492379170
		B(1)=-0.932469514203152
		B(2)=-0.661209386466265
		B(3)=-0.238619186083197
		B(4)=+0.238619186083197
		B(5)=+0.661209386466265
		B(6)=+0.932469514203152
	ELSE IF(N==7)THEN
		A(1)=0.129484966168870
		A(2)=0.279705391489277
		A(3)=0.381830050505119
		A(4)=0.417959183673469
		A(5)=0.381830050505119
		A(6)=0.279705391489277
		A(7)=0.129484966168870
		B(1)=-0.949107912342759
		B(2)=-0.741531185599394
		B(3)=-0.405845151377397
		B(4)=+0.000000000000000
		B(5)=+0.405845151377397
		B(6)=+0.741531185599394
		B(7)=+0.949107912342759
	END IF
	RETURN
END SUBROUTINE CONFIG


SUBROUTINE CHECI(IE,IR)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,IR,J,K,L,NI,NJ,NK,s
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB
	REAL(KIND=8)::GX,GY,GZ,DET
	IR=0
	s=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	DO NI=1,NGS
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GZ=WXYZ(NK)
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				IF(DET<=1.0E-12)THEN
					IR=1
					WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE CHECI(IE,IR):',IE,mat(ie)
					WRITE(40,*) 'FOR DET<=1.0E-12',det
				END IF
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE CHECI


SUBROUTINE CHECJ(IE,IR)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::DET
	INTEGER::IE,IR,J,K,L,JJ,KK
	DET=0.0
	IR=0
	DO J=1,NFS
		K=J+NFS
		JJ=NOD(IE,J)
		KK=NOD(IE,K)
		DO L=1,NDM
			DET=DET+ABS(COR(JJ,L)-COR(KK,L))
		END DO
	END DO
	DET=DET/REAL(NFS)
	IF(DET>0.5)THEN
		IR=0
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE CHECJ(IE,IR):',IE,mat(ie)
        WRITE(40,'(8(i6,4x))')(nod(ie,j),j=1,nds)
		WRITE(40,*) 'FOR DET>=5.0E-1',det
	END IF
	RETURN
END SUBROUTINE CHECJ


SUBROUTINE CHECK(IE,IR)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,IR,J,K,L,JJ,KK
	IR=0
	DO J=1,NFS
		K=J+1
		IF(K==NFS+1) K=1
		JJ=NOD(IE,J)
		KK=NOD(IE,K)
		IF(JJ/=KK) IR=1
	END DO
	DO L=1,NFS
		J=NFS+L
		K=J+1
		IF(K==NDS+1) K=NFS+1
		JJ=NOD(IE,J)
		KK=NOD(IE,K)
		IF(JJ/=KK) IR=1
	END DO
	IF(IR==1)THEN
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE CHECK(IE,IR):',IE
		WRITE(40,*) 'FOR ELEMENT NODES CONFLICTS!'
	END IF
	RETURN
END SUBROUTINE CHECK


SUBROUTINE SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
	USE COMM
	IMPLICIT NONE
	INTEGER::I,s
	REAL(KIND=8)::GX,GY,GZ
	REAL(KIND=8),DIMENSION(NDS)::SH,A,B,C
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV
!	DATA A/-1.0,-1.0,-1.0,-1.0,+1.0,+1.0,+1.0,+1.0/
!	DATA B/-1.0,+1.0,+1.0,-1.0,-1.0,+1.0,+1.0,-1.0/
!	DATA C/-1.0,-1.0,+1.0,+1.0,-1.0,-1.0,+1.0,+1.0/
	if(s.eq.1.or.s.eq.0)then
    a(1)=-1;a(2)=-1;a(3)=-1;a(4)=-1;a(5)=1;a(6)=1;a(7)=1;a(8)=1
    b(1)=1;b(2)=-1;b(3)=-1;b(4)=1;b(5)=1;b(6)=-1;b(7)=-1;b(8)=1
    c(1)=-1;c(2)=-1;c(3)=1;c(4)=1;c(5)=-1;c(6)=-1;c(7)=1;c(8)=1
    elseif(s.eq.2)then
    a(1)=-1;a(2)=1;a(3)=1;a(4)=-1;a(5)=-1;a(6)=1;a(7)=1;a(8)=-1
    b(1)=-1;b(2)=-1;b(3)=-1;b(4)=-1;b(5)=1;b(6)=1;b(7)=1;b(8)=1
    c(1)=-1;c(2)=-1;c(3)=1;c(4)=1;c(5)=-1;c(6)=-1;c(7)=1;c(8)=1
    elseif(s.eq.3)then
    a(1)=1;a(2)=-1;a(3)=-1;a(4)=1;a(5)=1;a(6)=-1;a(7)=-1;a(8)=1
    b(1)=-1;b(2)=-1;b(3)=1;b(4)=1;b(5)=-1;b(6)=-1;b(7)=1;b(8)=1
    c(1)=-1;c(2)=-1;c(3)=-1;c(4)=-1;c(5)=1;c(6)=1;c(7)=1;c(8)=1
	endif
	DO I=1,NDS
		SH(I)=0.125*(1.0+A(I)*GX)*(1.0+B(I)*GY)*(1.0+C(I)*GZ)
	END DO
	DO I=1,NDS
		DERIV(1,I)=0.125*(1.0+B(I)*GY)*(1.0+C(I)*GZ)*A(I)
		DERIV(2,I)=0.125*(1.0+A(I)*GX)*(1.0+C(I)*GZ)*B(I)
		DERIV(3,I)=0.125*(1.0+A(I)*GX)*(1.0+B(I)*GY)*C(I)
	END DO
	RETURN
END SUBROUTINE SHAPEFUNC


SUBROUTINE GRVT(IE,GV)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,KT,NI,NJ,NK,s
	REAL(KIND=8)::GHX,GHY,GHZ,GX,GY,GZ,DET,GMA,RR
	REAL(KIND=8),DIMENSION(NDS)::SH,GV
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB
	KT=MAT(IE)
	s=jplane(ie)
	IF(KGSS(KT,1)<0) RETURN
	GMA=CS(KT,15)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	DO NI=1,NGS
		GHX=WGH(NI)
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GHY=WGH(NJ)
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GHZ=WGH(NK)
				GZ=WXYZ(NK)
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				RR=GHX*GHY*GHZ*DET
				DO J=1,NDS
					GV(J)=GV(J)-RR*GMA*SH(J)
				END DO
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE GRVT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SRFP1(IP,ISF,SF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,ISF,IE,IS,IA,I,J,K,L,KK,KL,POI
	REAL(KIND=8)::VX,VY,VZ,HMX,HMN,ZG,DH,AA,RR
	REAL(KIND=8),DIMENSION(NDS*NDM)::SF
	REAL(KIND=8),DIMENSION(NGF,NGF,NDM)::LMN
	REAL(KIND=8),DIMENSION(NGF,NGF,NDS)::SHP
	REAL(KIND=8),DIMENSION(NGF,NGF)::WGS,EVT
	IE=KSP(ISF,1)
	IS=KSP(ISF,6)
	IA=KSP(ISF,7)
!	IF(IP==13)HWU(IP,1)=682
	IF(IA<0)THEN
		HMX=MAX(HWU(IP,1),HWU(IP,2))
		HMN=MIN(HWU(IP,1),HWU(IP,2))
		AA=1.0
		IF(HWU(IP,2)<HWU(IP,1)) AA=-1.0
	ELSE IF(IA>0)THEN
		HMX=MAX(HWD(IP,1),HWD(IP,2))
		HMN=MIN(HWD(IP,1),HWD(IP,2))
		AA=1.0
		IF(HWD(IP,2)<HWD(IP,1)) AA=-1.0
	END IF
	POI=ABS(IS)
	IF(POI==1)THEN
		IF(IS<0) VX=-1.0
		IF(IS>0) VX=+1.0
		CALL SRFI(IE,VX,LMN,SHP,WGS,EVT)
	ELSE IF(POI==2)THEN
		IF(IS<0) VY=-1.0
		IF(IS>0) VY=+1.0
		CALL SRFJ(IE,VY,LMN,SHP,WGS,EVT)
	ELSE IF(POI==3)THEN
		IF(IS<0) VZ=-1.0
		IF(IS>0) VZ=+1.0
		CALL SRFK(IE,VZ,LMN,SHP,WGS,EVT)
	END IF
	DO I=1,NGF
		DO J=1,NGF
			ZG=EVT(I,J)
			IF(ZG<=HMN)THEN
				DH=HMX-HMN
			ELSE IF(ZG<=HMX)THEN
				DH=HMX-ZG
			ELSE
				DH=0.0
			END IF
			RR=GMAW*DH*WGS(I,J)
			DO K=1,NDS
				KK=(K-1)*NDM
				DO L=1,NDM
					KL=KK+L
					SF(KL)=SF(KL)+SHP(I,J,K)*LMN(I,J,L)*RR
				END DO
			END DO
		END DO
	END DO
	DO KL=1,NDS*NDM
		SF(KL)=-AA*SF(KL)
	END DO
	RETURN
END SUBROUTINE SRFP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!增加计算心墙上游水压力的子程序

SUBROUTINE SRFP(IP,ISF,SF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,ISF,IE,IS,IA,I,J,K,L,KK,KL,POI
	REAL(KIND=8)::VX,VY,VZ,HMX,HMN,ZG,DH,AA,RR
	REAL(KIND=8),DIMENSION(NDS*NDM)::SF
	REAL(KIND=8),DIMENSION(NGF,NGF,NDM)::LMN
	REAL(KIND=8),DIMENSION(NGF,NGF,NDS)::SHP
	REAL(KIND=8),DIMENSION(NGF,NGF)::WGS,EVT
	IE=KSF(ISF,1)
	IS=KSF(ISF,6)
	IA=KSF(ISF,7)
	IF(IA<0)THEN
		HMX=MAX(HWU(IP,1),HWU(IP,2))
		HMN=MIN(HWU(IP,1),HWU(IP,2))
		AA=1.0
		IF(HWU(IP,2)<HWU(IP,1)) AA=-1.0
	ELSE IF(IA>0)THEN
		HMX=MAX(HWD(IP,1),HWD(IP,2))
		HMN=MIN(HWD(IP,1),HWD(IP,2))
		AA=1.0
		IF(HWD(IP,2)<HWD(IP,1)) AA=-1.0
	END IF
	POI=ABS(IS)
	IF(POI==1)THEN
		IF(IS<0) VX=-1.0
		IF(IS>0) VX=+1.0
		CALL SRFI(IE,VX,LMN,SHP,WGS,EVT)
	ELSE IF(POI==2)THEN
		IF(IS<0) VY=-1.0
		IF(IS>0) VY=+1.0
		CALL SRFJ(IE,VY,LMN,SHP,WGS,EVT)
	ELSE IF(POI==3)THEN
		IF(IS<0) VZ=-1.0
		IF(IS>0) VZ=+1.0
		CALL SRFK(IE,VZ,LMN,SHP,WGS,EVT)
	END IF
	DO I=1,NGF
		DO J=1,NGF
			ZG=EVT(I,J)
			IF(ZG<=HMN)THEN
				DH=HMX-HMN
			ELSE IF(ZG<=HMX)THEN
				DH=HMX-ZG
			ELSE
				DH=0.0
			END IF
			RR=GMAW*DH*WGS(I,J)
			DO K=1,NDS
				KK=(K-1)*NDM
				DO L=1,NDM
					KL=KK+L
					SF(KL)=SF(KL)+SHP(I,J,K)*LMN(I,J,L)*RR
				END DO
			END DO
		END DO
	END DO
	DO KL=1,NDS*NDM
		SF(KL)=-AA*SF(KL)
	END DO
	RETURN
END SUBROUTINE SRFP


SUBROUTINE SRFQ(IP,ISQ,QV)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,ISQ,IE,IS,I,J,K,POI
	REAL(KIND=8)::VX,VY,VZ,QA,R1,R2,R3,RR,AA
	REAL(KIND=8),DIMENSION(NDS)::QV
	REAL(KIND=8),DIMENSION(NGF,NGF,NDM)::LMN
	REAL(KIND=8),DIMENSION(NGF,NGF,NDS)::SHP
	REAL(KIND=8),DIMENSION(NGF,NGF)::WGS,EVT
	IE=KSQ(ISQ,1)
	IS=KSQ(ISQ,6)
	QA=DSQ(ISQ)
	POI=ABS(IS)
	IF(POI==1)THEN
		IF(IS<0) VX=-1.0
		IF(IS>0) VX=+1.0
		CALL SRFI(IE,VX,LMN,SHP,WGS,EVT)
	ELSE IF(POI==2)THEN
		IF(IS<0) VY=-1.0
		IF(IS>0) VY=+1.0
		CALL SRFJ(IE,VY,LMN,SHP,WGS,EVT)
	ELSE IF(POI==3)THEN
		IF(IS<0) VZ=-1.0
		IF(IS>0) VZ=+1.0
		CALL SRFK(IE,VZ,LMN,SHP,WGS,EVT)
	END IF
	DO I=1,NGF
		DO J=1,NGF
			R1=LMN(I,J,1)
			R2=LMN(I,J,2)
			R3=LMN(I,J,3)
			RR=SQRT(R1**2+R2**2+R3**2)
			AA=GMAW*WGS(I,J)*QA*RR
			DO K=1,NDS
				QV(K)=QV(K)+SHP(I,J,K)*AA
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SRFQ


SUBROUTINE SRFI(IE,VX,LMN,SHP,WGS,EVT)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,NJ,NK,s
	REAL(KIND=8)::VX,VY,VZ,VHY,VHZ,SGN
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB
	REAL(KIND=8),DIMENSION(NGF,NGF,NDM)::LMN
	REAL(KIND=8),DIMENSION(NGF,NGF,NDS)::SHP
	REAL(KIND=8),DIMENSION(NGF,NGF)::WGS,EVT
	s=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	IF(VX>0.0) SGN=1.0
	IF(VX<0.0) SGN=-1.0
	DO NJ=1,NGF
		VY=VXYZ(NJ)
		VHY=VGH(NJ)
		DO NK=1,NGF
			VZ=VXYZ(NK)
			VHZ=VGH(NK)
			CALL SHAPEFUNC(VX,VY,VZ,SH,DERIV,s)
			CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
			LMN(NJ,NK,1)=SGN*(JCB(2,2)*JCB(3,3)-JCB(2,3)*JCB(3,2))
			LMN(NJ,NK,2)=SGN*(JCB(2,3)*JCB(3,1)-JCB(2,1)*JCB(3,3))
			LMN(NJ,NK,3)=SGN*(JCB(2,1)*JCB(3,2)-JCB(2,2)*JCB(3,1))
			FORALL(L=1:NDS) SHP(NJ,NK,L)=SH(L)
			WGS(NJ,NK)=VHY*VHZ
			EVT(NJ,NK)=0.0
			DO L=1,NDS
				EVT(NJ,NK)=EVT(NJ,NK)+SH(L)*XYZ(L,NDM)
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SRFI


SUBROUTINE SRFJ(IE,VY,LMN,SHP,WGS,EVT)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,NI,NK,s
	REAL(KIND=8)::VX,VY,VZ,VHX,VHZ,SGN
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB
	REAL(KIND=8),DIMENSION(NGF,NGF,NDM)::LMN
	REAL(KIND=8),DIMENSION(NGF,NGF,NDS)::SHP
	REAL(KIND=8),DIMENSION(NGF,NGF)::WGS,EVT
	s=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	IF(VY>0.0) SGN=1.0
	IF(VY<0.0) SGN=-1.0
	DO NI=1,NGF
		VX=VXYZ(NI)
		VHX=VGH(NI)
		DO NK=1,NGF
			VZ=VXYZ(NK)
			VHZ=VGH(NK)
			CALL SHAPEFUNC(VX,VY,VZ,SH,DERIV,s)
			CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
			LMN(NI,NK,1)=SGN*(JCB(1,3)*JCB(3,2)-JCB(1,2)*JCB(3,3))
			LMN(NI,NK,2)=SGN*(JCB(1,1)*JCB(3,3)-JCB(1,3)*JCB(3,1))
			LMN(NI,NK,3)=SGN*(JCB(1,2)*JCB(3,1)-JCB(1,1)*JCB(3,2))
			FORALL(L=1:NDS) SHP(NI,NK,L)=SH(L)
			WGS(NI,NK)=VHX*VHZ
			EVT(NI,NK)=0.0
			DO L=1,NDS
				EVT(NI,NK)=EVT(NI,NK)+SH(L)*XYZ(L,NDM)
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SRFJ


SUBROUTINE SRFK(IE,VZ,LMN,SHP,WGS,EVT)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,NI,NJ,s
	REAL(KIND=8)::VX,VY,VZ,VHX,VHY,SGN
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB
	REAL(KIND=8),DIMENSION(NGF,NGF,NDM)::LMN
	REAL(KIND=8),DIMENSION(NGF,NGF,NDS)::SHP
	REAL(KIND=8),DIMENSION(NGF,NGF)::WGS,EVT
    s=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	IF(VZ>0.0) SGN=1.0
	IF(VZ<0.0) SGN=-1.0
	DO NI=1,NGF
		VX=VXYZ(NI)
		VHX=VGH(NI)
		DO NJ=1,NGF
			VY=VXYZ(NJ)
			VHY=VGH(NJ)
			CALL SHAPEFUNC(VX,VY,VZ,SH,DERIV,s)
			CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
			LMN(NI,NJ,1)=SGN*(JCB(1,2)*JCB(2,3)-JCB(1,3)*JCB(2,2))
			LMN(NI,NJ,2)=SGN*(JCB(1,3)*JCB(2,1)-JCB(1,1)*JCB(2,3))
			LMN(NI,NJ,3)=SGN*(JCB(1,1)*JCB(2,2)-JCB(1,2)*JCB(2,1))
			FORALL(L=1:NDS) SHP(NI,NJ,L)=SH(L)
			WGS(NI,NJ)=VHX*VHY
			EVT(NI,NJ)=0.0
			DO L=1,NDS
				EVT(NI,NJ)=EVT(NI,NJ)+SH(L)*XYZ(L,NDM)
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SRFK


SUBROUTINE VOLQ(IE,QV)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LP,J,K
	REAL(KIND=8),DIMENSION(NDS)::QV,TH
	REAL(KIND=8),DIMENSION(NDS,NDS)::KH
	LM=MAT(IE)
	LP=KGSS(LM,2)
	IF(LP==0) RETURN
	DO J=1,NDS
		K=NOD(IE,J)
		TH(J)=TWH(K)
	END DO
	DO J=1,NDS
		DO K=1,NDS
			KH(J,K)=0.0
		END DO
	END DO
	CALL SDKH(IE,KH)
	DO J=1,NDS
		DO K=1,NDS
			QV(J)=QV(J)+KH(J,K)*TH(K)
		END DO
	END DO
	RETURN
END SUBROUTINE VOLQ


SUBROUTINE INITCOND
	USE COMM
	IMPLICIT NONE
	INTEGER::I,J,K,LM,LT,S,ip,ix,le,ie,kt,me,KTT
	REAL(KIND=8)::ZC,ZH,I1,I2,I3,J2,J3,PP,QQ,QP,htg,density,hei,stsver,ht0
	REAL(KIND=8),DIMENSION(9)::RST
	REAL(KIND=8),DIMENSION(6)::RST0,rt
	REAL(KIND=8),DIMENSION(NDM,NDM)::SS,VE
!   initial stress
goto 111
do le=1,ne
    !forall(j=1:6) str(le,j)=0.0
end do
    do ip=1,np
	call mqq(ip)
	call mpp(ip)
	if(lq(ip,1).ne.1)cycle
    me=lq(ip,2)
	htg=conh(ip)
	ht0=hsur
	do ie=1,me
    le=lus(ie)
	kt=mat(le)
    KTT=KGSS(KT,1)
	density=cs(kt,15)
    if(ip.eq.1.and.ist(le).eq.1)then
    if(ht0.eq.0.0)then
      if(KTT<0)then
	str(le,1)=0
	str(le,2)=0
	str(le,3)=-str(le,3)
	str(le,4)=0
	str(le,5)=0
	str(le,6)=0
	  endif
	  cycle
	end if
	hei=0.0
	do j=1,nds
     hei=hei+cor(nod(le,j),ndm)/real(nds)
	end do
    stsver=(ht0-hei)*density
	if(stsver.le.100)stsver=100
	if(kTT>=0)then
	str(le,1)=stsver*0.5
	str(le,2)=str(le,1)
	str(le,3)=stsver
	str(le,4)=0
	str(le,5)=0
	str(le,6)=0
    else
	str(le,1)=0
	str(le,2)=0
	str(le,3)=-stsver
	str(le,4)=0
	str(le,5)=0
	str(le,6)=0
	end if
	end if

!	print*,ist(le)
    if(ist(le).eq.2)then
	hei=0.0
	do j=1,nds
     hei=hei+cor(nod(le,j),ndm)/real(nds)
	end do
    stsver=(htg-hei)*density
	if(stsver.le.100)stsver=100
    if(ktT>=0)then
	str(le,1)=stsver*0.5
	str(le,2)=str(le,1)
	str(le,3)=stsver
	str(le,4)=0
	str(le,5)=0
	str(le,6)=0
    else
	str(le,1)=0
	str(le,2)=0
	str(le,3)=-stsver
	str(le,4)=0
	str(le,5)=0
	str(le,6)=0
	end if
	end if
    end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 continue
	DO I=1,NX
		DO J=1,NDM
			DIS(I,J)=GG(I,J)
		END DO
		TWH(I)=GG(I,NFR)
	END DO
	IF(ICON(2)==0)THEN
		DO I=1,NX
			ID(I,NFR)=0
		END DO
		ELSE
		  DO J=1,NE
		  LM=MAT(J)
		  LT=KGSS(LM,2)
		  IF(LT==0)THEN
		    DO S=1,NDS
			  K=NOD(J,S)
			  ID(K,NFR)=0
			  END DO
          END IF
         END DO
	END IF
	DO I=1,NE
		LDG(I)=+1
		DO J=1,2
			RHS(I,J)=0.0
			PLS(I,J)=1.0
			YLD(I,J)=0.0
			YLDM(I,J)=0.0
		END DO
		DO J=1,4
			EST(I,J)=0.0
		END DO
		DO J=1,12
			STS(I,J)=0.0
			STN(I,J)=0.0
		END DO
		DO J=1,6
			STS(I,J)=STR(I,J)
			rhstn(i,j)=0.0
		END DO
	END DO
	DO I=1,NE
		PWS(I)=0.0
		LM=MAT(I)
		LT=KGSS(LM,2)
		IF(LM==0) CYCLE
		ZC=0.0
		ZH=0.0
		DO J=1,NDS
			K=NOD(I,J)
			ZC=ZC+COR(K,NDM)
			ZH=ZH+TWH(K)
		END DO
		PWS(I)=GMAW*(ZH-ZC)/REAL(NDS)
	END DO
    DO I=1,NE
    LM=MAT(I)
	SEEP(I,1)=CK(LM,1)*864
    SEEP(I,2)=CK(LM,5)*864
	SEEP(I,3)=CK(LM,9)*864
	END DO
	DO I=1,NE
		LM=MAT(I)
		LT=KGSS(LM,1)
		IF(LT>=0)THEN
			SS(1,1)=STS(I,1)
			SS(2,2)=STS(I,2)
			SS(3,3)=STS(I,3)
			SS(1,2)=STS(I,4)
			SS(1,3)=STS(I,6)
			SS(2,1)=STS(I,4)
			SS(2,3)=STS(I,5)
			SS(3,1)=STS(I,6)
			SS(3,2)=STS(I,5)
			CALL MNSTS(NDM,SS,VE)
			DO J=1,NDM
				DO K=1,NDM
					VCS(I,J,K)=VE(J,K)
				END DO
			END DO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do j=1,6
     rt(j)=sts(i,j)
	end do
!	print*,i,ss(3,3)
   ! call stsmod(i,ss,rt)
	do j=1,6
     sts(i,j)=rt(j)
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			STS(I,7)=SS(1,1)
			STS(I,8)=SS(2,2)
			STS(I,9)=SS(3,3)
			DO J=1,9
				RST(J)=STS(I,J)
			END DO
			DO K=1,6
				RST0(K)=STR(I,K)
			END DO
			CALL INVARS(RST,I1,I2,I3,J2,J3,PP,QQ)
			STS(I,10)=PP
			STS(I,11)=QQ
        !	print*,'4',i
		!	print*,i,-3*sqrt(3.0)/2.0*j3/j2**(1.5)
			if(abs(-3*sqrt(3.0)/2.0*j3/j2**(1.5))>=1.0)then
            print*,"error occur in the asin:",i,mat(i)
			 stop
			end if
			thata(i)=asin(-3*sqrt(3.0)/2.0*j3/j2**(1.5))/3.0
		!	print*,i,-3*sqrt(3.0)/2.0*j3/j2**(1.5)
		!	print*,'5',i
			CALL STSLVL(LM,LT,RST,RST0,QP)
			STS(I,12)=QP
			CALL YLDMCALC(I)
		ELSE
	!	print*,i
			DO J=1,NDM
				DO K=1,NDM
					SS(J,K)=0.0
				END DO
			END DO
			IF(LT==-1) CALL MATRTRANI(I,SS,PP)
			IF(LT==-2) CALL MATRTRANI(I,SS,PP)
			IF(LT==-3) CALL MATRTRANI(I,SS,PP)
			IF(LT==-4) CALL MATRTRANJ(I,SS,PP)
			EST(I,4)=PP
		!	write(40,*)i,pp
			DO J=1,NDM
			!if(i.eq.17815)then
             ! write(40,'(3f18.9)')(ss(j,k),k=1,ndm)
		!	end if
				DO K=1,NDM
					DMT(I,NDM+J,NDM+K)=SS(J,K)
				END DO
			END DO
		END IF
	END DO
!pause
	RETURN
END SUBROUTINE INITCOND


SUBROUTINE EVSNCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT
	REAL(KIND=8)::A,B,C,D
	KM=MAT(IE)
	KT=KGSS(KM,1)
	A=0.0
	B=0.0
	C=0.0
	D=0.0
	SELECT CASE(KT)
		CASE(-4)
			CALL KSKN4(IE,KM,A,B,C,D)
		CASE(-3)
			CALL KSKN3(IE,KM,A,B,C,D)
		CASE(-2)
			CALL KSKN2(IE,KM,A,B,C,D)
		CASE(-1)
			CALL KSKN1(IE,KM,A,B,C,D)
		CASE(0)
			CALL EVKG0(IE,KM,A,B,C,D)
		CASE(1)
			CALL EVKG1(IE,KM,A,B,C,D)
		CASE(2)
			CALL EVKG2(IE,KM,A,B,C,D)
		CASE(3)
			CALL EVKG3(IE,KM,A,B,C,D)
		CASE(11)
			CALL EVKG11(IE,KM,A,B,C,D)
		CASE(12)
			CALL EVKG12(IE,KM,A,B,C,D)
		CASE(13)
			CALL EVKG13(IE,KM,A,B,C,D)
		CASE(14)
			CALL EVKG14(IE,KM,A,B,C,D)
		CASE(15)
			CALL EVKG15(IE,KM,A,B,C,D)
		CASE DEFAULT
			WRITE(40,*) 'ERRORS MAY OCCUR IN SUBROUTINE EVSNCALC(IE)!'
			WRITE(40,*) 'IE=',IE
			STOP
	END SELECT
	EST(IE,1)=A
	EST(IE,2)=B
	EST(IE,3)=C
	EST(IE,4)=D
	RETURN
END SUBROUTINE EVSNCALC


SUBROUTINE EVKG0(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D
	A=CS(KM,2)
	B=CS(KM,4)
	IF(B>0.49) B=0.49
	IF(B<0.01) B=0.01
	C=A/(1.0-2.0*B)/3.0
	D=A/(1.0+B)/2.0
	RETURN
END SUBROUTINE EVKG0


SUBROUTINE EVKG1(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,SG3,SL,RW1,RW2,s,phi,dphi
	REAL(KIND=8)::AR,AK,AN,AG,AF,AD
	AR=CS(KM,1)
	AK=CS(KM,2)
	AN=CS(KM,3)
	AG=CS(KM,4)
	AF=CS(KM,5)
	AD=CS(KM,6)
    phi=cs(km,12)
	dphi=cs(km,13)
	SG3=STS(IE,10)
	sl=sts(ie,12)
	IF(SG3<PAR) SG3=PAR
	IF(SL>0.90) SL=0.90
	IF(LDG(IE)==-1)THEN
		AK=CS(KM,7)
		AN=CS(KM,8)
		A=AK*PAR*(SG3/PAR)**AN
	ELSE
		A=ak*par*(sg3/par)**an*(1.0-AR*SL)**2
	END IF
	if(sg3<par)then
    b=0.3
	goto 20
	end if
	RW1=(1.0-AR*SL)*AK*PAR*(SG3/PAR)**AN
	RW2=AD*(STS(IE,7)-STS(IE,9))/RW1
	B=(AG-AF*LOG10(SG3/PAR))/(1.0-RW2)**2
	IF(B>0.49) B=0.49
	IF(B<0.01) B=0.01
    20 continue
	C=A/(1.0-2.0*B)/3.0
	D=A/(1.0+B)/2.0
	RETURN
END SUBROUTINE EVKG1


SUBROUTINE EVKG2(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,SG3,SL,RW1,RW2
	REAL(KIND=8)::AR,AK,AN,AB,AM
	AR=CS(KM,1)
	AK=CS(KM,2)
	AN=CS(KM,3)
	AB=CS(KM,4)
	AM=CS(KM,5)
	SG3=STS(IE,10)
	SL=STS(IE,12)
	IF(SG3<PAR) SG3=PAR
	IF(SL>0.90) SL=0.90
	IF(LDG(IE)==-1)THEN
		AK=CS(KM,7)
		AN=CS(KM,8)
		A=AK*PAR*(SG3/PAR)**AN
	ELSE
		a=ak*par*(sg3/par)**an*(1.0-AR*SL)**2
	END IF
   if(sg3<par)then
   b=0.3
   goto 20
   end if
	RW2=AB*PAR*(SG3/PAR)**AM
!	B=0.5-RW1*(1.0-AR*SL)/RW2/6.0
	B=0.5-A/RW2/6.0
	IF(B>0.49) B=0.49
	IF(B<0.01) B=0.01
20 continue
	C=A/(1.0-2.0*B)/3.0
	D=A/(1.0+B)/2.0
	RETURN
END SUBROUTINE EVKG2


SUBROUTINE EVKG3(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,PP,QQ,K,G
	REAL(KIND=8)::KI,AK,GI,AG,BG
	KI=CS(KM,1)
	AK=CS(KM,2)
	GI=CS(KM,3)
	AG=CS(KM,4)
	BG=CS(KM,5)
	PP=STS(IE,10)
	QQ=STS(IE,11)
	IF(PP<PAR) PP=PAR
	K=KI+AK*PP
	G=GI+AG*PP+BG*QQ
	A=9.0*K*G/(3.0*K+G)
	B=(3.0*K-2.0*G)/(6.0*K+2.0*G)
	IF(B>0.49) B=0.49
	IF(B<0.01) B=0.01
	C=A/(1.0-2.0*B)/3.0
	D=A/(1.0+B)/2.0
	RETURN
END SUBROUTINE EVKG3


SUBROUTINE EVKG11(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,SG3,SL,RW0,RW1,RW2,RW3,WR,WS,PP,QQ,YT
	REAL(KIND=8)::AR,AK,AN,BC,BN,BR,ET,VV,UT,K,G
	WR=2.0
	WS=2.0
	AR=CS(KM,1)
	AK=CS(KM,2)
	AN=CS(KM,3)
	VV=CS(KM,6)
	BC=CS(KM,9)
	BN=CS(KM,10)
	BR=CS(KM,11)
	SG3=STS(IE,10)
	PP=STS(IE,10)
	QQ=STS(IE,11)
	SL=STS(IE,12)
	YT=QQ/PP
	IF(SG3<PAR) SG3=PAR
	IF(SL>0.90) SL=0.90
	RW0=AK*PAR*((SG3/PAR)**AN)
	AK=CS(KM,7)
	AN=CS(KM,8)
	A=AK*PAR*((SG3/PAR)**AN)
	B=VV
	K=A/(3.0*(1.0-2.0*B))
	G=A/(2.0*(1.0+B))
	ET=RW0*(1.0-AR*SL)**2
	RW1=BC*(SG3/PAR)**BN
	RW2=(AR*SL*(1.0-BR))/((STS(IE,7)-STS(IE,9))*BR)
	RW3=1.0-AR*SL*(1.0-BR)/((1.0-AR*SL)*BR)
	UT=2.0*RW0*RW1*RW2*RW3
	RW0=WS+WR*WR*YT*YT
	RW1=(9.0/ET)-(3.0*UT/ET)-(3.0/G)
	RW2=(3.0*UT/ET)-(1.0/K)
	C=(YT*RW1+3.0*WS*RW2)/(3.0*(1.0+3.0*WR*WR*YT)*RW0)
	D=(RW1-3.0*WR*WR*YT*RW2)/(3.0*(3.0*WS-YT)*RW0)
	C=0.25*C/(PP*PP)
	D=D*PP*PP*QQ*QQ/(QQ**(2.0*WS))
	IF(C<0.0) C=0.0
	IF(D<0.0) D=0.0
	RETURN
END SUBROUTINE EVKG11


SUBROUTINE EVKG12(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D

	RETURN
END SUBROUTINE EVKG12


SUBROUTINE EVKG13(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D

	RETURN
END SUBROUTINE EVKG13


SUBROUTINE EVKG14(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D

	RETURN
END SUBROUTINE EVKG14


SUBROUTINE EVKG15(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D

	RETURN
END SUBROUTINE EVKG15


SUBROUTINE KSKN1(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,AR,AK,AN,AF,AC,SZX,SZY,SZZ,SLX,SLY,RW
	AR=CS(KM,1)
	AK=CS(KM,2)
	AN=CS(KM,3)
	AF=CS(KM,12)/RAD
	AC=CS(KM,14)
	SZZ=-STS(IE,3)
	IF(SZZ<PAR) SZZ=PAR
	RW=AK*GMAW*(SZZ/PAR)**AN
	SZX=STS(IE,6)
	SZY=STS(IE,5)
	SLX=ABS(SZX/(AC+SZZ*TAN(AF)))
	SLY=ABS(SZY/(AC+SZZ*TAN(AF)))
	IF(SLX>0.90) SLX=0.90
	IF(SLY>0.90) SLY=0.90
	A=RW*(1.0-AR*SLX)**2
	B=RW*(1.0-AR*SLY)**2
	C=1.0E8
	IF(LDG(IE)==-1)THEN
		A=1.0E2
		B=1.0E2
		C=1.0E3
	END IF
	D=EST(IE,4)
	RETURN
END SUBROUTINE KSKN1


SUBROUTINE KSKN2(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,DLT,CA,CB
	DLT=STN(IE,3)
	IF(DLT<=0.0125)THEN
		A=608.0+1400.0
	ELSE
		A=560.0+1400.0
	END IF
	CA=225.0
	CB=40.0
	B=CA/(1.0-CB*DLT)/(1.0-CB*DLT)
	IF(LDG(IE)==-1)THEN
		CA=175.0
		CB=47.6
		C=CA/(1.0-CB*DLT)/(1.0-CB*DLT)
		IF(DLT>0.0115) C=C+600.0
		IF(DLT<=0.0115) C=C+4000.0
	ELSE
		C=1.0E8
	END IF
	D=EST(IE,4)
	RETURN
END SUBROUTINE KSKN2


SUBROUTINE KSKN3(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D,DLT,CA,CB
	DLT=STN(IE,3)
	CA=225.0
	CB=40.0
	A=CA/(1.0-CB*DLT)/(1.0-CB*DLT)
	IF(DLT<=0.0125)THEN
		B=608.0+1400.0
	ELSE
		B=560.0+1400.0
	END IF
	IF(LDG(IE)==-1)THEN
		CA=175.0
		CB=47.6
		C=CA/(1.0-CB*DLT)/(1.0-CB*DLT)
		IF(DLT>0.0115) C=C+600.0
		IF(DLT<=0.0115) C=C+4000.0
	ELSE
		C=1.0E8
	END IF
	D=EST(IE,4)
	RETURN
END SUBROUTINE KSKN3


SUBROUTINE KSKN4(IE,KM,A,B,C,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM
	REAL(KIND=8)::A,B,C,D
	A=CS(KM,2)
	B=CS(KM,3)
	C=CS(KM,4)
	D=EST(IE,4)
	RETURN
END SUBROUTINE KSKN4


SUBROUTINE MATRIX_M2(A,M1,N1,B,M2,N2,AB)
	IMPLICIT NONE
	INTEGER::M1,N1,M2,N2,I,J,K
	REAL(KIND=8),DIMENSION(M1,N1)::A
	REAL(KIND=8),DIMENSION(M2,N2)::B
	REAL(KIND=8),DIMENSION(M1,N2)::AB
	IF(N1/=M2)THEN
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE MATRIX_M2:'
		WRITE(40,*) 'FOR N1/=M2'
		STOP
	END IF
	DO I=1,M1
		DO J=1,N2
			AB(I,J)=0.0
			DO K=1,N1
				AB(I,J)=AB(I,J)+A(I,K)*B(K,J)
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE MATRIX_M2


SUBROUTINE MATRIX_M3(A,M1,N1,B,M2,N2,C,M3,N3,ABC)
	IMPLICIT NONE
	INTEGER::M1,N1,M2,N2,M3,N3
	REAL(KIND=8),DIMENSION(M1,N1)::A
	REAL(KIND=8),DIMENSION(M2,N2)::B
	REAL(KIND=8),DIMENSION(M3,N3)::C
	REAL(KIND=8),DIMENSION(M1,N2)::AB
	REAL(KIND=8),DIMENSION(M1,N3)::ABC
	IF(N1/=M2.OR.N2/=M3)THEN
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE MATRIX_M3:'
		WRITE(40,*) 'FOR N1/=M2 OR N2/=M3'
		STOP
	END IF
	CALL MATRIX_M2(A,M1,N1,B,M2,N2,AB)
	CALL MATRIX_M2(AB,M1,N2,C,M3,N3,ABC)
	RETURN
END SUBROUTINE MATRIX_M3


SUBROUTINE MATRIX_DET(A,DET)
	IMPLICIT NONE
	REAL(KIND=8)::DET
	REAL(KIND=8),DIMENSION(3,3)::A
	DET=0.0
	DET=DET+A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)
	DET=DET-A(3,1)*A(2,2)*A(1,3)-A(1,2)*A(2,1)*A(3,3)-A(1,1)*A(2,3)*A(3,2)
	RETURN
END SUBROUTINE MATRIX_DET


SUBROUTINE MATRIX_I3(DET,A,B)
	IMPLICIT NONE
	INTEGER::I,J
	REAL(KIND=8)::DET
	REAL(KIND=8),DIMENSION(3,3)::A,B
	B(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
	B(2,1)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
	B(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
	B(1,2)=A(1,3)*A(3,2)-A(1,2)*A(3,3)
	B(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
	B(3,2)=A(1,2)*A(3,1)-A(1,1)*A(3,2)
	B(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
	B(2,3)=A(1,3)*A(2,1)-A(1,1)*A(2,3)
	B(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
	DO I=1,3
		DO J=1,3
			B(I,J)=B(I,J)/DET
		END DO
	END DO
	RETURN
END SUBROUTINE MATRIX_I3


SUBROUTINE MATRIX_I6(N,A,C)
	IMPLICIT NONE
	INTEGER::I,J,K,L,N
	REAL(KIND=8)::EPS,R,TM
	INTEGER,DIMENSION(N)::LP
	REAL(KIND=8),DIMENSION(N)::B
	REAL(KIND=8),DIMENSION(N,N)::A,C
	EPS=1.0E-10
	DO I=1,N
		R=A(I,I)
		K=I
		FORALL(J=1:N) B(J)=0.0
		DO J=I,N
			IF(ABS(A(J,I))>ABS(R)) K=J
			R=A(J,I)
		END DO
		LP(I)=K
		B(I)=1.0
		IF(I/=K)THEN
			DO L=1,N
				TM=A(I,L)
				A(I,L)=A(K,L)
				A(K,L)=TM
			END DO
		END IF
		R=A(I,I)
		IF(ABS(R)<=EPS)THEN
			R=SIGN(EPS,R)
			WRITE(40,*) 'ERROR MAY OCCUR IN SUBROUTINE MATRIX_I6:'
			WRITE(40,*) 'FOR A(I,I)<EPS(1.0E-10)'
		END IF
		B(I)=B(I)/R
		DO K=1,N
			A(I,K)=A(I,K)/R
		END DO
		DO J=1,N
			IF(J==I) CYCLE
			TM=A(J,I)
			B(J)=B(J)-TM*B(I)
			DO K=1,N
				A(J,K)=A(J,K)-TM*A(I,K)
			END DO
		END DO
		FORALL(J=1:N) A(J,I)=B(J)
	END DO
	DO I=N,1,-1
		K=LP(I)
		IF(I==K) CYCLE
		DO J=1,N
			TM=A(J,I)
			A(J,I)=A(J,K)
			A(J,K)=TM
		END DO
	END DO
	DO I=1,N
		DO J=1,N
			C(I,J)=A(I,J)
		END DO
	END DO
	RETURN
END SUBROUTINE MATRIX_I6


SUBROUTINE MATRIX_BB(N,A,B)
	IMPLICIT NONE
	INTEGER::J,K,L,N
	REAL(KIND=8),DIMENSION(3,N)::A
	REAL(KIND=8),DIMENSION(6,24)::B
	DO J=1,N
		L=(J-1)*3
		K=L+1
		B(1,K)=-A(1,J)
		B(2,K)=0.0
		B(3,K)=0.0
		B(4,K)=-A(2,J)
		B(5,K)=0.0
		B(6,K)=-A(3,J)
		K=L+2
		B(1,K)=0.0
		B(2,K)=-A(2,J)
		B(3,K)=0.0
		B(4,K)=-A(1,J)
		B(5,K)=-A(3,J)
		B(6,K)=0.0
		K=L+3
		B(1,K)=0.0
		B(2,K)=0.0
		B(3,K)=-A(3,J)
		B(4,K)=0.0
		B(5,K)=-A(2,J)
		B(6,K)=-A(1,J)
	END DO
	RETURN
END SUBROUTINE MATRIX_BB


SUBROUTINE STSLVL(LM,LT,SS,SS0,QP)
	USE COMM
	IMPLICIT NONE
	INTEGER::LM,LT,J
	REAL(KIND=8)::QP
	REAL(KIND=8),DIMENSION(3)::SG
	REAL(KIND=8),DIMENSION(9)::SS
	REAL(KIND=8),DIMENSION(6)::ST,SS0,ST0
	REAL(KIND=8),DIMENSION(15)::CC
	DO J=1,15
		CC(J)=CS(LM,J)
	END DO
	DO J=1,6
		ST(J)=SS(J)
		ST0(J)=SS0(J)
	END DO
	DO J=1,3
		SG(J)=SS(J+6)
	END DO
	IF(LT==0) QP=0.10
	IF(LT==1) CALL STSLVL3(ST,SG,CC,QP)
	IF(LT==2) CALL STSLVL3(ST,SG,CC,QP)
	IF(LT==3) CALL STSLVL3(ST,SG,CC,QP)
	IF(LT==11) CALL STSLVL3(ST,SG,CC,QP)
!	IF(LT==12) CALL STSLVL2(SG,CC,QP)
	IF(LT==13) CALL STSLVL13(ST,CC,QP)
	IF(LT==14) CALL STSLVL14(ST,ST0,CC,QP)
	RETURN
END SUBROUTINE STSLVL


SUBROUTINE STSLVL2(SG,CC,QP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::QP,AF,AC,S3,QP1,QP2
	REAL(KIND=8),DIMENSION(3)::SG
	REAL(KIND=8),DIMENSION(15)::CC
	S3=(SG(1)+SG(2)+SG(3))/3.0
	IF(S3<PAR) S3=PAR
	AF=(CC(12)-CC(13)*LOG10(S3/PAR))/RAD
	AC=CC(14)
	QP1=0.5*(1.0-SIN(AF))*(SG(1)-SG(3))
	QP2=AC*COS(AF)+S3*SIN(AF)
	IF(QP2<0.1*PAR) QP2=0.1*PAR
	QP=QP1/QP2
	IF(QP>0.95) QP=0.95
	RETURN
END SUBROUTINE STSLVL2


SUBROUTINE STSLVL3(SS,SG,CC,QP)
	USE COMM
	IMPLICIT NONE
	REAL(KIND=8)::QP,AF,AC,S3,BF,I1,I2,I3,J2,J3,	&
					&	PP,QQ,QP1,QP2
	REAL(KIND=8),DIMENSION(3)::SG
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8),DIMENSION(15)::CC
	S3=(SG(1)+SG(2)+SG(3))/3.0
	IF(S3<PAR) S3=PAR
	AF=(CC(12)-CC(13)*LOG10(S3/PAR))/RAD
	AC=CC(14)
	CALL INVARS(SS,I1,I2,I3,J2,J3,PP,QQ)
	if(j2.lt.0)then
	print*,j2
	end if
	BF=-1.5*SQRT(3.0)*J3/(J2*SQRT(J2))
	BF=ASIN(BF)/3.0
	QP1=(COS(BF)+SIN(BF)*SIN(AF)/SQRT(3.0))*SQRT(J2)
	QP2=I1*SIN(AF)/3.0+AC*COS(AF)
	IF(QP2<0.1*PAR) QP2=0.1*PAR
	QP=QP1/QP2
	IF(QP>0.95) QP=0.95
	RETURN
END SUBROUTINE STSLVL3


SUBROUTINE STSLVL13(SS,CC,QP)
	USE COMM
	IMPLICIT NONE
	INTEGER::L
	REAL(KIND=8)::QP,AM,AF,AC,I1,I2,I3,J2,J3,PP,QQ,S0
	REAL(KIND=8),DIMENSION(6)::SS,SW
	REAL(KIND=8),DIMENSION(15)::CC
	REAL(KIND=8),DIMENSION(6,6)::TP
	AM=CC(5)
	AF=CC(12)/RAD
	AC=CC(14)
	S0=AC/TAN(AF)
	DO L=1,3
		SS(L)=SS(L)+S0
	END DO
	CALL TRNSFMS(SS,SW,TP)
	CALL INVARS(SW,I1,I2,I3,J2,J3,PP,QQ)
	QP=QQ/PP
	IF(QP<0.01) QP=0.01
	IF(QP>AM) QP=0.99*AM
	RETURN
END SUBROUTINE STSLVL13


SUBROUTINE STSLVL14(SS,SS0,CC,QP)
	USE COMM
	IMPLICIT NONE
	INTEGER::L
	REAL(KIND=8)::QP,AM,AL,AF,AC,I1,I2,I3,J2,J3,PP,PP0,QQ,S0,DLT
	REAL(KIND=8),DIMENSION(6)::SS,SS0,SW,SW0,YT
	REAL(KIND=8),DIMENSION(15)::CC
	REAL(KIND=8),DIMENSION(6,6)::TP
	AM=CC(5)
	AL=CC(7)
	AF=CC(12)/RAD
	AC=CC(14)
	S0=AC/TAN(AF)
	DO L=1,3
		SS0(L)=SS0(L)+S0
		SS(L)=SS(L)+S0
	END DO
	CALL TRNSFMS(SS0,SW0,TP)
	CALL TRNSFMS(SS,SW,TP)
	PP=(SW(1)+SW(2)+SW(3))/3.0
	PP0=(SW0(1)+SW0(2)+SW0(3))/3.0
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=0.0
		SW(L)=SW(L)/PP-DLT
		SW0(L)=SW0(L)/PP0-DLT
	END DO
	DO L=1,6
		YT(L)=SW(L)-AL*SW0(L)
	END DO
	CALL INVARS(YT,I1,I2,I3,J2,J3,PP,QQ)
	PP=0.0
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=2.0
		PP=PP+DLT*SW(L)*YT(L)
	END DO
	QP=1.5*PP/QQ
	IF(QP<0.01) QP=0.01
	IF(QP>AM) QP=0.99*AM
	RETURN
END SUBROUTINE STSLVL14


SUBROUTINE MNSTS(N,A,B)
	IMPLICIT NONE
	INTEGER::N,I,J,K,LT
	REAL(KIND=8)::AII,AIJ,AJJ,AKI,AKJ,BKI,BKJ
	REAL(KIND=8)::AA,RR,CC,SS,GG,EPS
	REAL(KIND=8),DIMENSION(N,N)::A,B
	LT=0
	EPS=1.0E-10
	DO I=1,N
		FORALL(J=1:N) B(I,J)=0.0
		B(I,I)=1.0
	END DO
	DO
		DO I=1,N-1
			DO J=I+1,N
				AII=A(I,I)
				AIJ=A(I,J)
				AJJ=A(J,J)
				AA=AII-AJJ
				IF(ABS(AA)<=EPS) AA=SIGN(EPS,AA)
				RR=2.0*AIJ/AA
				RR=0.5*ATAN(RR)
				CC=COS(RR)
				SS=SIN(RR)
				DO K=1,N
					AKI=A(K,I)*CC+A(K,J)*SS
					AKJ=A(K,J)*CC-A(K,I)*SS
					A(K,I)=AKI
					A(K,J)=AKJ
					A(I,K)=AKI
					A(J,K)=AKJ
					BKI=B(K,I)*CC+B(K,J)*SS
					BKJ=B(K,J)*CC-B(K,I)*SS
					B(K,I)=BKI
					B(K,J)=BKJ
				END DO
				A(I,J)=0.0
				A(J,I)=0.0
				A(I,I)=AII*CC*CC+AJJ*SS*SS+2.0*AIJ*SS*CC
				A(J,J)=AII*SS*SS+AJJ*CC*CC-2.0*AIJ*SS*CC
			END DO
		END DO
		GG=0.0
		DO I=1,N-1
			DO J=I+1,N
				GG=GG+ABS(A(I,J))
			END DO
		END DO
		LT=LT+1
		IF(GG<EPS) EXIT
		IF(LT==50)THEN
			WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE MNSTS(N,A,B):'
			WRITE(40,*) 'FOR LT==50'
			STOP
		END IF
	END DO
	DO J=1,N
		GG=0.0
		DO I=1,N
			GG=GG+B(I,J)*B(I,J)
		END DO
		GG=SQRT(GG)
		DO I=1,N
			B(I,J)=B(I,J)/GG
		END DO
	END DO
	DO I=1,N-1
		K=I
		RR=A(I,I)
		DO J=I+1,N
			IF(A(J,J)<RR) CYCLE
			K=J
			RR=A(J,J)
		END DO
		IF(K/=I)THEN
			RR=A(I,I)
			A(I,I)=A(K,K)
			A(K,K)=RR
			DO J=1,N
				RR=B(J,I)
				B(J,I)=B(J,K)
				B(J,K)=RR
			END DO
		END IF
	END DO
	RETURN
END SUBROUTINE MNSTS


SUBROUTINE MPP(IP)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,ME,IE,LE,IX,J,K
	ME=LQ(IP,2)
	DO IX=1,NX
		JST(IX)=0
	END DO
	DO IE=1,ME
		LE=LUS(IE)
		IF(IST(LE)==2)THEN
			DO J=1,NDS
				K=NOD(LE,J)
				JST(K)=2
			END DO
		END IF
	END DO
	DO IE=1,ME
		LE=LUS(IE)
		IF(IST(LE)==1)THEN
			DO J=1,NDS
				K=NOD(LE,J)
				JST(K)=1
			END DO
		END IF
	END DO
	DO IE=1,ME
		LE=LUS(IE)
		IF(IST(LE)==-1)THEN
			DO J=1,NDS
				K=NOD(LE,J)
				JST(K)=-1
			END DO
		END IF
	END DO
	DO IX=1,NX
		IF(JST(IX)==2)THEN
			NODIS(IX,NFR)=0.0
		ELSE IF(JST(IX)==-1)THEN
			NODIS(IX,NFR)=COR(IX,NDM)-TWH(IX)
		END IF
	END DO
	RETURN
END SUBROUTINE MPP


SUBROUTINE MQQ(IP)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,IP,JJ,KK,LL
	DO IE=1,NE
		LUS(IE)=0
		IST(IE)=0
	END DO
	KK=0
	DO IE=1,NE
		IF(JAD(IE)<IP)THEN
			KK=KK+1
			LUS(KK)=IE
			IST(IE)=1
		ELSE IF(JAD(IE)==IP)THEN
			KK=KK+1
			LUS(KK)=IE
			IST(IE)=2
		END IF
	END DO
	LL=KK
	KK=0
	DO IE=1,LL
		JJ=LUS(IE)
		IF(JRD(JJ)>IP)THEN
			KK=KK+1
			LUS(KK)=LUS(IE)
		END IF
		IF(JRD(JJ)==IP)THEN
			KK=KK+1
			LUS(KK)=LUS(IE)
			IST(JJ)=-1
		END IF
	END DO
	lq(ip,2)=kk
	print*,kk
	IF(LQ(IP,2)/=KK)THEN
		WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE MQQ(IP):'
		WRITE(40,*) 'FOR LQ(IP,2)/=KK','	STEP=',IP
		STOP
	END IF
	RETURN
END SUBROUTINE MQQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!增加子程序MMM！！ADD BY WANG
!!!!!!!!!!!!!!该程序用于计算变带宽存储劲度矩阵 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MMM(IP)
   USE COMM
   IMPLICIT NONE
   INTEGER::IP,I,J,K,L,ME,NN,NV,IW,LL,LE
  !INTEGER,DIMENSION(NDS*NFR)::KE
 !  INTEGER,DIMENSION(NE,NDS)::NEMM
   INTEGER,ALLOCATABLE::NEMM(:,:)
!   INTEGER,DIMENSION(NTEM)::NDEX
!print*,ip

   ME=LQ(IP,2)
   ALLOCATE(NEMM(ne,8))
   K=0
   DO I=1,NX
      JR(I)=0
	  DO J=1,NFR
       JRR(I,J)=0
       K=K+1
	   MA(K)=0
	   MB(K)=0
	  END DO
   END DO
   DO I=1,ME
		LE=LUS(I)
		IF(IST(LE)==-1) CYCLE
		DO J=1,NDS
			K=NOD(LE,J)
			JR(K)=1
		END DO
	END DO
   K=0
   DO I=1,NX
     IF(JR(I)/=1)CYCLE
	 DO J=1,NDM
			IF(ID(I,J)==0) CYCLE
			K=K+1
			JRR(I,J)=K
		END DO
		!	PRINT*,JST(I)
		IF(ID(I,NFR)/=0.AND.JST(I)==1)THEN
			K=K+1

			JRR(I,NFR)=K
		END IF
   END DO
   NN=K
   K=0
   DO I=1,ME
     LE=LUS(I)
     IF(IST(LE)==-1)CYCLE
	  K=K+1
	 DO J=1,NDS
      NEMM(K,J)=NOD(LE,J)
	 END DO
!	 PRINT*,NEMM(K,1)
   END DO
   do i=1,ntem
!    ndex(i)=0
   end do
   CALL MMK(NEMM,k,NN)
   NV=MA(NN)
   NDOF(IP)=NN
   NVOL(IP)=NV
 !  PRINT*,NV,nn
   IF(NV.GT.NTEM)STOP 'THE NTEM IS TOO LITTLE'
   DEALLOCATE(NEMM)
   RETURN

END SUBROUTINE MMM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
SUBROUTINE MMK(KF,ME,NN)
   USE COMM
   IMPLICIT  NONE
   INTEGER::I,J,K,L,NN,NV,LE,LL,ME,IP1,IP2,IPM,ITOV,NTOV
!   INTEGER,DIMENSION(NTEM)::NDEX
  ! INTEGER,DIMENSION(80,NX*NFR)::LAND
 !  INTEGER,DIMENSION(NX*NFR)::KAND
   INTEGER,ALLOCATABLE::LAND(:,:),KAND(:)
   INTEGER,DIMENSION(ne,8)::KF

   ALLOCATE(LAND(80,NX),KAND(NX))
  DO I=1,ME
  ! PRINT*,KF(I,1)
  END DO
 LAND=0
 KAND(:)=0
   DO I=1,ME
!  PRINT*,KF(I,1)
      DO 101 J=1,NDS
          IP1=KF(I,J)
	     DO 102 K=1,NDS
               IP2=KF(I,K)
	            IF(IP2.GT.IP1)CYCLE
                   DO L=1,KAND(IP1)
	                 IF(IP2.EQ.LAND(L,IP1))GOTO 102
	               END DO 
				   KAND(IP1)=KAND(IP1)+1
				   IF(KAND(IP1).GT.80)STOP 'TOO MANY REALTION NODES'
				   LAND(KAND(IP1),IP1)=IP2
	         102 CONTINUE
	      101 CONTINUE 
   END DO
! 调整节点顺序!
  LL=0
  DO I=1,NX
     IF(JR(I).EQ.0)CYCLE
     LE=KAND(I)
	 IF(LE.EQ.0)THEN
         IF(SUM(JRR(I,1:NFR)).NE.0) LL=1
	 ENDIF
     DO J=1,LE
       IPM=LAND(J,I)
	    L=J
	   DO K=J,LE
       IF(LAND(K,I).LT.IPM)THEN
	    L=K
	   IPM=LAND(L,I)
	   END IF 
	   END DO
	   LAND(L,I)=LAND(J,I)
	   LAND(J,I)=IPM
	 END DO
  END DO
IF(LL.EQ.1) STOP 'THE MODEL IS INCORRECT'
LL=0
     DO I=1,NX
	     IF(JR(I).EQ.0)CYCLE
         DO J=1,NFR
          ITOV=JRR(I,J)
          IF(ITOV.EQ.0)CYCLE
          DO K=1,KAND(I)
		  IP2=LAND(K,I)
	          DO L=1,NFR
              NTOV=JRR(IP2,L)
			  IF(NTOV.GT.ITOV)CYCLE
			  IF(NTOV.EQ.0)CYCLE
			  LL=LL+1
	!		  NDEX(LL)=NTOV
             END DO 
          END DO 
		  MA(ITOV)=LL
         END DO
     END DO
!	   PRINT*,'OK'
	 DEALLOCATE(LAND,KAND)
RETURN
END SUBROUTINE MMK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
SUBROUTINE MMKK(IP,NN,NV,NDEX)
   USE COMM
   IMPLICIT  NONE
   INTEGER::I,J,K,L,NN,NV,LE,LL,ME,IP1,IP2,IPM,ITOV,NTOV,IP
   INTEGER,DIMENSION(NV)::NDEX
  ! INTEGER,DIMENSION(80,NX)::LAND
 !  INTEGER,DIMENSION(NX)::KAND
 !  INTEGER,DIMENSION(NE,8)::KF
   INTEGER,ALLOCATABLE::LAND(:,:),KAND(:)
   
ALLOCATE(LAND(80,NX),KAND(NX))
   ME=LQ(IP,2)
LAND=0
KAND(:)=0
   DO I=1,ME
     LE=LUS(I)
	 IF(IST(LE)==-1)CYCLE
      DO 101 J=1,NDS
          IP1=NOD(LE,J)
	     DO 102 K=1,NDS
               IP2=NOD(LE,K)
	            IF(IP2.GT.IP1)CYCLE
                   DO L=1,KAND(IP1)
	                 IF(IP2.EQ.LAND(L,IP1))GOTO 102
	               END DO 
				   KAND(IP1)=KAND(IP1)+1
				   IF(KAND(IP1).GT.80)STOP 'TOO MANY REALTION NODES'
				   LAND(KAND(IP1),IP1)=IP2
	         102 CONTINUE
	      101 CONTINUE 
   END DO
! 调整节点顺序!
  LL=0
  DO I=1,NX
     IF(JR(I).EQ.0)CYCLE
     LE=KAND(I)
	 IF(LE.EQ.0)THEN
         IF(SUM(JRR(I,1:NFR)).NE.0) LL=1
	 ENDIF
     DO J=1,LE
       IPM=LAND(J,I)
	   L=J
	   DO K=J,LE
       IF(LAND(K,I).LT.IPM)THEN
	   L=K
	   IPM=LAND(L,I)
	   END IF 
	   END DO
	   LAND(L,I)=LAND(J,I)
	   LAND(J,I)=IPM
	 END DO
  END DO
IF(LL.EQ.1) STOP 'THE MODEL IS INCORRECT'
LL=0
     DO I=1,NX
	     IF(JR(I).EQ.0)CYCLE
         DO J=1,NFR
          ITOV=JRR(I,J)
          IF(ITOV.EQ.0)CYCLE
      DO K=1,KAND(I)
	  IP2=LAND(K,I)
         DO L=1,NFR
              NTOV=JRR(IP2,L)
			  IF(NTOV.GT.ITOV)CYCLE
			  IF(NTOV.EQ.0)CYCLE
			  LL=LL+1
			  NDEX(LL)=NTOV
             END DO 
          END DO 
!		  MA(ITOV)=LL
         END DO
     END DO
	 DO I=1,NV
!      WRITE(40,*)I,NDEX(I)
	 END DO
	 DEALLOCATE(LAND,KAND)
RETURN
END SUBROUTINE MMKK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

SUBROUTINE MRR(IP)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,I,J,K,L,ME,NN,NV,IW,LL,LE
	INTEGER,DIMENSION(NDS*NFR)::KE
	ME=LQ(IP,2)
	K=0
	DO I=1,NX
		JR(I)=0
		DO J=1,NFR
			JRR(I,J)=0
			K=K+1
			MA(K)=0
			MB(K)=0
		END DO
	END DO
	DO I=1,ME
		LE=LUS(I)
		IF(IST(LE)==-1) CYCLE
		DO J=1,NDS
			K=NOD(LE,J)
			JR(K)=1
		END DO
	END DO
	K=0
	DO I=1,NX
		IF(JR(I)/=1) CYCLE
		DO J=1,NDM
			IF(ID(I,J)==0) CYCLE
			K=K+1
			JRR(I,J)=K
		END DO
		IF(ID(I,NFR)/=0.AND.JST(I)==1)THEN
			K=K+1
			JRR(I,NFR)=K
		END IF
	END DO
	NN=K
	DO I=1,ME
		LE=LUS(I)
		IF(IST(LE)==-1) CYCLE
		IW=NN
		DO J=1,NDS
			L=NOD(LE,J)
			LL=(J-1)*NFR
			DO K=1,NFR
				KE(LL+K)=JRR(L,K)
			END DO
		END DO
		DO J=1,NDS*NFR
			IF(KE(J)==0) CYCLE
			IF(KE(J)<IW) IW=KE(J)
		END DO
		DO J=1,NDS*NFR
			K=KE(J)
			IF(K==0) CYCLE
			LL=K-IW+1
			IF(LL>MA(K)) MA(K)=LL
		END DO
	END DO
	MA(1)=1
	DO K=2,NN
		MA(K)=MA(K)+MA(K-1)
	END DO
	NV=MA(NN)
	NDOF(IP)=NN
	NVOL(IP)=NV
!	print*,nv,nn
	RETURN
END SUBROUTINE MRR


SUBROUTINE DECOMP(NN,NV,SK)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,K,IJ,IK,JK,II,JJ,MI,MJ,IE,IS,JE,JS
	REAL(KIND=8),DIMENSION(NV)::SK
LP1:DO I=2,NN
		II=MA(I)
		IE=I-1
		MI=MA(IE)+I+1-II
		IS=1+MI
		IF(IS<=IE)THEN
LP2:		DO J=IS,IE
				IJ=II-I+J
				MJ=MA(J-1)+J+1-MA(J)
				JS=MAX(MI,MJ)
				JE=J-1
				IF(JS<=JE)THEN
LP3:				DO K=JS,JE
						IK=MA(I)-I+K
						JK=MA(J)-J+K
						SK(IJ)=SK(IJ)-SK(IK)*SK(JK)
					END DO LP3
				END IF
			END DO LP2
		END IF
		IF(MI>IE) CYCLE LP1
LP4:	DO J=MI,IE
			IJ=MA(I)-I+J
			JJ=MA(J)
			SK(IJ)=SK(IJ)/SK(JJ)
			SK(II)=SK(II)-SK(IJ)*SK(IJ)*SK(JJ)
		END DO LP4
	END DO LP1
	RETURN
END SUBROUTINE DECOMP


SUBROUTINE FORBACK(NV,SK,NN,RD)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,IJ,II,MI,IE
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NN)::RD
LP1:DO I=2,NN
		MI=I+MA(I-1)-MA(I)+1
		IE=I-1
		IF(MI>IE)CYCLE LP1
LP2:	DO J=MI,IE
			IJ=MA(I)-I+J
			RD(I)=RD(I)-SK(IJ)*RD(J)
		END DO LP2
	END DO LP1
LP3:DO I=1,NN
		II=MA(I)
		RD(I)=RD(I)/SK(II)
	END DO LP3
LP4:DO I=NN,2,-1
		MI=I+MA(I-1)-MA(I)+1
		IE=I-1
		IF(MI>IE) CYCLE LP4
LP5:	DO J=MI,IE
			IJ=MA(I)-I+J
			RD(J)=RD(J)-SK(IJ)*RD(I)
		END DO LP5
	END DO LP4
	RETURN
END SUBROUTINE FORBACK


SUBROUTINE MSS(IP)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,I,J,K,L,ME,NN,NV,IW,LL,LE
	INTEGER,DIMENSION(NDS*NFR)::KE
	ME=LQ(IP,2)
	K=0
	DO I=1,NX
		JR(I)=0
		DO J=1,NFR
			JRR(I,J)=0
			K=K+1
			MA(K)=0
			MB(K)=0
		END DO
	END DO
	DO I=1,ME
		LE=LUS(I)
		IF(IST(LE)==-1) CYCLE
		DO J=1,NDS
			K=NOD(LE,J)
			JR(K)=1
		END DO
	END DO
	K=0
	DO I=1,NX
		IF(JR(I)/=1) CYCLE
		DO J=1,NDM
			IF(ID(I,J)==0) CYCLE
			K=K+1
			JRR(I,J)=K
		END DO
		IF(ID(I,NFR)/=0.AND.JST(I)==1)THEN
			K=K+1
			JRR(I,NFR)=K
		END IF
	END DO
	NN=K
	DO I=1,ME
		LE=LUS(I)
		IF(IST(LE)==-1) CYCLE
		IW=NN
		DO J=1,NDS
			L=NOD(LE,J)
			LL=(J-1)*NFR
			DO K=1,NFR
				KE(LL+K)=JRR(L,K)
			END DO
		END DO
		DO J=1,NDS*NFR
			IF(KE(J)==0) CYCLE
			IF(KE(J)<IW) IW=KE(J)
		END DO
		DO J=1,NDS*NFR
			K=KE(J)
			IF(K==0) CYCLE
			LL=K-IW+1
			IF(LL>MA(K)) MA(K)=LL
		END DO
	END DO
	MA(1)=1
	DO K=2,NN
		MA(K)=MA(K)+MA(K-1)
	END DO
	NDOF(IP)=NN
	MB(1)=MA(1)
	MB(2)=MA(2)
	NV=MA(NN)-MA(NN-1)-1
	DO I=3,NN
		MB(I)=MA(I-1)-MA(I-2)-1
	END DO
	K=0
	DO I=3,NN
		K=K+MB(I)
		MB(I)=MA(I)+K
	END DO
	NV=NV+MB(NN)
	NDOF(IP)=NN
	NVOL(IP)=NV
	RETURN
END SUBROUTINE MSS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!ADD 迭代法求解器 BY WANG
SUBROUTINE SOL(SK,NV,NN,RD,NEX)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,K,L,II,JI,JS,KS,JE,KE
	REAL(KIND=8)::EPS,YY,RI,RJ
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NN)::RD
	REAL(KIND=8),ALLOCATABLE::ADISP(:)
	INTEGER,DIMENSION(NV)::NEX

	ALLOCATE(ADISP(NN))
    
	DO I=1,NN
       ADISP(I)=0
	END DO
	
	CALL SSOR(SK,NV,NN,RD,NEX,ADISP)

	DO I=1,NN
    RD(I)=ADISP(I)
	END DO
	DEALLOCATE(ADISP)
	RETURN
END SUBROUTINE SOL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SSOR(SK,NV,NN,RD,NEX,ADISP)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,K,L,II,JI,JS,KS,JE,KE,LABLE
	REAL(KIND=8)::EPS,YY,RI,RJ
	REAL(KIND=8),DIMENSION(NV)::SK
!	REAL(KIND=8),DIMENSION(NN)::RD,ADISP,S,G0,H0,D0
	REAL(KIND=8),DIMENSION(NN)::RD,ADISP
	INTEGER,DIMENSION(NV)::NEX
	REAL(KIND=8)::MAX,ZERO,W,TOLER,TOLER2,C,T,ERR,R,BETA
	REAL(KIND=8),ALLOCATABLE::S(:),G0(:),H0(:),D0(:)

	ALLOCATE(S(NN),G0(NN),H0(NN),D0(NN))
	
	MAX=2000
	ZERO=1.0E-7
	W=1.0
	K=0
	LABLE=1
CALL MUL(S,SK,NV,NN,NEX,ADISP,LABLE)
    DO I=1,NN
      G0(I)=S(I)-RD(I)
	  H0(I)=G0(I)
	END DO
CALL TRIA(H0,SK,NV,NN,NEX,LABLE,W)
     DO I=1,NN
       D0(I)=-H0(I)
	 END DO
	 TOLER=0
!	 PRINT*,TOLER
DO I=1,NV
!WRITE(40,*)SK(I)
!IF(SK(I)==0)PRINT*,I
END DO
	 DO I=1,NN
!	 WRITE(40,*)H0(I),G0(I)
     TOLER=TOLER+G0(I)*H0(I)
	 END DO
!    PRINT*,TOLER
100	 K=K+1
!PRINT*,K
CALL MUL(S,SK,NV,NN,NEX,D0,LABLE)
     C=0.0
	 DO I=1,NN
     C=C+D0(I)*S(I)
	 END DO
	 T=TOLER/C
	 DO I=1,NN
     ADISP(I)=ADISP(I)+T*D0(I)
	 END DO
	 IF(MOD(K,1)==0)THEN
     C=0.0
	 ERR=0.0
	 DO I=1,NN
     C=C+ADISP(I)*ADISP(I)
	 ERR=ERR+D0(I)*D0(I)
	 END DO
	 ERR=T*SQRT(ERR)/SQRT(C)
	 IF(ABS(ERR).LE.ZERO.OR.K.GT.MAX) THEN
     CALL MUL(S,SK,NV,NN,NEX,ADISP,LABLE)
     R=0.0
	 DO I=1,NN
      R=R+(S(I)-RD(I))**2
	 END DO
	 R=SQRT(R)
	! PRINT*,R,K
	 RETURN
	 ENDIF
	 END IF
	 
	 DO I=1,NN
     G0(I)=G0(I)+T*S(I)
	 H0(I)=G0(I)
	 END DO
CALL TRIA(H0,SK,NV,NN,NEX,LABLE,W)
     TOLER2=0.0
	 DO I=1,NN
      TOLER2=TOLER2+G0(I)*H0(I)
	 END DO
	  IF(ABS(TOLER2).LE.1.0E-18)RETURN
	  BETA=TOLER2/TOLER
	  DO I=1,NN
	  D0(I)=-H0(I)+BETA*D0(I)
	  END DO
	  TOLER=TOLER2
	  GOTO 100
	  DEALLOCATE(S,G0,H0,D0)
	  RETURN
END SUBROUTINE SSOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MUL(S,SK,NV,NN,NEX,D0,LABLE)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,K,L,II,JI,JS,KS,JE,KE,LABLE,LOW,UPW,NJ
	REAL(KIND=8)::EPS,YY,RI,RJ
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NN)::RD,D0,S
	INTEGER,DIMENSION(NV)::NEX
    REAL(KIND=8)::T
	DO I=1,NN
    S(I)=SK(MA(I))*D0(I)
	END DO
  
  IF(LABLE.EQ.0)THEN
   DO I=1,NN-1
   LOW=MA(I)+1
   UPW=MA(I+1)-1
   T=D0(I)
   DO J=LOW,UPW
   NJ=NEX(J)
   S(I)=S(I)+SK(J)*D0(NJ)
   S(NJ)=S(NJ)+SK(J)*T
   END DO
   END DO
 !  RETURN
   ELSEIF(LABLE.EQ.1)THEN
   DO I=2,NN
   LOW=MA(I-1)+1
   UPW=MA(I)-1
   T=D0(I)
   DO J=LOW,UPW
   NJ=NEX(J)
   S(I)=S(I)+SK(J)*D0(NJ)
   S(NJ)=S(NJ)+SK(J)*T
   END DO

   END DO

  END IF
RETURN
END SUBROUTINE MUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TRIA(HH,SK,NV,NN,NEX,LABLE,W)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,K,L,II,JI,JS,KS,JE,KE,LABLE,LOW,UPW,NJ
	REAL(KIND=8)::EPS,YY,RI,RJ
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NN)::RD,ADISP,HH
	INTEGER,DIMENSION(NV)::NEX
    REAL(KIND=8)::T,W

	DO I=1,NN
    !  WRITE(40,*)HH(I)
	END DO 

    IF(LABLE.EQ.0)THEN
        DO I=1,NN-1
        HH(I)=HH(I)/SK(MA(I))*W
		LOW=MA(I)+1
		UPW=MA(I+1)-1
		T=HH(I)
		DO J=LOW,UPW
        NJ=NEX(J)
		HH(NJ)=HH(NJ)-SK(J)*T
		END DO
		END DO
		HH(NN)=HH(NN)/SK(MA(NN))*W
		DO I=1,NN
        HH(I)=HH(I)*SK(MA(I))/W
		END DO
		HH(NN)=HH(NN)/SK(MA(NN))*W
		DO I=NN-1,1,-1
        LOW=MA(I)+1
		UPW=MA(I+1)-1
		DO J=LOW,UPW
         NJ=NEX(J)
		 HH(I)=HH(I)-SK(J)*HH(NJ)
		END DO
         HH(I)=HH(I)/SK(MA(I))*W
		END DO

    	ELSE IF(LABLE.EQ.1)THEN
	DO I=1,NN
	 low=ma(i)+1
	 upw=ma(i+1)-1
!	 write(40,*)low,upw
	 do j=low,upw
       nj=nex(j)
!	   write(40,*)nj
	 end do
     ! WRITE(40,*)HH(I),MA(I)
	END DO 
	DO I=1,NV
     ! WRITE(40,*)NEX(I)
	END DO 

        HH(1)=HH(1)/SK(1)*W
		DO I=2,NN
          LOW=MA(I-1)+1
		  UPW=MA(I)-1
	      DO J=LOW,UPW
          NJ=NEX(J)
		 ! IF(I==6)WRITE(40,*)HH(I),NJ,SK(J),HH(NJ)
		  HH(I)=HH(I)-SK(J)*HH(NJ)
		  END DO
		 ! IF(I==6)WRITE(40,*)HH(I),SK(MA(I))
		!  WRITE(40,*)HH(I),SK(MA(I))
		  HH(I)=HH(I)/SK(MA(I))*W
		END DO
		DO I=1,NN
		HH(I)=HH(I)*SK(MA(I))/W
        END DO
	DO I=1,NN
  !   WRITE(40,*)HH(I)
	END DO         
        DO I=NN,2,-1
		HH(I)=HH(I)/SK(MA(I))*W
          LOW=MA(I-1)+1
		  UPW=MA(I)-1
         DO J=LOW,UPW
          NJ=NEX(J)
		  HH(NJ)=HH(NJ)-SK(J)*HH(I)
		 END DO
    	END DO
        HH(1)=HH(1)/SK(1)*W   
	END IF
RETURN

END SUBROUTINE TRIA




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SOLVE(NV,SK,NN,RD)
	USE COMM
	IMPLICIT NONE
	INTEGER::NN,NV,I,J,K,L,II,JI,JS,KS,JE,KE
	REAL(KIND=8)::EPS,YY,RI,RJ
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NN)::RD
	EPS=1.0E-10
LP1:DO I=1,NN-1
		II=MB(I)
		RI=SK(II)
		IF(ABS(RI)<=EPS)THEN
			RI=SIGN(EPS,RI)
			WRITE(40,*) 'ERROR MAY OCCUR IN SUBROUTINE SOLVE(NV,SK,NN,RD):'
			WRITE(40,*) 'FOR SK(I,I)<EPS(1.0E-10)'
		END IF
		RI=1.0/RI
LP2:	DO J=I+1,NN
			JS=J-MA(J)+MA(J-1)+1
			IF(I<JS) CYCLE LP2
			JI=MB(J)-(J-I)
			RJ=SK(JI)*RI
LP3:		DO K=I+1,NN
				KS=K-MA(K)+MA(K-1)+1
				IF(I<KS) CYCLE LP3
				L=J
				IF(K>L) L=K
				JE=MB(L)-(J-K)
				KE=MB(K)+(K-I)
				SK(JE)=SK(JE)-SK(KE)*RJ
			END DO LP3
			RD(J)=RD(J)-RD(I)*RJ
		END DO LP2
	END DO LP1
	II=MB(NN)
	RD(NN)=RD(NN)/SK(II)
LP4:DO I=NN-1,1,-1
		YY=0.0
		II=MB(I)
LP5:	DO J=I+1,NN
			JS=J-MA(J)+MA(J-1)+1
			IF(I<JS) CYCLE LP5
			JE=MB(J)+J-I
			YY=YY+RD(J)*SK(JE)
		END DO LP5
		RD(I)=(RD(I)-YY)/SK(II)
	END DO LP4
	RETURN
END SUBROUTINE SOLVE


SUBROUTINE UNBALFRC(LE,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,LE,LM,LT,J,K,L,KJ
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF
	LM=MAT(LE)
	LT=KGSS(LM,1)
	DO J=1,6
		SS(J)=STS(LE,J)
	END DO
	IF(LT>=0)THEN
		IF(ICON(2)/=0)THEN
			DO K=1,3
				SS(K)=SS(K)+PWS(LE)
			END DO
		END IF
		CALL NDFRCI(LE,SS,FF)
	ELSE
		IF(LT==-1.OR.LT==-2.OR.LT==-3) THEN
			CALL NDFRCJ(LE,SS,FF)
		END IF
		IF(LT==-4)THEN
			CALL NDFRCK(LE,SS,FF)
		END IF
	END IF
	RETURN
END SUBROUTINE UNBALFRC


SUBROUTINE NDFRCI(IE,SS,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,I,J,K,L,NI,NJ,NK,s
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF,DF
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV,DERIX
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB,INV
	REAL(KIND=8),DIMENSION(6,NDS*NDM)::BB
	REAL(KIND=8)::GHX,GHY,GHZ,GX,GY,GZ,DET,HH
	s=jplane(ie)
	DO I=1,NDS*NDM
		FF(I)=0.0
	END DO
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	DO NI=1,NGS
		GHX=WGH(NI)
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GHY=WGH(NJ)
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GHZ=WGH(NK)
				GZ=WXYZ(NK)
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				CALL MATRIX_I3(DET,JCB,INV)
				CALL MATRIX_M2(INV,NDM,NDM,DERIV,NDM,NDS,DERIX)
				CALL MATRIX_BB(NDS,DERIX,BB)
				DO I=1,NDS*NDM
					DF(I)=0.0
					DO J=1,6
						DF(I)=DF(I)+BB(J,I)*SS(J)
					END DO
				END DO
				HH=GHX*GHY*GHZ*DET
				DO I=1,NDS*NDM
					FF(I)=FF(I)+DF(I)*HH
				END DO
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE NDFRCI


SUBROUTINE NDFRCJ(IE,SS,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,I,J,K,L
	REAL(KIND=8)::DLT,S,BL
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8),DIMENSION(NDM)::ST,DF
	REAL(KIND=8),DIMENSION(NDM,NDM)::TR,RT
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF
	DO I=1,NDS*NDM
		FF(I)=0.0
	END DO
	ST(1)=SS(6)
	ST(2)=SS(5)
	ST(3)=SS(3)
	IF(ST(3)>0.0) RETURN
	DO J=1,NDM
		DO K=1,NDM
			TR(J,K)=DMT(IE,NDM+J,NDM+K)
			RT(K,J)=TR(J,K)
		END DO
	END DO
	S=EST(IE,4)
	BL=1.0/REAL(NFS)
	DO L=1,NDM
		ST(L)=ST(L)*S*BL
	END DO
	DO J=1,NDM
		DF(J)=0.0
		DO K=1,NDM
			DF(J)=DF(J)+RT(J,K)*ST(K)
		END DO
	END DO
	DO J=1,NDS
		K=(J-1)*NDM
		DLT=1.0
		IF(J<=NFS) DLT=-1.0
		DO L=1,NDM
			FF(K+L)=DLT*DF(L)
		END DO
	END DO
	RETURN
END SUBROUTINE NDFRCJ


SUBROUTINE NDFRCK(IE,SS,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,I,J,K,L
	REAL(KIND=8)::SZ,SE,SA,SL,PP,BL,DLT
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8),DIMENSION(NDM)::XLMN,DF
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF
	DO I=1,NDS*NDM
		FF(I)=0.0
	END DO
	SZ=SS(3)
	SE=EST(IE,1)
	SA=EST(IE,3)
	SL=EST(IE,4)
	XLMN(1)=DMT(IE,4,4)
	XLMN(2)=DMT(IE,5,5)
	XLMN(3)=DMT(IE,6,6)
	PP=SZ*SA
	BL=1.0/REAL(NFS)
	DO L=1,NDM
		DF(L)=PP*XLMN(L)*BL
	END DO
	DO J=1,NDS
		K=(J-1)*NDM
		DLT=1.0
		IF(J<=NFS) DLT=-1.0
		DO L=1,NDM
			FF(K+L)=DLT*DF(L)
		END DO
	END DO
	RETURN
END SUBROUTINE NDFRCK


SUBROUTINE NDFRCL(IE,SS,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF

	RETURN
END SUBROUTINE NDFRCL


SUBROUTINE BNDAPP(IE,STF,RR)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,N,JJ,JL
	INTEGER,DIMENSION(NDS*NFR)::KE
	REAL(KIND=8),DIMENSION(NDS*NFR)::DS,RR
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	N=0
	DO J=1,NDS
		K=NOD(IE,J)
		JJ=NFR*(J-1)
		DO L=1,NFR
			JL=JJ+L
			IF(JRR(K,L)==0) N=N+1
			KE(JL)=JRR(K,L)
			DS(JL)=NODIS(K,L)
		END DO
	END DO
	IF(N==0) RETURN
	DO J=1,NDS*NFR
		IF(KE(J)==0) CYCLE
		DO K=1,NDS*NFR
			IF(KE(K)/=0) CYCLE
			RR(J)=RR(J)+STF(J,K)*DS(K)
		END DO
	END DO
	RETURN
END SUBROUTINE BNDAPP


SUBROUTINE INVARS(SS,H1,H2,H3,O2,O3,PP,QQ)
	IMPLICIT NONE
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8)::H1,H2,H3,O2,O3,PP,QQ
	REAL(KIND=8)::SXX,SYY,SZZ,TXY,TYZ,TZX
	SXX=SS(1)
	SYY=SS(2)
	SZZ=SS(3)
	TXY=SS(4)
	TYZ=SS(5)
	TZX=SS(6)
	H1=SXX+SYY+SZZ
	H2=-SXX*SYY-SYY*SZZ-SZZ*SXX+TXY*TXY+TYZ*TYZ+TZX*TZX
	H3=SXX*SYY*SZZ+2.0*TXY*TYZ*TZX-SXX*TYZ*TYZ-SYY*TZX*TZX-SZZ*TXY*TXY
	O2=H2+H1*H1/3.0
	O3=H3+H1*H2/3.0+2.0*H1*H1*H1/27.0
	PP=H1/3.0
	QQ=SQRT(3.0*O2)
	RETURN
END SUBROUTINE INVARS


SUBROUTINE INVARN(SS,VV,RR)
	IMPLICIT NONE
	REAL(KIND=8),DIMENSION(6)::SS
	REAL(KIND=8)::VV,RR,SXX,SYY,SZZ,TXY,TYZ,TZX,R1,R2,R3,R4
	SXX=SS(1)
	SYY=SS(2)
	SZZ=SS(3)
	TXY=SS(4)
	TYZ=SS(5)
	TZX=SS(6)
	VV=SXX+SYY+SZZ
	R1=(SXX-SYY)**2
	R2=(SYY-SZZ)**2
	R3=(SZZ-SXX)**2
	R4=1.5*(TXY**2+TYZ**2+TZX**2)
	RR=SQRT(2.0*(R1+R2+R3+R4)/9.0)
	RETURN
END SUBROUTINE INVARN


SUBROUTINE DSGIJ(SS,DI1,DI2,DI3,DJ2,DJ3)
	IMPLICIT NONE
	REAL(KIND=8),DIMENSION(6)::SS,DI1,DI2,DI3,DJ2,DJ3
	REAL(KIND=8)::I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::SXX,SYY,SZZ,TXY,TYZ,TZX
	CALL INVARS(SS,I1,I2,I3,J2,J3,PP,QQ)
	SXX=SS(1)
	SYY=SS(2)
	SZZ=SS(3)
	TXY=SS(4)
	TYZ=SS(5)
	TZX=SS(6)
	DI1(1)=1.0
	DI1(2)=1.0
	DI1(3)=1.0
	DI1(4)=0.0
	DI1(5)=0.0
	DI1(6)=0.0
	DI2(1)=-(SYY+SZZ)
	DI2(2)=-(SXX+SZZ)
	DI2(3)=-(SXX+SYY)
	DI2(4)=2.0*TXY
	DI2(5)=2.0*TYZ
	DI2(6)=2.0*TZX
	DI3(1)=SYY*SZZ-TYZ*TYZ
	DI3(2)=SXX*SZZ-TZX*TZX
	DI3(3)=SXX*SYY-TXY*TXY
	DI3(4)=2.0*TYZ*TZX-2.0*SZZ*TXY
	DI3(5)=2.0*TXY*TZX-2.0*SXX*TYZ
	DI3(6)=2.0*TXY*TYZ-2.0*SYY*TZX
	DJ2(1)=SXX-PP
	DJ2(2)=SYY-PP
	DJ2(3)=SZZ-PP
	DJ2(4)=2.0*TXY
	DJ2(5)=2.0*TYZ
	DJ2(6)=2.0*TZX
	DJ3(1)=DI3(1)+I1*DI2(1)/3.0+I2*DI1(1)/3.0+2.0*I1*I1*DI1(1)/9.0
	DJ3(2)=DI3(2)+I1*DI2(2)/3.0+I2*DI1(2)/3.0+2.0*I1*I1*DI1(2)/9.0
	DJ3(3)=DI3(3)+I1*DI2(3)/3.0+I2*DI1(3)/3.0+2.0*I1*I1*DI1(3)/9.0
	DJ3(4)=DI3(4)+I1*DI2(4)/3.0+I2*DI1(4)/3.0+2.0*I1*I1*DI1(4)/9.0
	DJ3(5)=DI3(5)+I1*DI2(5)/3.0+I2*DI1(5)/3.0+2.0*I1*I1*DI1(5)/9.0
	DJ3(6)=DI3(6)+I1*DI2(6)/3.0+I2*DI1(6)/3.0+2.0*I1*I1*DI1(6)/9.0
	RETURN
END SUBROUTINE DSGIJ


SUBROUTINE DSGPQ(SS,DP,DQ)
	IMPLICIT NONE
	REAL(KIND=8),DIMENSION(6)::SS,DP,DQ
	REAL(KIND=8)::I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::SXX,SYY,SZZ,TXY,TYZ,TZX
	CALL INVARS(SS,I1,I2,I3,J2,J3,PP,QQ)
	SXX=SS(1)
	SYY=SS(2)
	SZZ=SS(3)
	TXY=SS(4)
	TYZ=SS(5)
	TZX=SS(6)
	DP(1)=0.3333333333333333
	DP(2)=0.3333333333333333
	DP(3)=0.3333333333333333
	DP(4)=0.0
	DP(5)=0.0
	DP(6)=0.0
	DQ(1)=1.5*(SXX-PP)/QQ
	DQ(2)=1.5*(SYY-PP)/QQ
	DQ(3)=1.5*(SZZ-PP)/QQ
	DQ(4)=3.0*TXY/QQ
	DQ(5)=3.0*TYZ/QQ
	DQ(6)=3.0*TZX/QQ
	RETURN
END SUBROUTINE DSGPQ


SUBROUTINE DSGPY(ALF,SS,SS0,DP,DY)
	IMPLICIT NONE
	INTEGER::L
	REAL(KIND=8),DIMENSION(6)::SS,SS0,YT,DP,DY
	REAL(KIND=8)::ALF,PP,QQ,PP0,DLT,R
	REAL(KIND=8)::I1,I2,I3,J2,J3
	PP=(SS(1)+SS(2)+SS(3))/3.0
	PP0=(SS0(1)+SS0(2)+SS0(3))/3.0
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=0.0
		SS(L)=SS(L)/PP-DLT
		SS0(L)=SS0(L)/PP0-DLT
	END DO
	DO L=1,6
		YT(L)=SS(L)-ALF*SS0(L)
	END DO
	CALL INVARS(YT,I1,I2,I3,J2,J3,DLT,QQ)
	DP(1)=0.3333333333333333
	DP(2)=0.3333333333333333
	DP(3)=0.3333333333333333
	DP(4)=0.0
	DP(5)=0.0
	DP(6)=0.0
	R=0.0
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=2.0
		R=R+DLT*YT(L)*SS(L)
	END DO
	DY(1)=(3.0*YT(1)-R)/(2.0*PP*QQ)
	DY(2)=(3.0*YT(2)-R)/(2.0*PP*QQ)
	DY(3)=(3.0*YT(3)-R)/(2.0*PP*QQ)
	DY(4)=3.0*YT(4)/(PP*QQ)
	DY(5)=3.0*YT(5)/(PP*QQ)
	DY(6)=3.0*YT(6)/(PP*QQ)
	RETURN
END SUBROUTINE DSGPY


SUBROUTINE MORCLM(RFI,A,B)
	IMPLICIT NONE
	REAL(KIND=8)::RFI,ALF,I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::FI1,FJ2,FJ3,EPS,CS3,A1,A2,A3,B1,B2,B3
	REAL(KIND=8),DIMENSION(6)::A,B,DI1,DI2,DI3,DJ2,DJ3
	EPS=1.0E-3
	CALL INVARS(A,I1,I2,I3,J2,J3,PP,QQ)
	CALL DSGIJ(A,DI1,DI2,DI3,DJ2,DJ3)
	ALF=-1.5*SQRT(3.0)*J3/(J2*SQRT(J2))
	ALF=ASIN(ALF)/3.0
	CS3=COS(3.0*ALF)
	IF(CS3<EPS) CS3=EPS
	A1=SIN(RFI)
	A2=COS(ALF)+SIN(ALF)*SIN(RFI)/SQRT(3.0)
	A3=-SIN(ALF)+COS(ALF)*SIN(RFI)/SQRT(3.0)
	B1=0.5/SQRT(J2)
	B2=SQRT(27.0)*J3/(4.0*J2*SQRT(J2)*CS3)
	B3=SQRT(3.0)/(2.0*SQRT(J2)*CS3)
	FI1=-A1/3.0
	FJ2=A2*B1+A3*B2
	FJ3=-A3*B3
	B(1)=FI1*DI1(1)+FJ2*DJ2(1)+FJ3*DJ3(1)
	B(2)=FI1*DI1(2)+FJ2*DJ2(2)+FJ3*DJ3(2)
	B(3)=FI1*DI1(3)+FJ2*DJ2(3)+FJ3*DJ3(3)
	B(4)=FI1*DI1(4)+FJ2*DJ2(4)+FJ3*DJ3(4)
	B(5)=FI1*DI1(5)+FJ2*DJ2(5)+FJ3*DJ3(5)
	B(6)=FI1*DI1(6)+FJ2*DJ2(6)+FJ3*DJ3(6)
	RETURN
END SUBROUTINE MORCLM


SUBROUTINE DRKPRG(ALF,A,B)
	IMPLICIT NONE
	REAL(KIND=8)::ALF,I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::SXX,SYY,SZZ,TXY,TYZ,TZX
	REAL(KIND=8),DIMENSION(6)::A,B
	CALL INVARS(A,I1,I2,I3,J2,J3,PP,QQ)
	J2=SQRT(J2)
	SXX=A(1)-PP
	SYY=A(2)-PP
	SZZ=A(3)-PP
	TXY=A(4)
	TYZ=A(5)
	TZX=A(6)
	B(1)=0.5*SXX/J2-ALF
	B(2)=0.5*SYY/J2-ALF
	B(3)=0.5*SZZ/J2-ALF
	B(4)=TXY/J2
	B(5)=TYZ/J2
	B(6)=TZX/J2
	RETURN
END SUBROUTINE DRKPRG


SUBROUTINE STHMOD(A,B,C)
	IMPLICIT NONE
	INTEGER::L
	REAL(KIND=8)::I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::FR,FS,DFP1,DFQ1,DFP2,DFQ2
	REAL(KIND=8),DIMENSION(6)::A,B,C,DP,DQ
	CALL INVARS(A,I1,I2,I3,J2,J3,PP,QQ)
	FR=2.0
	FS=2.0
	DFP1=2.0*PP
	DFQ1=2.0*FR*FR*QQ
	DFP2=-QQ**FS/(PP*PP)
	DFQ2=FS*QQ**(FS-1.0)/PP
	CALL DSGPQ(A,DP,DQ)
	DO L=1,6
		B(L)=DFP1*DP(L)+DFQ1*DQ(L)
		C(L)=DFP2*DP(L)+DFQ2*DQ(L)
	END DO
	RETURN
END SUBROUTINE STHMOD


SUBROUTINE CAMBRG(FL,FK,FE,FM,A,B)
	IMPLICIT NONE
	INTEGER::L
	REAL(KIND=8)::I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::FL,FK,FE,FM,RR,DFP,DFQ
	REAL(KIND=8),DIMENSION(6)::A,B,DP,DQ
	CALL INVARS(A,I1,I2,I3,J2,J3,PP,QQ)
	RR=(FL-FK)/(1.0+FE)
	DFP=RR*(FM*FM*PP*PP-QQ*QQ)/(FM*FM*PP*PP+QQ*QQ)/PP
	DFQ=2.0*RR*QQ/(FM*FM*PP*PP+QQ*QQ)
	CALL DSGPQ(A,DP,DQ)
	DO L=1,6
		B(L)=DFP*DP(L)+DFQ*DQ(L)
	END DO
	RETURN
END SUBROUTINE CAMBRG


SUBROUTINE SEKIOHTA(ALF,FL,FK,FE,FM,A,B,C)
	IMPLICIT NONE
	INTEGER::L
	REAL(KIND=8)::PP,RR
	REAL(KIND=8)::ALF,FL,FK,FE,FM,DFP,DFY
	REAL(KIND=8),DIMENSION(6)::A,B,C,DP,DY
	PP=A(1)+A(2)+A(3)
	RR=(FL-FK)/(1.0+FE)
	DFP=RR/PP
	DFY=RR/FM
	CALL DSGPY(ALF,A,B,DP,DY)
	DO L=1,6
		C(L)=DFP*DP(L)+DFY*DY(L)
	END DO
	RETURN
END SUBROUTINE SEKIOHTA


SUBROUTINE HOHAI
	IMPLICIT NONE
	INTEGER::L

	RETURN
END SUBROUTINE HOHAI


!!!!!!!!!!!!!!!!采用变带宽法存储劲度矩阵!!!!!!!!!!!!!!!!!!!!!!!!!!!ADD BY WANG 
 SUBROUTINE SUMKK(IE,EK,NV,SK,NEX)
	USE COMM
	IMPLICIT NONE
	INTEGER::J,K,L,LL,IE,NV,JW,KW,JJ,JK
	INTEGER,DIMENSION(32)::KE
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::EK
	INTEGER,DIMENSION(NV)::NEX
	DO J=1,NDS
		L=NOD(IE,J)
		LL=NFR*(J-1)
		DO K=1,NFR
			KE(LL+K)=JRR(L,K)
		END DO
	END DO

    DO J=1,NDS*NFR
	JW=KE(J)
	IF(JW.EQ.0)CYCLE
         DO L=1,NDS*NFR
           JK=KE(L)
		   IF(JK.EQ.0)CYCLE
           IF(KE(J)<KE(L))CYCLE
		   IF(KE(J).EQ.1)THEN
              SK(1)=SK(1)+EK(J,L)
			  CYCLE
		   END IF
	!	   PRINT*,KE(J)
		   DO K=MA(KE(J)-1)+1,MA(KE(J))
           IF(KE(L).EQ.NEX(K))THEN
              SK(K)=SK(K)+EK(J,L)
		   END IF
		   END DO
         END DO
    END DO
	DO IE=1,NV
	
	END DO 
    RETURN
END SUBROUTINE SUMKK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SUMRR(IE,EK,NV,SK)
	USE COMM
	IMPLICIT NONE
	INTEGER::J,K,L,LL,IE,NV,JW,KW,JJ,JK
	INTEGER,DIMENSION(32)::KE
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::EK
	DO J=1,NDS
		L=NOD(IE,J)
		LL=NFR*(J-1)
		DO K=1,NFR
			KE(LL+K)=JRR(L,K)
		END DO
	END DO
LPJ:DO J=1,NDS*NFR
		JW=KE(J)
		IF(JW==0) CYCLE LPJ
		JJ=MA(JW)
LPK:	DO K=1,NDS*NFR
			KW=KE(K)
			IF(KW==0) CYCLE LPK
			L=JW-KW
			IF(L<0) CYCLE LPK
			JK=JJ-L
			SK(JK)=SK(JK)+EK(J,K)
		END DO LPK
	END DO LPJ
	DO IE=1,NV
!    IF(ABS(SK(IE))>1E-7)WRITE(40,*)SK(IE)
	END DO
	RETURN
END SUBROUTINE SUMRR


SUBROUTINE SUMSS(IE,EK,NV,SK)
	USE COMM
	IMPLICIT NONE
	INTEGER::J,K,L,LL,IE,NV,JW,JJ,KW,KK,JK
	INTEGER,DIMENSION(NDS*NFR)::KE
	REAL(KIND=8),DIMENSION(NV)::SK
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::EK
	DO J=1,NDS
		L=NOD(IE,J)
		LL=NFR*(J-1)
		DO K=1,NFR
			KE(LL+K)=JRR(L,K)
		END DO
	END DO
LPJ:DO J=1,NDS*NFR
		JW=KE(J)
		IF(JW==0) CYCLE LPJ
		JJ=MB(JW)
LPK:	DO K=1,NDS*NFR
			KW=KE(K)
			IF(KW==0) CYCLE LPK
			KK=MB(KW)
			L=JW-KW
			IF(L>=0)THEN
				JK=JJ-L
			ELSE
				JK=KK-L
			END IF
			SK(JK)=SK(JK)+EK(J,K)
		END DO LPK
	END DO LPJ
	RETURN
END SUBROUTINE SUMSS


SUBROUTINE SUMLD(IE,GQ,NN,RD)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,NN,I,J,K,L,LL,ME
	INTEGER,DIMENSION(NDS*NFR)::KE
	REAL(KIND=8),DIMENSION(NN)::RD
	REAL(KIND=8),DIMENSION(NDS*NFR)::GQ
	DO J=1,NDS
		L=NOD(IE,J)
		LL=NFR*(J-1)
		DO K=1,NFR
			KE(LL+K)=JRR(L,K)
		END DO
	END DO
	DO J=1,NDS*NFR
		LL=KE(J)
		IF(LL==0) CYCLE
		RD(LL)=RD(LL)+GQ(J)
	END DO
	RETURN
END SUBROUTINE SUMLD


SUBROUTINE STFNS(IE,STF,DETA)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,NI,NJ,LM,LP,JJ,KK,LL,NN,JK1,JK2,LN1,LN2
	REAL(KIND=8)::DETA
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS*NDM)::KS
	REAL(KIND=8),DIMENSION(NDS,NDS)::KH
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS)::KC
	REAL(KIND=8),DIMENSION(NDS*NFR,NDS*NFR)::STF
	LM=MAT(IE)
	LP=KGSS(LM,2)
	DO J=1,NDS*NFR
		DO K=1,NDS*NFR
			STF(J,K)=0.0
		END DO
	END DO
	DO J=1,NDS*NDM
		DO K=1,NDS*NDM
			KS(J,K)=0.0	
		END DO
	END DO
	DO J=1,NDS*NDM
		DO K=1,NDS
			KC(J,K)=0.0	
		END DO
	END DO
	DO J=1,NDS
		DO K=1,NDS
			KH(J,K)=0.0
		END DO
	END DO
	CALL SDKS(IE,KS)
	DO J=1,NDS
		JJ=NDM*(J-1)
		LL=NFR*(J-1)
		DO K=1,NDS
			KK=NDM*(K-1)
			NN=NFR*(K-1)
			DO NI=1,NDM
				JK1=JJ+NI
				LN1=LL+NI
				DO NJ=1,NDM
					JK2=KK+NJ
					LN2=NN+NJ
					STF(LN1,LN2)=KS(JK1,JK2)
				END DO
			END DO
		END DO
	END DO
	IF(LP>0)THEN
		if(kgss(lm,1).ge.0)	CALL SDKC(IE,KC)
		CALL SDKH(IE,KH)
		DO J=1,NDS
			JJ=J*NFR
			DO K=1,NDS
				KK=K*NFR
				STF(JJ,KK)=-GATA*DETA*KH(J,K)
			END DO
		END DO
		DO J=1,NDS
			JJ=NDM*(J-1)
			LL=NFR*(J-1)
			DO K=1,NDM
				JK1=JJ+K
				LN1=LL+K			
				DO L=1,NDS
					JK2=L
					LN2=L*NFR
					STF(LN1,LN2)=KC(JK1,JK2)
					STF(LN2,LN1)=KC(JK1,JK2)
				END DO
			END DO
		END DO
	END IF
	RETURN
END SUBROUTINE STFNS


SUBROUTINE DMAT(IE,KM,KT,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT,J,K
	REAL(KIND=8)::A,B
	REAL(KIND=8),DIMENSION(6,6)::D
	IF(KT<=10)THEN
		A=EST(IE,1)
		B=EST(IE,2)
		CALL DMATE(A,B,D)
	ELSE
		IF(KT==11) CALL DMATP11(IE,KM,D)
		IF(KT==12) CALL DMATP12(IE,KM,D)
		IF(KT==13) CALL DMATP13(IE,KM,D)
		IF(KT==14) CALL DMATP14(IE,KM,D)
	END IF
	DO J=1,6
		DO K=1,6
			DMT(IE,J,K)=D(J,K)
		END DO
	END DO
	RETURN
END SUBROUTINE DMAT


SUBROUTINE CMATE(A,B,C)
	USE COMM
	IMPLICIT NONE
	INTEGER::I,J
	REAL(KIND=8)::A,B
	REAL(KIND=8),DIMENSION(6,6)::C
	DO I=1,6
		DO J=1,6
			C(I,J)=0.0
		END DO
	END DO
	IF(B>0.49) B=0.49
	IF(B<0.01) B=0.01
	C(1,1)=1.0
	C(2,2)=1.0
	C(3,3)=1.0
	C(1,2)=-1.0*B
	C(1,3)=-1.0*B
	C(2,1)=-1.0*B
	C(2,3)=-1.0*B
	C(3,1)=-1.0*B
	C(3,2)=-1.0*B
	C(4,4)=2.0*(1.0+B)
	C(5,5)=2.0*(1.0+B)
	C(6,6)=2.0*(1.0+B)
	DO I=1,6
		DO J=1,6
			C(I,J)=C(I,J)/A
		END DO
	END DO
	RETURN
END SUBROUTINE CMATE


SUBROUTINE DMATE(A,B,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::J,K
	REAL(KIND=8)::A,B,R
	REAL(KIND=8),DIMENSION(6,6)::D
	DO J=1,6
		DO K=1,6
			D(J,K)=0.0
		END DO
	END DO
	IF(B>0.49) B=0.49
	IF(B<0.01) B=0.01
	R=(A/(1.0+B))/(1.0-2.0*B)
	D(1,1)=R*(1.0-B)
	D(2,2)=D(1,1)
	D(3,3)=D(1,1)
	D(1,2)=R*B
	D(1,3)=R*B
	D(2,1)=R*B
	D(2,3)=R*B
	D(3,1)=R*B
	D(3,2)=R*B
	D(4,4)=0.5*R*(1.0-2.0*B)
	D(5,5)=0.5*R*(1.0-2.0*B)
	D(6,6)=0.5*R*(1.0-2.0*B)
	RETURN
END SUBROUTINE DMATE


SUBROUTINE DMATP11(IE,KM,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,J,K,L,LK
	REAL(KIND=8)::EE,VV,A1,A2,X1,X2,ALF,RW
	REAL(KIND=8),DIMENSION(6)::ST,DF1,DF2
	REAL(KIND=8),DIMENSION(6,6)::CE,CP,C,D
	DO J=1,6
		DO K=1,6
			CE(J,K)=0.0
			CP(J,K)=0.0
			C(J,K)=0.0
			D(J,K)=0.0
		END DO
	END DO
	EE=EST(IE,1)
	VV=EST(IE,2)
	A1=EST(IE,3)
	A2=EST(IE,4)
	X1=PLS(IE,1)
	X2=PLS(IE,2)
	IF(LDG(IE)==-1)THEN
		CALL DMATE(EE,VV,D)
		RETURN
	END IF
	CALL CMATE(EE,VV,CE)
	DO L=1,6
		DF1(L)=0.0
		DF2(L)=0.0
		ST(L)=STS(IE,L)
	END DO
	CALL STHMOD(ST,DF1,DF2)
	X1=1.0-X1
	X2=1.0-X2
	DO J=1,6
		DO K=1,6
			CP(J,K)=X1*A1*DF1(J)*DF1(K)+X2*A2*DF2(J)*DF2(K)
		END DO
	END DO
	ALF=1.0
	DO
		LK=0
		DO L=1,6
			RW=CE(L,L)+ALF*CP(L,L)
			IF(RW>0.0) CYCLE
			LK=1
			EXIT
		END DO
		IF(LK==0) EXIT
		ALF=0.9*ALF
	END DO
	DO J=1,6
		DO K=1,6
			C(J,K)=CE(J,K)+ALF*CP(J,K)
		END DO
	END DO
	CALL MATRIX_I6(6,C,D)
	RETURN
END SUBROUTINE DMATP11


SUBROUTINE DMATP12(IE,KM,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,J,K,L,LK
	REAL(KIND=8)::EE,VV,A1,A2,X1,X2,ALF,RW
	REAL(KIND=8),DIMENSION(6)::ST,DF1,DF2
	REAL(KIND=8),DIMENSION(6,6)::CE,CP,C,D
	DO J=1,6
		DO K=1,6
			CE(J,K)=0.0
			CP(J,K)=0.0
			C(J,K)=0.0
			D(J,K)=0.0
		END DO
	END DO
	RETURN
END SUBROUTINE DMATP12


SUBROUTINE DMATP13(IE,KM,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,J,K,L,LK
	REAL(KIND=8),DIMENSION(6)::ST,SS,DF
	REAL(KIND=8),DIMENSION(6,6)::CE,CP,CT,TP,C,D
	REAL(KIND=8)::BL,BK,BE,BM,BMF,BC,BF,EE,VV,S0,	&
						R1,R2,R3,QP,DL,X1,ALF,RW
	DO J=1,6
		DO K=1,6
			CE(J,K)=0.0
			CP(J,K)=0.0
			CT(J,K)=0.0
			C(J,K)=0.0
			D(J,K)=0.0
		END DO
	END DO
	BL=CS(KM,1)
	BK=CS(KM,2)
	BE=CS(KM,3)
	BM=CS(KM,4)
	BMF=CS(KM,5)
	VV=CS(KM,6)
	BF=CS(KM,12)/RAD
	BC=CS(KM,14)
	X1=PLS(IE,1)
	S0=STS(IE,10)
	IF(S0<0.1*PAR) S0=0.1*PAR
	EE=3.0*(1.0-2.0*VV)*(1.0+BE)*S0/BK
	IF(LDG(IE)==-1)THEN
		CALL DMATE(EE,VV,D)
		RETURN
	END IF
	CALL CMATE(EE,VV,CE)
	DO L=1,6
		DF(L)=0.0
		ST(L)=STS(IE,L)
	END DO
	S0=BC/TAN(BF)
	DO L=1,3
		ST(L)=ST(L)+S0
	END DO
	CALL TRNSFMS(ST,SS,TP)
	CALL CAMBRG(BL,BK,BE,BM,SS,DF)
	QP=STS(IE,12)
	IF(QP>=BMF) QP=0.99*BMF
	R1=(BMF/BM)**4
	R2=(BM**4-QP**4)/(BMF**4-QP**4)
	R3=DF(1)+DF(2)+DF(3)
	DL=R1*R2/R3
	DO J=1,6
		DO K=1,6
			CT(J,K)=DL*DF(J)*DF(K)
		END DO
	END DO
	CALL MATRIX_M2(CT,6,6,TP,6,6,CP)
	X1=1.0-X1
	DO J=1,6
		DO K=1,6
			CP(J,K)=X1*CP(J,K)
		END DO
	END DO
	ALF=1.0
	DO
		LK=0
		DO L=1,6
			RW=CE(L,L)+ALF*CP(L,L)
			IF(RW>0.0) CYCLE
			LK=1
			EXIT
		END DO
		IF(LK==0) EXIT
		ALF=0.9*ALF
	END DO
	DO J=1,6
		DO K=1,6
			C(J,K)=CE(J,K)+ALF*CP(J,K)
		END DO
	END DO
	CALL MATRIX_I6(6,C,D)
	RETURN
END SUBROUTINE DMATP13


SUBROUTINE DMATP14(IE,KM,D)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,J,K,L,LK
	REAL(KIND=8),DIMENSION(6,6)::CE,CP,CT,TP,C,D
	REAL(KIND=8),DIMENSION(6)::ST,ST0,SS,SS0,DF
	REAL(KIND=8)::BL,BK,BE,BM,BMF,BC,BF,AL,EE,VV,	&
					S0,R1,R2,R3,QP,DL,X1,ALF,RW
	DO J=1,6
		DO K=1,6
			CE(J,K)=0.0
			CP(J,K)=0.0
			CT(J,K)=0.0
			C(J,K)=0.0
			D(J,K)=0.0
		END DO
	END DO
	BL=CS(KM,1)
	BK=CS(KM,2)
	BE=CS(KM,3)
	BM=CS(KM,4)
	BMF=CS(KM,5)
	VV=CS(KM,6)
	AL=CS(KM,7)
	BF=CS(KM,12)/RAD
	BC=CS(KM,14)
	X1=PLS(IE,1)
	S0=STS(IE,10)
	IF(S0<0.1*PAR) S0=0.1*PAR
	EE=3.0*(1.0-2.0*VV)*(1.0+BE)*S0/BK
	IF(LDG(IE)==-1)THEN
		CALL DMATE(EE,VV,D)
		RETURN
	END IF
	CALL CMATE(EE,VV,CE)
	DO L=1,6
		DF(L)=0.0
		ST(L)=STS(IE,L)
		ST0(L)=STR(IE,L)
	END DO
	S0=BC/TAN(BF)
	DO L=1,3
		ST(L)=ST(L)+S0
		ST0(L)=ST0(L)+S0
	END DO
	CALL TRNSFMS(ST0,SS0,TP)
	CALL TRNSFMS(ST,SS,TP)
	CALL SEKIOHTA(AL,BL,BK,BE,BM,SS,SS0,DF)
	QP=STS(IE,12)
	IF(QP>=BMF) QP=0.99*BMF
	R1=(BMF/BM)**4
	R2=(BM**4-QP**4)/(BMF**4-QP**4)
	R3=DF(1)+DF(2)+DF(3)
	DL=R1*R2/R3
	DO J=1,6
		DO K=1,6
			CT(J,K)=DL*DF(J)*DF(K)
		END DO
	END DO
	CALL MATRIX_M2(CT,6,6,TP,6,6,CP)
	X1=1.0-X1
	DO J=1,6
		DO K=1,6
			CP(J,K)=X1*CP(J,K)
		END DO
	END DO
	ALF=1.0
	DO
		LK=0
		DO L=1,6
			RW=CE(L,L)+ALF*CP(L,L)
			IF(RW>0.0) CYCLE
			LK=1
			EXIT
		END DO
		IF(LK==0) EXIT
		ALF=0.9*ALF
	END DO
	DO J=1,6
		DO K=1,6
			C(J,K)=CE(J,K)+ALF*CP(J,K)
		END DO
	END DO
	CALL MATRIX_I6(6,C,D)
	RETURN
END SUBROUTINE DMATP14


SUBROUTINE KMAT(IE,AK)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KT,LT,I,J
	REAL(KIND=8),DIMENSION(NDM,NDM)::AK
	KT=MAT(IE)
	LT=KGSS(KT,2)
!	if(ie==1514)print*,lt
	DO I=1,NDM
		DO J=1,NDM
			AK(I,J)=0.0
		END DO
	END DO
	SELECT CASE(LT)
		CASE(0)
			RETURN
		CASE(1)
			CALL KMAT1(IE,KT,AK)
		CASE(2)
			CALL KMAT2(IE,KT,AK)
		CASE DEFAULT
			WRITE(40,*) 'ERROR MAY OCCUR IN SUBROUTINE KMAT(IE,AK)'
			WRITE(40,*) 'IE=',IE
			STOP
	END SELECT
	RETURN
END SUBROUTINE KMAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1增加渗透系数变化子程序
SUBROUTINE KCAGE(IE,KT,SKF)
   USE COMM
   IMPLICIT NONE 
   INTEGER::IE,KT
   REAL(KIND=8)::EV,EN,EKS,ED
   REAL(KIND=8),DIMENSION(NDM)::SKF
   EV=STN(IE,10)
   EN=CK(KT,15)
   EKS=CK(KT,1)
   EN=EN/(1+EN)
   ED=(EN+EV)*(1-EN)/EN/(1-EN-EV)
   ed=ed**3
   if(ed>=100)ed=100
   if(ed<=1)ed=1
   SKF(1)=EKS/ED
   SKF(2)=SKF(1)
   SKF(3)=SKF(1)
   RETURN
END SUBROUTINE KCAGE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE KMAT2(IE,KT,AK)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KT
	REAL(KIND=8),DIMENSION(NDM,NDM)::AK
    REAL(KIND=8),DIMENSION(NDM)::SKF
	CALL KCAGE(IE,KT,SKF)
	AK(1,1)=skf(1)
	AK(1,2)=0
	AK(1,3)=0
	AK(2,1)=0
	AK(2,2)=skf(2)
	AK(2,3)=0
	AK(3,1)=0
	AK(3,2)=0
	AK(3,3)=skf(3)
	RETURN
END SUBROUTINE KMAT2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111111
SUBROUTINE KMAT1(IE,KT,AK)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KT
	REAL(KIND=8),DIMENSION(NDM,NDM)::AK
	AK(1,1)=seep(ie,1)
	AK(1,2)=0
	AK(1,3)=0
	AK(2,1)=0
	AK(2,2)=seep(ie,2)
	AK(2,3)=0
	AK(3,1)=0
	AK(3,2)=0
	AK(3,3)=seep(ie,3)
	RETURN
END SUBROUTINE KMAT1


SUBROUTINE SDKS(IE,KS)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS*NDM)::KS
	KM=MAT(IE)
	KT=KGSS(KM,1)
	IF(KT>=0)THEN
		CALL SDKSI(IE,KM,KT,KS)
	ELSE
		CALL SDKSJ(IE,KM,KT,KS)
	END IF
	RETURN
END SUBROUTINE SDKS


SUBROUTINE SDKSI(IE,KM,KT,KS)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT,I,J,K,L,NI,NJ,NK,s
	REAL(KIND=8)::GHX,GHY,GHZ,GX,GY,GZ,DET,RR
	REAL(KIND=8),DIMENSION(6,6)::D
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS*NDM)::KS,KSD
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV,DERIX
	REAL(KIND=8),DIMENSION(6,NDS*NDM)::BB
	REAL(KIND=8),DIMENSION(NDS*NDM,6)::BT
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB,INV
	s=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	CALL DMAT(IE,KM,KT,D)
	DO NI=1,NGS
		GHX=WGH(NI)
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GHY=WGH(NJ)
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GHZ=WGH(NK)
				GZ=WXYZ(NK)
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				CALL MATRIX_I3(DET,JCB,INV)
				CALL MATRIX_M2(INV,NDM,NDM,DERIV,NDM,NDS,DERIX)
				CALL MATRIX_BB(NDS,DERIX,BB)
				FORALL(J=1:6,K=1:NDS*NDM) BT(K,J)=BB(J,K)
				CALL MATRIX_M3(BT,NDS*NDM,6,D,6,6,BB,6,NDS*NDM,KSD)
				RR=GHX*GHY*GHZ*DET
				DO I=1,NDS*NDM
					DO J=1,NDS*NDM
						KS(I,J)=KS(I,J)+KSD(I,J)*RR
					END DO
				END DO
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SDKSI


SUBROUTINE SDKSJ(IE,KM,KT,KS)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT,J,K,L
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS*NDM)::KS
	SELECT CASE(KT)
		CASE(-1)
			CALL SDKSJ3(IE,KM,KS)
		CASE(-2)
			CALL SDKSJ3(IE,KM,KS)
		CASE(-3)
			CALL SDKSJ3(IE,KM,KS)
		CASE(-4)
			CALL SDKSJ4(IE,KM,KS)
		CASE DEFAULT
			WRITE(40,*) 'ERROR MAY OCCUR IN SUBROUTINE SDKSJ(IE,KM,KT,KS)'
			WRITE(40,*) 'IE=',IE
			STOP
	END SELECT
	RETURN
END SUBROUTINE SDKSJ


SUBROUTINE SDKC(IE,KC)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,I,J,K,L,NI,NJ,NK,s
	REAL(KIND=8)::GHX,GHY,GHZ,GX,GY,GZ,DET,RR
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS)::KC,KCD
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDS*NDM)::BU
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV,DERIX
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB,INV
	s=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	DO NI=1,NGS
		GHX=WGH(NI)
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GHY=WGH(NJ)
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GHZ=WGH(NK)
				GZ=WXYZ(NK)
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				CALL MATRIX_I3(DET,JCB,INV)
				CALL MATRIX_M2(INV,NDM,NDM,DERIV,NDM,NDS,DERIX)
				DO I=1,NDS
					K=NDM*(I-1)
					DO J=1,NDM
						L=K+J
						BU(L)=-DERIX(J,I)
					END DO
				END DO
				RR=GHX*GHY*GHZ*DET
				DO I=1,NDS*NDM
					DO J=1,NDS
						KCD(I,J)=BU(I)*SH(J)*GMAW
					END DO
				END DO
				DO I=1,NDS*NDM
					DO J=1,NDS
						KC(I,J)=KC(I,J)+RR*KCD(I,J)
					END DO
				END DO
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SDKC


SUBROUTINE SDKH(IE,KH)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KT,I,J,K,L,NI,NJ,NK,mt,jj,jk,k1,l1,l2,lj,lk,jp
	REAL(KIND=8)::GHX,GHY,GHZ,GX,GY,GZ,DET,RR,b,bl,k2,k3,k4
	REAL(KIND=8),DIMENSION(NDM,NDM)::KK,RRK,TR,rt
	REAL(KIND=8),DIMENSION(NDS,NDS)::KH,KHD
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV,DERIX
	REAL(KIND=8),DIMENSION(NDS,NDM)::DERIY
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB,INV
	mt=mat(ie)
	kt=kgss(mt,1)
    jp=jplane(ie)
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			XYZ(J,K)=COR(L,K)
		END DO
	END DO
	CALL KMAT(IE,KK)
	if(kt>=0)then
	do j=1,ndm
    do k=1,ndm
    rrk(j,k)=kk(j,k)
	end do
	end do
    else
	DO J=1,NDM
		DO K=1,NDM
			TR(J,K)=DMT(IE,NDM+J,NDM+K)
			RT(K,J)=TR(J,K)
		END DO
	END DO
	CALL MATRIX_M3(RT,NDM,NDM,kk,NDM,NDM,TR,NDM,NDM,RRk)
	endif
	DO NI=1,NGS
		GHX=WGH(NI)
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GHY=WGH(NJ)
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GHZ=WGH(NK)
				GZ=WXYZ(NK)
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,jp)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				CALL MATRIX_I3(DET,JCB,INV)
				CALL MATRIX_M2(INV,NDM,NDM,DERIV,NDM,NDS,DERIX)
				FORALL(J=1:NDM,K=1:NDS) DERIY(K,J)=DERIX(J,K)
				CALL MATRIX_M3(DERIY,NDS,NDM,rrk,NDM,NDM,DERIX,NDM,NDS,KHD)
				RR=GHX*GHY*GHZ*DET
				DO I=1,NDS
					DO J=1,NDS
						KH(I,J)=KH(I,J)+KHD(I,J)*RR*GMAW
					END DO
				END DO
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SDKH


SUBROUTINE ADRSLT(IP)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,ME,IX,IE,LE,J,KT
    REAL(KIND=8),DIMENSION(NDM)::SKF
	ME=LQ(IP,2)
	DO IX=1,NX
		DO J=1,NDM
			DIS(IX,J)=DIS(IX,J)+DDIS(IX,J)
		END DO
		TWH(IX)=TWH(IX)+DTWH(IX)
		IF(TWH(IX)<COR(IX,NDM)) TWH(IX)=COR(IX,NDM)
	END DO
	DO IE=1,NE
	  KT=MAT(IE)
	  if(kgss(kt,2).ne.2)cycle
      CALL KCAGE(IE,KT,SKF)
      SEEP(IE,1)=SKF(1)
	  SEEP(IE,2)=SKF(2)
	  SEEP(IE,3)=SKF(3)
	END DO
	DO IE=1,ME
		LE=LUS(IE)
		IF(IST(LE)==-1) CYCLE
		CALL STNCALC(LE)
		CALL STSCALC(LE)
	END DO
	IF(ICON(2)/=0)THEN
		DO IE=1,ME
			LE=LUS(IE)
			CALL WPCALC(LE)
		END DO
	END IF
	RETURN
END SUBROUTINE ADRSLT


SUBROUTINE DDISCALC(IX,NN,RD)
	USE COMM
	IMPLICIT NONE
	INTEGER::IX,J,NN,JJ
	REAL(KIND=8)::DS,TH
	REAL(KIND=8),DIMENSION(NN)::RD
	DO J=1,NFR
		JJ=JRR(IX,J)
		DS=0.0
		IF(JJ>0)THEN
			DS=RD(JJ)
		ELSE
			DS=NODIS(IX,J)
		END IF
		IF(J<NFR)THEN
			DDIS(IX,J)=DS
		ELSE
			DTWH(IX)=DS
		END IF
	END DO
	RETURN
END SUBROUTINE DDISCALC


SUBROUTINE DSTNCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT
	LM=MAT(IE)
	LT=KGSS(LM,1)
	IF(LT>=0)THEN
		CALL DSTNCALCI(IE)
	ELSE
		IF(LT==-1) CALL DSTNCALCJ(IE)
		IF(LT==-2) CALL DSTNCALCJ(IE)
		IF(LT==-3) CALL DSTNCALCJ(IE)
		IF(LT==-4) CALL DSTNCALCK(IE)
	END IF
	RETURN
END SUBROUTINE DSTNCALC


SUBROUTINE DSTNCALCI(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,LL,NI,NJ,NK,s,k1,k2,k3,i,ic,ic1,ic2,ic3,ic4,ic5,ic6,ii,jj,jl
	REAL(KIND=8)::GHX,GHY,GHZ,GX,GY,GZ,DET,RR
	REAL(KIND=8),DIMENSION(6)::DST
	REAL(KIND=8),DIMENSION(NDM*NDS)::UU
	REAL(KIND=8),DIMENSION(NDS)::SH
	REAL(KIND=8),DIMENSION(NDM,NDS)::DERIV,DERIX
	REAL(KIND=8),DIMENSION(NDS,NDM)::XYZ
	REAL(KIND=8),DIMENSION(NDM,NDM)::JCB,INV
	REAL(KIND=8),DIMENSION(6,NDS*NDM)::B
	REAL(KIND=8),DIMENSION(48,NDS*NDM)::Bz
	REAL(KIND=8),DIMENSION(6,NDS)::strain
	s=jplane(ie)
	DO J=1,6
	!	DST(J)=0.0
	END DO
	LL=0
	DO J=1,NDS
		L=NOD(IE,J)
		DO K=1,NDM
			LL=LL+1
			XYZ(J,K)=COR(L,K)
			UU(LL)=DDIS(L,K)
		END DO
	END DO
!!!!!!!!!!!!!!!!
ic=-6
!!!!!!!!!!!!!!!!
	DO NI=1,NGS
		GHX=WGH(NI)
		GX=WXYZ(NI)
		DO NJ=1,NGS
			GHY=WGH(NJ)
			GY=WXYZ(NJ)
			DO NK=1,NGS
				GHZ=WGH(NK)
				GZ=WXYZ(NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!
ic=ic+6
!!!!!!!!!!!!!!!!!!!!!!!!
				CALL SHAPEFUNC(GX,GY,GZ,SH,DERIV,s)
				CALL MATRIX_M2(DERIV,NDM,NDS,XYZ,NDS,NDM,JCB)
				CALL MATRIX_DET(JCB,DET)
				CALL MATRIX_I3(DET,JCB,INV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
	K3=0
DO K=1,8
	K3=K3+3
	K2=K3-1
	K1=K2-1
	DO L=1,3
	  B(L,K1)=0
	  B(L,K2)=0
	  B(L,K3)=0 
    end do
	  DO I=1,3
	     B(1,K1)=B(1,K1)-INV(1,I)*DERIV(I,K)
	     B(2,K2)=B(2,K2)-INV(2,I)*DERIV(I,K)
	     B(3,K3)=B(3,K3)-INV(3,I)*DERIV(I,K)
      end do
	    B(4,K1)=B(2,K2)
	    B(4,K2)=B(1,K1)
	    B(4,K3)=0
	    B(5,K1)=0
	    B(5,K2)=B(3,K3)
	    B(5,K3)=B(2,K2)
	    B(6,K1)=B(3,K3)
	    B(6,K2)=0
	    B(6,K3)=B(1,K1)
end	do		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IC1=IC+1
	IC2=IC+2
	IC3=IC+3
	IC4=IC+4
	IC5=IC+5
	IC6=IC+6
	DO I=1,24
	BZ(IC1,I)=B(1,I)
	BZ(IC2,I)=B(2,I)
	BZ(IC3,I)=B(3,I)
	BZ(IC4,I)=B(4,I)
	BZ(IC5,I)=B(5,I)
	BZ(IC6,I)=B(6,I)			
	end do					
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
	!			CALL MATRIX_M2(INV,NDM,NDM,DERIV,NDM,NDS,DERIX)
	!			CALL MATRIX_BB(NDS,DERIX,BB)
!				DO J=1,6
!					DO K=1,NDM*NDS
!						DST(J)=DST(J)+BB(J,K)*UU(K)
!					END DO
!				END DO
			END DO
		END DO
	END DO
!	RR=REAL(NGS*NGS*NGS)

DO  JJ=1,8
	DO II=1,6
	  STRAIN(II,JJ)=0.0 
    end do
end do

DO  J=1,8
	JL=(J-1)*6
  DO K=1,24
	DO L=1,6
	  LL=JL+L
	  STRAIN(L,J)=STRAIN(L,J)+BZ(LL,K)*uu(k) !计算得到总的应变增量STRAIN(6,8)
    end do
  end do
end do

	DO J=1,6
	    dst(j)=0
	  do l=1,8
        dst(j)=dst(j)+strain(j,l)/real(nds)
	  end do
	    DSTN(IE,J)=DST(J)               !/RR
	END DO
	RETURN
END SUBROUTINE DSTNCALCI


SUBROUTINE DSTNCALCJ(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,JJ,KK
	REAL(KIND=8),DIMENSION(NDM,NDM)::TR
	REAL(KIND=8),DIMENSION(NDM)::A,B
	REAL(KIND=8)::DU,DV,DW
	DU=0.0
	DV=0.0
	DW=0.0  
	DO J=1,NDM
		DO K=1,ndm 
			TR(J,K)=DMT(IE,NDM+J,NDM+K)
		END DO
	END DO
	DO J=1,NFS
		JJ=NOD(IE,J)
		K=J+NFS
		KK=NOD(IE,K)
		DU=DU+DDIS(KK,1)-DDIS(JJ,1)
		DV=DV+DDIS(KK,2)-DDIS(JJ,2)
		DW=DW+DDIS(KK,3)-DDIS(JJ,3) 
	END DO 
	DU=DU/REAL(NFS)
	DV=DV/REAL(NFS)
	DW=DW/REAL(NFS)
	A(1)=DU
	A(2)=DV
	A(3)=DW
	DO J=1,NDM
		B(J)=0.0
		DO K=1,NDM
			B(J)=B(J)+TR(J,K)*A(K)
		END DO
	END DO
	DU=B(1)
	DV=B(2)
	DW=B(3)
	DSTN(IE,3)=DW
	DSTN(IE,5)=DV
	DSTN(IE,6)=DU
	RETURN
END SUBROUTINE DSTNCALCJ


SUBROUTINE DSTNCALCK(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,I,KL
	REAL(KIND=8)::DLT,AE,AA,AL,BL
	REAL(KIND=8),DIMENSION(NDM)::XLMN
	REAL(KIND=8),DIMENSION(NDS*NDM)::BB,UU
	XLMN(1)=DMT(IE,4,4)
	XLMN(2)=DMT(IE,5,5)
	XLMN(3)=DMT(IE,6,6)
	AE=EST(IE,1)
	AA=EST(IE,3)
	AL=EST(IE,4)
	BL=1.0/REAL(NFS)
	DO J=1,NDS
		I=NOD(IE,J)
		K=(J-1)*NDM
		DLT=1.0
		IF(J<=NFS) DLT=-1.0
		DO L=1,NDM
			KL=K+L
			BB(KL)=DLT*XLMN(L)
			UU(KL)=DDIS(I,L)
		END DO
	END DO
	DLT=0.0
	DO J=1,NDS*NDM
		DLT=DLT+BB(J)*UU(J)
	END DO
	DSTN(IE,3)=DLT*BL/AL
	RETURN
END SUBROUTINE DSTNCALCK


SUBROUTINE DSTSCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT
	LM=MAT(IE)
	LT=KGSS(LM,1)
	IF(LT>=0)THEN
		CALL DSTSCALCI(IE)
	ELSE
		IF(LT==-1) CALL DSTSCALCJ(IE)
		IF(LT==-2) CALL DSTSCALCJ(IE)
		IF(LT==-3) CALL DSTSCALCJ(IE)
		IF(LT==-4) CALL DSTSCALCK(IE)
	END IF
	RETURN
END SUBROUTINE DSTSCALC


SUBROUTINE DSTSCALCI(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,J,K,L,lm,mt
	REAL(KIND=8),DIMENSION(6,6)::DP
	REAL(KIND=8),DIMENSION(6)::ST,SS
	mt=mat(ie)
	lm=kgss(mt,6)
	DO J=1,6
		DO K=1,6
			DP(J,K)=DMT(IE,J,K)
		END DO
	END DO
	DO K=1,6
	if(lm.eq.0)then
		ST(K)=DSTN(IE,K)
    else
	    st(k)=dstn(ie,k)-rhstn(ie,k)      !减去流变产生的应变
	end if
	END DO
	DO J=1,6
		SS(J)=0.0
		DO K=1,6
			SS(J)=SS(J)+DP(J,K)*ST(K)
		END DO
	END DO
	DO L=1,6
		DSTS(IE,L)=SS(L)
	END DO
	RETURN
END SUBROUTINE DSTSCALCI


SUBROUTINE DSTSCALCJ(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,K,L
	REAL(KIND=8),DIMENSION(NDM,NDM)::D
	DO K=1,NDM
		DO L=1,NDM
			D(K,L)=DMT(IE,K,L)
		END DO
	END DO
	DSTS(IE,3)=D(3,3)*DSTN(IE,3)
	DSTS(IE,6)=D(1,1)*DSTN(IE,6)
	DSTS(IE,5)=D(2,2)*DSTN(IE,5)	
	RETURN
END SUBROUTINE DSTSCALCJ


SUBROUTINE DSTSCALCK(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE
	REAL(KIND=8)::AE,AA,AL
	AE=EST(IE,1)
	AA=EST(IE,3)
	AL=DSTN(IE,3)
	DSTS(IE,3)=AE*AL
	RETURN
END SUBROUTINE DSTSCALCK


SUBROUTINE STNCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT
	KM=MAT(IE)
	KT=KGSS(KM,1)
	IF(KT>=0)THEN
		CALL STNCALCI(IE)
	ELSE
		IF(KT==-1) CALL STNCALCJ(IE)
		IF(KT==-2) CALL STNCALCJ(IE)
		IF(KT==-3) CALL STNCALCJ(IE)
		IF(KT==-4) CALL STNCALCK(IE)
	END IF
	RETURN
END SUBROUTINE STNCALC


SUBROUTINE STNCALCI(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,K,L
	REAL(KIND=8)::VV,RR
	REAL(KIND=8),DIMENSION(6)::RSN
	REAL(KIND=8),DIMENSION(NDM,NDM)::ST,VEST
	DO L=1,6
		STN(IE,L)=STN(IE,L)+DSTN(IE,L)
		RSN(L)=STN(IE,L)
	END DO
	ST(1,1)=STN(IE,1)
	ST(2,2)=STN(IE,2)
	ST(3,3)=STN(IE,3)
	ST(1,2)=0.5*STN(IE,4)
	ST(1,3)=0.5*STN(IE,6)
	ST(2,1)=0.5*STN(IE,4)
	ST(2,3)=0.5*STN(IE,5)
	ST(3,1)=0.5*STN(IE,6)
	ST(3,2)=0.5*STN(IE,5)
	CALL MNSTS(NDM,ST,VEST)
	DO K=1,NDM
		DO L=1,NDM
			VCN(IE,K,L)=VEST(K,L)
		END DO
	END DO
	STN(IE,7)=ST(1,1)
	STN(IE,8)=ST(2,2)
	STN(IE,9)=ST(3,3)
	CALL INVARN(RSN,VV,RR)
	STN(IE,10)=VV
	STN(IE,11)=RR
	IF(ABS(VV)<=1.0E-6) VV=SIGN(1.0E-6,VV)
	STN(IE,12)=RR/VV
	RETURN
END SUBROUTINE STNCALCI


SUBROUTINE STNCALCJ(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,L
	DO L=1,6
		STN(IE,L)=STN(IE,L)+DSTN(IE,L)
	END DO
	IF(STN(IE,3)<0.0) STN(IE,3)=0.0
	RETURN
END SUBROUTINE STNCALCJ


SUBROUTINE STNCALCK(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE
	STN(IE,3)=STN(IE,3)+DSTN(IE,3)
	RETURN
END SUBROUTINE STNCALCK


SUBROUTINE STSCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT
	KM=MAT(IE)
	KT=KGSS(KM,1)
	IF(KT>=0)THEN
		CALL STSCALCI(IE)
	ELSE
		IF(KT==-1) CALL STSCALCJ(IE)
		IF(KT==-2) CALL STSCALCJ(IE)
		IF(KT==-3) CALL STSCALCJ(IE)
		IF(KT==-4) CALL STSCALCK(IE)
	END IF
	RETURN
END SUBROUTINE STSCALC


SUBROUTINE STSCALCI(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,K,L
	REAL(KIND=8)::I1,I2,I3,J2,J3,PP,QQ,QP
	REAL(KIND=8),DIMENSION(9)::RST
	REAL(KIND=8),DIMENSION(6)::RS,RST0,rt
	REAL(KIND=8),DIMENSION(NDM,NDM)::SS,VESS
	LM=MAT(IE)
	LT=KGSS(LM,1)
	DO L=1,6
		STS(IE,L)=STS(IE,L)+DSTS(IE,L)
	END DO
	SS(1,1)=STS(IE,1)
	SS(2,2)=STS(IE,2)
	SS(3,3)=STS(IE,3)
	SS(1,2)=STS(IE,4)
	SS(1,3)=STS(IE,6)
	SS(2,1)=STS(IE,4)
	SS(2,3)=STS(IE,5)
	SS(3,1)=STS(IE,6)
	SS(3,2)=STS(IE,5)
	CALL MNSTS(NDM,SS,VESS)
	DO K=1,NDM
		DO L=1,NDM
			VCS(IE,K,L)=VESS(K,L)
		END DO
	END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,6
    rt(k)=sts(ie,k)
	end do
    !call stsmod(ie,ss,rt)
    do k=1,6
    sts(ie,k)=rt(k)
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	STS(IE,7)=SS(1,1)
	STS(IE,8)=SS(2,2)
	STS(IE,9)=SS(3,3)
	DO K=1,9
		RST(K)=STS(IE,K)
	END DO
	DO L=1,6
		RS(L)=STS(IE,L)
		RST0(L)=STR(IE,L)
	END DO
	CALL INVARS(RS,I1,I2,I3,J2,J3,PP,QQ)
	STS(IE,10)=PP
	STS(IE,11)=QQ
    thata(ie)=asin(-3*sqrt(3.0)/2.0*j3/j2**(1.5))/3.0
	CALL STSLVL(LM,LT,RST,RST0,QP)
	STS(IE,12)=QP
	RETURN
END SUBROUTINE STSCALCI

subroutine stsmod(ie,ss,rt)
use comm
implicit none
real(kind=8)::ra,rb,rc,pa,pb,dhh,bbb,pc,rm,a
integer::ie,lm,lt,i
real(kind=8),dimension(ndm,ndm)::ss
real(kind=8),dimension(6)::rt(6)
!print*,ie,ss(3,3)
	lm=mat(ie)
	lt=kgss(lm,1)
	pa=cs(lm,12)
	pb=cs(lm,14)
    pc=cs(lm,13)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!破坏单元应力修正
    rm=ss(3,3)*0.01
!	print*,ie,ss(3,3)
	if(rm.lt.1.0)rm=1.0
!	print*,rm
	!if(ie==632)print*,pa,pc,rm,lm
	pa=(pa-pc*log10(rm))/57.3
    if(pa.le.0.0) then
	write(40,33) ie
    stop
	endif
33	format(3x,'error occured in element',i8)
    a=-pb/tan(pa)
	if(ss(3,3).gt.a.or.lt.eq.0)then
	return
	end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ra=ss(1,1)
	if(ra.lt.20)ra=str(ie,1)
	rc=0.7854-0.5*pa
	rc=tan(rc)
    rc=rc*(rc*ra-2.0*pb)
    dhh=ss(1,1)-ss(3,3)
	if(dhh.lt.0.1) dhh=0.1
	bbb=(ss(1,1)-ss(2,2))/dhh
	if(bbb.lt.0.0)bbb=0.0
	if(bbb.gt.0.9)bbb=0.9
    rb=ra-(ra-rc)*bbb
    call sfai(ss,ra,rb,rc,ie,rt)
    ss(1,1)=ra
	ss(2,2)=rb
	ss(3,3)=rc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
return
end subroutine


   SUBROUTINE SFAI(HH,RA,RB,RC,L,rt)
	 use comm
	 implicit none
	 integer::i,j,k,l
	 real(kind=8)::tcc,ta,tb,tc,ra,rb,rc
	 real(kind=8),dimension(ndm,ndm)::hh,dk
	 real(kind=8),dimension(6)::rt
	 tcc=abs(sts(l,4))+abs(sts(l,5))+abs(sts(l,6))
	 do 46 i=1,3
	 do 46 j=1,3
46	 dk(i,j)=0.0
     if(tcc.lt.1.0e-4)then
     dk(1,3)=1.0
	 dk(2,2)=1.0
	 dk(3,1)=1.0
	 goto 761
	 end if
	 do 76 j=1,3
	 tc=(hh(j,j)-sts(l,2))*sts(l,5)+sts(l,4)*sts(l,6)
	 if(abs(tc).lt.10.0) tc=sign(10.0,tc)
	 ta=(hh(j,j)-sts(l,1))*sts(l,4)+sts(l,5)*sts(l,6)/tc
	 tb=(hh(j,j)-sts(l,1))*(hh(j,j)-sts(l,2))-sts(l,6)**2/tc
	 tc=1.0+ta**2+tb**2
     dk(j,1)=1.0/sqrt(tc)
	 dk(j,2)=ta*dk(j,1)
76	 dk(j,3)=tb*dk(j,1)
761  do 77 j=1,3
77   rt(j)=DK(1,J)*DK(1,J)*RA+DK(2,J)*DK(2,J)*RB+DK(3,J)*DK(3,J)*RC
     rt(4)=DK(1,2)*DK(1,3)*RA+DK(2,2)*DK(2,3)*RB+DK(3,2)*DK(3,3)*RC
     rt(5)=DK(1,3)*DK(1,1)*RA+DK(2,3)*DK(2,1)*RB+DK(3,3)*DK(3,1)*RC
     rt(6)=DK(1,1)*DK(1,2)*RA+DK(2,1)*DK(2,2)*RB+DK(3,1)*DK(3,2)*RC
     RETURN
   END SUBROUTINE SFAI


SUBROUTINE STSCALCJ(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,L
	DO L=1,6
		STS(IE,L)=STS(IE,L)+DSTS(IE,L)
	END DO
	IF(STS(IE,3)>-PAR) STS(IE,3)=-PAR
	RETURN
END SUBROUTINE STSCALCJ


SUBROUTINE STSCALCK(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM
	LM=MAT(IE)
	STS(IE,3)=STS(IE,3)+DSTS(IE,3)
	IF(STS(IE,3)>CS(LM,14)) STS(IE,3)=CS(LM,14)
	IF(STS(IE,3)<-CS(LM,14)) STS(IE,3)=-CS(LM,14)
	RETURN
END SUBROUTINE STSCALCK


SUBROUTINE WPCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KP,J,L
	REAL(KIND=8)::PP
	REAL(KIND=8),DIMENSION(NDS)::DU
	KM=MAT(IE)
	KP=KGSS(KM,2)
	IF(KP/=0)THEN
		DO J=1,NDS
			L=NOD(IE,J)
			DU(J)=DTWH(L)
		END DO
		PP=0.0
		DO L=1,NDS
			PP=PP+DU(L)
		END DO
		PP=GMAW*PP/REAL(NDS)
		PWS(IE)=PWS(IE)+PP
	ELSE
		PWS(IE)=0.0
	END IF
	RETURN
END SUBROUTINE WPCALC


SUBROUTINE WETDFS(IE,LM,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,J,K,L
	REAL(KIND=8),DIMENSION(6)::SS,SN,DS
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF
	REAL(KIND=8),DIMENSION(6,6)::DD
	REAL(KIND=8)::SG3,PP,QQ,SL,WB,WC,WD,BN,CN,DN,VS,RS,DLT1,DLT2,SK
	DO J=1,NDS*NDM
		FF(J)=0.0
	END DO
	DO K=1,6
		SN(K)=0.0
	END DO
	DO L=1,6
		SS(L)=STS(IE,L)
	END DO
	WB=CW(LM,1)
	WC=CW(LM,2)
	WD=CW(LM,3)
	BN=CW(LM,4)
	CN=CW(LM,5)
	DN=CW(LM,6)
	SG3=STS(IE,9)
	PP=STS(IE,10)
	QQ=STS(IE,11)
	SL=STS(IE,12)
	IF(SG3<PAR) SG3=PAR
	IF(SL>0.90) SL=0.90
	SL=SL/(1.0-SL)
	VS=WB*(SG3/PAR)**BN+WC*(QQ/PAR)**CN
	RS=WD*SL**DN
	DO L=1,3
		SS(L)=SS(L)-PP
	END DO
	DO K=1,6
		SK=SS(K)
		DLT1=1.0
		DLT2=1.0
		IF(K>3)THEN
			DLT1=0.0
			DLT2=2.0
		END IF
		SN(K)=(1.0/3.0)*VS*DLT1+(SK/QQ)*RS*DLT2
	END DO
	LT=KGSS(LM,1)
	CALL EVSNCALC(IE)
	CALL DMAT(IE,LM,LT,DD)
	DO J=1,6
		DS(J)=0.0
		DO K=1,6
			DS(J)=DS(J)+DD(J,K)*SN(K)
		END DO
	END DO
	CALL NDFRCI(IE,DS,FF)
	RETURN
END SUBROUTINE WETDFS


SUBROUTINE WETDFD(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::I,J,K,IE

	RETURN
END SUBROUTINE WETDFD


SUBROUTINE TRNSFMS(A,B,C)
	IMPLICIT NONE
	INTEGER::K,L
	REAL(KIND=8)::I1,I2,I3,J2,J3,PP,QQ
	REAL(KIND=8)::FF,LO,LS,DLT,RF,LI1,LI2,LI3
	REAL(KIND=8),DIMENSION(6,6)::C
	REAL(KIND=8),DIMENSION(6)::A,B,DI1,DI2,DI3,DJ2,DJ3,DL
	CALL INVARS(A,I1,I2,I3,J2,J3,PP,QQ)
	CALL DSGIJ(A,DI1,DI2,DI3,DJ2,DJ3)
	FF=(I1*I2+I3)/(I1*I2+9.0*I3)
	FF=SQRT(FF)
	RF=((3.0*FF-1.0)*(I1*I2+9.0*I3))**2
	RF=24.0/(FF*RF)
	LO=SQRT(2.0/3.0)*2.0*I1/(3.0*FF-1.0)
	LS=SQRT(2.0/3.0)*QQ
	LI1=SQRT(2.0/3.0)*(2.0/(3.0*FF-1.0)-RF*I1*I2*I3)
	LI2=-SQRT(2.0/3.0)*RF*I1*I1*I3
	LI3=SQRT(2.0/3.0)*RF*I1*I1*I2
	DO L=1,6
		DL(L)=LI1*DI1(L)+LI2*DI2(L)+LI3*DI3(L)
	END DO
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=0.0
		A(L)=A(L)-PP*DLT
	END DO
	DO K=1,6
		DO L=1,6
			C(K,L)=0.0
		END DO
	END DO
	DO K=1,3
		DO L=1,3
			C(K,L)=(1.0-LO/LS)/3.0
		END DO
	END DO
	DO K=1,6
		DLT=1.0
		C(K,K)=C(K,K)+DLT*LO/LS
	END DO
	DO K=1,6
		DO L=1,6
			DLT=1.0
			IF(L>3) DLT=2.0
			C(K,L)=C(K,L)-DLT*A(K)*A(L)*LO/(LS**3)
		END DO
	END DO
	DO K=1,6
		DO L=1,6
			C(K,L)=C(K,L)+A(K)*DL(L)/LS
		END DO
	END DO
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=0.0
		B(L)=A(L)*LO/LS+PP*DLT
	END DO
	RETURN
END SUBROUTINE TRNSFMS


SUBROUTINE SGMTS(N,A,B,C)
	IMPLICIT NONE		!A-MAIN;B-LMN;C-NEW	[C]=[B][A][BT]
	INTEGER::I,J,K,N
	REAL(KIND=8),DIMENSION(N,N)::A,B,C,D
	DO I=1,N
		DO J=1,N
			D(I,J)=0.0
			DO K=1,N
				D(I,J)=D(I,J)+B(I,K)*A(K,J)
			END DO
		END DO
	END DO
	DO I=1,N
		DO J=1,N
			C(I,J)=0.0
			DO K=1,N
				C(I,J)=C(I,J)+D(I,K)*B(J,K)
			END DO
		END DO
	END DO
	RETURN
END SUBROUTINE SGMTS


SUBROUTINE LDING(IE)
	USE COMM
	INTEGER::IE,KM,KT,KL
	KM=MAT(IE)
	KT=KGSS(KM,1)
	KL=KGSS(KM,7)
	SELECT CASE(KL)
		CASE(0)
			CALL LDING0(IE)
		CASE(1)
			CALL LDING1(IE)
		CASE(2)
			CALL LDING2(IE)
		CASE(11)
			CALL LDING11(IE)
		CASE(-1)
			CALL JNTLDING1(IE)
		CASE DEFAULT
			WRITE(40,*) 'ERROR OCCUR IN SUBROUTINE LDING(IE):'
			WRITE(40,*) 'IE=',IE
			STOP
	END SELECT
	RETURN
END SUBROUTINE LDING


SUBROUTINE LDING0(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE
	LDG(IE)=1
	PLS(IE,1)=0.0
	PLS(IE,2)=0.0
	RETURN
END SUBROUTINE LDING0


SUBROUTINE LDING1(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,L,i,j
	REAL(KIND=8)::QP,i1,i2,i3,j2,j3,pp,qq
	REAL(KIND=8),DIMENSION(9)::ST
	REAL(KIND=8),DIMENSION(NDM,NDM)::A,B
	REAL(KIND=8),DIMENSION(6)::SS,DS,ST0,rt
	LM=MAT(IE)
	LT=KGSS(LM,1)
	DO L=1,6
		SS(L)=STS(IE,L)
		DS(L)=DSTS(IE,L)
		ST(L)=SS(L)+DS(L)
		ST0(L)=STR(IE,L)
		SS(L)=ST(L)
	END DO
	A(1,1)=ST(1)
	A(2,2)=ST(2)
	A(3,3)=ST(3)
	A(1,2)=ST(4)
	A(1,3)=ST(6)
	A(2,1)=ST(4)
	A(2,3)=ST(5)
	A(3,1)=ST(6)
	A(3,2)=ST(5)
	CALL MNSTS(NDM,A,B)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,6
     rt(i)=st(i)
	end do
    !call stsmod(ie,a,rt)
	do i=1,6
     st(i)=rt(i)
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ST(7)=A(1,1)
	ST(8)=A(2,2)
	ST(9)=A(3,3)
    forall(j=1:6)ss(j)=st(j)
    CALL INVARS(ss,I1,I2,I3,J2,J3,PP,QQ)
	CALL STSLVL(LM,LT,ST,ST0,QP)
	IF(QQ<=0.95*YLDM(IE,1).AND.QP<=0.95*YLDM(IE,2))THEN
		LDG(IE)=-1
	ELSE
		LDG(IE)=1
	END IF
	RETURN
END SUBROUTINE LDING1


SUBROUTINE LDING2(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,L,i
	REAL(KIND=8)::QP,DF,SG3
	REAL(KIND=8),DIMENSION(9)::ST
	REAL(KIND=8),DIMENSION(NDM,NDM)::A,B
	REAL(KIND=8),DIMENSION(6)::SS,DS,ST0,rt
	LM=MAT(IE)
	LT=KGSS(LM,1)
	DO L=1,6
		SS(L)=STS(IE,L)
		DS(L)=DSTS(IE,L)
		ST(L)=SS(L)+DS(L)
		ST0(L)=STR(IE,L)
		SS(L)=ST(L)
	END DO
	A(1,1)=ST(1)
	A(2,2)=ST(2)
	A(3,3)=ST(3)
	A(1,2)=ST(4)
	A(1,3)=ST(6)
	A(2,1)=ST(4)
	A(2,3)=ST(5)
	A(3,1)=ST(6)
	A(3,2)=ST(5)
	CALL MNSTS(NDM,A,B)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,6
     rt(i)=st(i)
	end do
   ! call stsmod(ie,a,rt)
	do i=1,6
     st(i)=rt(i)
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ST(7)=A(1,1)
	ST(8)=A(2,2)
	ST(9)=A(3,3)
	CALL STSLVL(LM,LT,ST,ST0,QP)
	SG3=(st(7)+st(8)+st(9))/3.0
	IF(SG3<PAR) SG3=PAR
	DF=QP*(SG3/PAR)**0.25
	IF(DF<0.95*YLDM(IE,1).AND.QP<=0.95*YLDM(IE,2))THEN
		LDG(IE)=-1
	ELSE
		LDG(IE)=1
	END IF
	RETURN
END SUBROUTINE LDING2


SUBROUTINE YLDMCALC(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,KT,KL,L
	REAL(KIND=8)::YLD1,YLD2,I1,I2,I3,J2,J3,PP,QQ,QP
	REAL(KIND=8),DIMENSION(3)::SG
	REAL(KIND=8),DIMENSION(6)::SS,SS0
	REAL(KIND=8),DIMENSION(NDM,NDM)::A,B
	KM=MAT(IE)
	KT=KGSS(KM,1)
	KL=KGSS(KM,7)
	DO L=1,6
		SS(L)=STS(IE,L)
		SS0(L)=STR(IE,L)
	END DO
	SELECT CASE(KL)
		CASE(0)
			RETURN
		CASE(1)
			CALL YLDM1(ie,KM,KT,SS,SS0,PP,QP)
			YLD(IE,1)=PP
			YLD(IE,2)=QP
			IF(YLDM(IE,1)<PP) YLDM(IE,1)=PP
			IF(YLDM(IE,2)<QP) YLDM(IE,2)=QP
		CASE(2)
			CALL YLDM2(ie,KM,KT,SS,SS0,YLD1,YLD2)
			YLD(IE,1)=YLD1
			YLD(IE,2)=YLD2
			IF(YLDM(IE,1)<YLD1) YLDM(IE,1)=YLD1
			IF(YLDM(IE,2)<YLD2) YLDM(IE,2)=YLD2
		CASE(11)
			IF(KT==11)THEN
				CALL YLDM11(KM,SS,SS0,YLD1,YLD2)
				YLD(IE,1)=YLD1
				YLD(IE,2)=YLD2
				IF(YLDM(IE,1)<YLD1) YLDM(IE,1)=YLD1
				IF(YLDM(IE,2)<YLD2) YLDM(IE,2)=YLD2
			END IF
			IF(KT==12)THEN
				CALL YLDM12(KM,SS,SS0,YLD1,YLD2)
				YLD(IE,1)=YLD1
				YLD(IE,2)=YLD2
				IF(YLDM(IE,1)<YLD1) YLDM(IE,1)=YLD1
				IF(YLDM(IE,2)<YLD2) YLDM(IE,2)=YLD2
			END IF
			IF(KT==13)THEN
				CALL YLDM13(KM,SS,SS0,YLD1)
				YLD(IE,1)=YLD1
				IF(YLDM(IE,1)<YLD1) YLDM(IE,1)=YLD1
			END IF
			IF(KT==14)THEN
				CALL YLDM14(KM,SS,SS0,YLD1)
				YLD(IE,1)=YLD1
				IF(YLDM(IE,1)<YLD1) YLDM(IE,1)=YLD1
			END IF
		CASE DEFAULT
			RETURN
	END SELECT
	RETURN
END SUBROUTINE YLDMCALC


SUBROUTINE YLDM1(ie,KM,KT,SS,SS0,PP,QP)
	USE COMM
	IMPLICIT NONE
	INTEGER::KM,KT,L,i,ie,j
	REAL(KIND=8)::PP,QP,qq,i1,i2,i3,j2,j3
	REAL(KIND=8),DIMENSION(9)::ST
	REAL(KIND=8),DIMENSION(6)::SS,SS0,rt
	REAL(KIND=8),DIMENSION(NDM,NDM)::A,B
	A(1,1)=SS(1)
	A(2,2)=SS(2)
	A(3,3)=SS(3)
	A(1,2)=SS(4)
	A(1,3)=SS(6)
	A(2,1)=SS(4)
	A(2,3)=SS(5)
	A(3,1)=SS(6)
	A(3,2)=SS(5)
	CALL MNSTS(NDM,A,B)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,6
     rt(i)=ss(i)
	end do
   ! call stsmod(ie,a,rt)
	do i=1,6
     ss(i)=rt(i)
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO L=1,6
		ST(L)=SS(L)
	END DO
	ST(7)=A(1,1)
	ST(8)=A(2,2)
	ST(9)=A(3,3)
    forall(j=1:6)ss(j)=st(j)
	call INVARS(SS,I1,I2,I3,J2,J3,qq,pp)
	CALL STSLVL(KM,KT,ST,SS0,QP)
	RETURN
END SUBROUTINE YLDM1


SUBROUTINE YLDM2(ie,KM,KT,SS,SS0,DF,QP)
	USE COMM
	IMPLICIT NONE
	INTEGER::KM,KT,L,i,ie
	REAL(KIND=8)::QP,SG3,DF
	REAL(KIND=8),DIMENSION(9)::ST
	REAL(KIND=8),DIMENSION(6)::SS,SS0,rt
	REAL(KIND=8),DIMENSION(NDM,NDM)::A,B
	A(1,1)=SS(1)
	A(2,2)=SS(2)
	A(3,3)=SS(3)
	A(1,2)=SS(4)
	A(1,3)=SS(6)
	A(2,1)=SS(4)
	A(2,3)=SS(5)
	A(3,1)=SS(6)
	A(3,2)=SS(5)
	CALL MNSTS(NDM,A,B)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,6
     rt(i)=ss(i)
	end do
   ! call stsmod(ie,a,rt)
	do i=1,6
     ss(i)=rt(i)
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO L=1,6
		ST(L)=SS(L)
	END DO
	ST(7)=A(1,1)
	ST(8)=A(2,2)
	ST(9)=A(3,3)
	CALL STSLVL(KM,KT,ST,SS0,QP)
	SG3=(ST(7)+ST(8)+ST(9))/3.0
	IF(SG3<PAR) SG3=PAR
	DF=QP*(SG3/PAR)**0.25
	RETURN
END SUBROUTINE YLDM2


SUBROUTINE YLDM11(KM,SS,SS0,YLD1,YLD2)
	USE COMM
	IMPLICIT NONE
	INTEGER::KM
	REAL(KIND=8),DIMENSION(6)::SS,SS0
	REAL(KIND=8)::AR,AS,PP,QQ,YLD1,YLD2
	REAL(KIND=8)::I1,I2,I3,J2,J3
	AR=2.0
	AS=2.0
	CALL INVARS(SS,I1,I2,I3,J2,J3,PP,QQ)
	YLD1=PP*PP+AR*AR*QQ*QQ
	YLD2=QQ**AS/PP
	RETURN
END SUBROUTINE YLDM11


SUBROUTINE YLDM12(KM,SS,SS0,YLD1,YLD2)
	USE COMM
	IMPLICIT NONE
	INTEGER::KM
	REAL(KIND=8),DIMENSION(6)::SS,SS0
	REAL(KIND=8)::AR,AS,PP,QQ,YLD1,YLD2
	REAL(KIND=8)::I1,I2,I3,J2,J3
	YLD1=0.0
	YLD2=0.0
	RETURN
END SUBROUTINE YLDM12


SUBROUTINE YLDM13(KM,ST,ST0,YLD1)
	USE COMM
	IMPLICIT NONE
	INTEGER::KM,L
	REAL(KIND=8)::BL,BK,BE,BM,BC,BF,S0,YLD1,CP,	&
					I1,I2,I3,J2,J3,PP0,QQ0,PP,QQ
	REAL(KIND=8),DIMENSION(6,6)::TP
	REAL(KIND=8),DIMENSION(6)::ST,ST0,SS,SS0
	BL=CS(KM,1)
	BK=CS(KM,2)
	BE=CS(KM,3)
	BM=CS(KM,4)
	BF=CS(KM,12)/RAD
	BC=CS(KM,14)
	CP=(BL-BK)/(1.0+BE)
	S0=BC/TAN(BF)
	DO L=1,3
		ST(L)=ST(L)+S0
		ST0(L)=ST0(L)+S0
	END DO
	CALL TRNSFMS(ST0,SS0,TP)
	CALL TRNSFMS(ST,SS,TP)
	PP0=(SS0(1)+SS0(2)+SS0(3))/3.0
	CALL INVARS(SS,I1,I2,I3,J2,J3,PP,QQ)
	YLD1=LOG(PP/PP0)+LOG(1.0+QQ*QQ/(BM*BM*PP*PP))
	YLD1=YLD1*CP
	RETURN
END SUBROUTINE YLDM13


SUBROUTINE YLDM14(KM,ST,ST0,YLD1)
	USE COMM
	IMPLICIT NONE
	INTEGER::KM,L
	REAL(KIND=8)::BL,BK,BE,BM,ALF,BC,BF,S0,YLD1,CP,	&
					I1,I2,I3,J2,J3,PP0,QQ0,PP,QQ,DLT
	REAL(KIND=8),DIMENSION(6,6)::TP
	REAL(KIND=8),DIMENSION(6)::ST,ST0,SS,SS0,YT
	BL=CS(KM,1)
	BK=CS(KM,2)
	BE=CS(KM,3)
	BM=CS(KM,4)
	ALF=CS(KM,7)
	BF=CS(KM,12)/RAD
	BC=CS(KM,14)
	CP=(BL-BK)/(1.0+BE)
	S0=BC/TAN(BF)
	DO L=1,3
		ST(L)=ST(L)+S0
		ST0(L)=ST0(L)+S0
	END DO
	CALL TRNSFMS(ST0,SS0,TP)
	CALL TRNSFMS(ST,SS,TP)
	PP=(SS(1)+SS(2)+SS(3))/3.0
	PP0=(SS0(1)+SS0(2)+SS0(3))/3.0
	DO L=1,6
		DLT=1.0
		IF(L>3) DLT=0.0
		SS(L)=SS(L)/PP-DLT
		SS0(L)=SS0(L)/PP0-DLT
	END DO
	DO L=1,6
		YT(L)=SS(L)-ALF*SS0(L)
	END DO
	CALL INVARS(YT,I1,I2,I3,J2,J3,YLD1,QQ)
	YLD1=LOG(PP/PP0)+QQ/BM
	YLD1=YLD1*CP
	RETURN
END SUBROUTINE YLDM14


SUBROUTINE LDING11(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,L
	REAL(KIND=8),DIMENSION(6)::SS,SS0,DS
	REAL(KIND=8)::A1,A2
	LM=MAT(IE)
	LT=KGSS(LM,1)
	A1=1.0
	A2=1.0
	LDG(IE)=1
	DO L=1,6
		SS(L)=STS(IE,L)
		DS(L)=DSTS(IE,L)
		SS0(L)=STR(IE,L)
	END DO
	IF(LT==11.OR.LT==12)THEN
		CALL BISECTD(IE,LM,SS,SS0,DS,A1,A2)
		PLS(IE,1)=A1
		PLS(IE,2)=A2
	END IF
	IF(LT==13.OR.LT==14)THEN
		CALL BISECTS(IE,LM,SS,SS0,DS,A1)
		PLS(IE,1)=A1
	END IF
	RETURN
END SUBROUTINE LDING11


SUBROUTINE BISECTD(IE,LM,SS,SS0,DS,A1,A2)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,L,LD1,LD2
	REAL(KIND=8)::EPS,A1,A2,YLDA1,YLDB1,YLDO1,YLDA2,YLDB2,YLDO2
	REAL(KIND=8),DIMENSION(6)::SS,SS0,DS,RR
	LT=KGSS(LM,1)
	EPS=1.0E-5
	YLDA1=YLD(IE,1)
	YLDO1=YLDM(IE,1)
	YLDA2=YLD(IE,2)
	YLDO2=YLDM(IE,2)
	DO L=1,6
		RR(L)=SS(L)+DS(L)
	END DO
	IF(LT==11) CALL YLDM11(LM,RR,SS0,YLDB1,YLDB2)
	IF(LT==12) CALL YLDM12(LM,RR,SS0,YLDB1,YLDB2)
	IF(ABS(YLDO1-YLDA1)<=EPS)THEN
		IF(YLDB1>YLDO1)THEN
			A1=0.0
			LD1=1
		ELSE
			A1=1.0
			LD1=-1
		END IF
	ELSE
		IF(YLDB1<YLDO1)THEN
			A1=1.0
			LD1=-1
		ELSE
			A1=(YLDO1-YLDA1)/(YLDB1-YLDA1)
			LD1=1
		END IF
	END IF
	IF(ABS(YLDO2-YLDA2)<=EPS)THEN
		IF(YLDB2>YLDO2)THEN
			A2=0.0
			LD2=1
		ELSE
			A2=1.0
			LD2=-1
		END IF
	ELSE
		IF(YLDB2<YLDO2)THEN
			A2=1.0
			LD2=-1
		ELSE
			A2=(YLDO2-YLDA2)/(YLDB2-YLDA2)
			LD2=1
		END IF
	END IF
	IF(A1>1.0) A1=1.0
	IF(A1<0.0) A1=0.0
	IF(A2>1.0) A2=1.0
	IF(A2<0.0) A2=0.0
	IF(LD1==-1.AND.LD2==-1) LDG(IE)=-1
	RETURN
END SUBROUTINE BISECTD


SUBROUTINE BISECTS(IE,LM,SS,SS0,DS,A1)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,L,LD1
	REAL(KIND=8)::EPS,A1,YLDA1,YLDB1,YLDO1
	REAL(KIND=8),DIMENSION(6)::SS,SS0,DS,RR
	LT=KGSS(LM,1)
	EPS=1.0E-5
	YLDA1=YLD(IE,1)
	YLDO1=YLDM(IE,1)
	DO L=1,6
		RR(L)=SS(L)+DS(L)
	END DO
	IF(LT==13) CALL YLDM13(LM,RR,SS0,YLDB1)
	IF(LT==14) CALL YLDM14(LM,RR,SS0,YLDB1)
	IF(ABS(YLDO1-YLDA1)<=EPS)THEN
		IF(YLDB1>YLDO1)THEN
			A1=0.0
			LD1=1
		ELSE
			A1=1.0
			LD1=-1
		END IF
	ELSE
		IF(YLDB1<YLDO1)THEN
			A1=1.0
			LD1=-1
		ELSE
			A1=(YLDO1-YLDA1)/(YLDB1-YLDA1)
			LD1=1
		END IF
	END IF
	IF(A1>1.0) A1=1.0
	IF(A1<0.0) A1=0.0
	IF(LD1==-1) LDG(IE)=-1
	RETURN
END SUBROUTINE BISECTS


SUBROUTINE MATRTRANI(IE,A,S)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,K1,K2,K3,K4,j
	REAL(KIND=8)::DX1,DY1,DZ1,DX2,DY2,DZ2,DX,DY,DZ,DD,S
	REAL(KIND=8),DIMENSION(NDM,NDM)::A
!	if(ie==8822)print*,(nod(ie,j),j=1,4)
	K1=NOD(IE,1)
	K2=NOD(IE,2)
	K3=NOD(IE,3)
	K4=NOD(IE,4)
!	if(ie==8822) print*,k1,k2
!	if(ie.eq.18641)print*,NOD(IE,1),NOD(IE,2),NOD(IE,3),NOD(IE,4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if((k1.eq.k2.and.k3.eq.k4).or.(k1.eq.k3.and.k2.eq.k4).or.(k1.eq.k4.and.k2.eq.k3))then
	K1=NOD(IE,5)
	K2=NOD(IE,6)
	K3=NOD(IE,7)
	K4=NOD(IE,8)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DX1=COR(K3,1)-COR(K1,1)
	DY1=COR(K3,2)-COR(K1,2)
	DZ1=COR(K3,3)-COR(K1,3)
	DX2=COR(K2,1)-COR(K4,1)
	DY2=COR(K2,2)-COR(K4,2)
	DZ2=COR(K2,3)-COR(K4,3)
	DX=DY1*DZ2-DZ1*DY2
	DY=DX2*DZ1-DZ2*DX1
	DZ=DX1*DY2-DY1*DX2
	DD=SQRT(DX**2+DY**2+DZ**2)
	!if(ie.eq.18641)print*,dx,dy,dz
	A(3,1)=DX/DD
	A(3,2)=DY/DD
	A(3,3)=DZ/DD
	DX=0.5*(COR(K3,1)+COR(K4,1)-COR(K2,1)-COR(K1,1))
	DY=0.5*(COR(K3,2)+COR(K4,2)-COR(K2,2)-COR(K1,2))
	DZ=0.5*(COR(K3,3)+COR(K4,3)-COR(K2,3)-COR(K1,3))
	DD=SQRT(DX**2+DY**2+DZ**2)
	A(1,1)=DX/DD
	A(1,2)=DY/DD
	A(1,3)=DZ/DD
	DX=A(1,3)*A(3,2)-A(1,2)*A(3,3)
	DY=A(1,1)*A(3,3)-A(1,3)*A(3,1)
	DZ=A(1,2)*A(3,1)-A(1,1)*A(3,2)
	DD=SQRT(DX**2+DY**2+DZ**2)
	A(2,1)=DX/DD
	A(2,2)=DY/DD
	A(2,3)=DZ/DD
	S=0.0

	if((k1.ne.k4).and.(k2.ne.k3).and.(k3.ne.k4).and.(k1.ne.k2))then
	!	print*,ie
	DX1=SQRT((COR(K1,1)-COR(K2,1))**2+(COR(K1,2)-COR(K2,2))**2+(COR(K1,3)-COR(K2,3))**2)
	DY1=SQRT((COR(K3,1)-COR(K1,1))**2+(COR(K3,2)-COR(K1,2))**2+(COR(K3,3)-COR(K1,3))**2)
	DZ1=SQRT((COR(K2,1)-COR(K3,1))**2+(COR(K2,2)-COR(K3,2))**2+(COR(K2,3)-COR(K3,3))**2)
	DD=(DX1**2+DY1**2-DZ1**2)/(2.0*DX1*DY1)
	DD=ACOS(DD)
	S=S+0.5*DX1*DY1*SIN(DD)
	DX2=SQRT((COR(K1,1)-COR(K3,1))**2+(COR(K1,2)-COR(K3,2))**2+(COR(K1,3)-COR(K3,3))**2)
	DY2=SQRT((COR(K1,1)-COR(K4,1))**2+(COR(K1,2)-COR(K4,2))**2+(COR(K1,3)-COR(K4,3))**2)
	DZ2=SQRT((COR(K3,1)-COR(K4,1))**2+(COR(K3,2)-COR(K4,2))**2+(COR(K3,3)-COR(K4,3))**2)
	DD=(DX2**2+DY2**2-DZ2**2)/(2.0*DX2*DY2)
	DD=ACOS(DD)
	S=S+0.5*DX2*DY2*SIN(DD)
	elseif((k1.eq.k4).or.(k3.eq.k4))then
	!	print*,ie
	DX1=SQRT((COR(K1,1)-COR(K2,1))**2+(COR(K1,2)-COR(K2,2))**2+(COR(K1,3)-COR(K2,3))**2)
	DY1=SQRT((COR(K3,1)-COR(K1,1))**2+(COR(K3,2)-COR(K1,2))**2+(COR(K3,3)-COR(K1,3))**2)
	DZ1=SQRT((COR(K2,1)-COR(K3,1))**2+(COR(K2,2)-COR(K3,2))**2+(COR(K2,3)-COR(K3,3))**2)
	DD=(DX1**2+DY1**2-DZ1**2)/(2.0*DX1*DY1)
	DD=ACOS(DD)
	S=S+0.5*DX1*DY1*SIN(DD)
	elseif(k2.eq.k3)then
	!	print*,ie
	DX1=SQRT((COR(K1,1)-COR(K2,1))**2+(COR(K1,2)-COR(K2,2))**2+(COR(K1,3)-COR(K2,3))**2)
	DY1=SQRT((COR(K4,1)-COR(K1,1))**2+(COR(K4,2)-COR(K1,2))**2+(COR(K4,3)-COR(K1,3))**2)
	DZ1=SQRT((COR(K2,1)-COR(K4,1))**2+(COR(K2,2)-COR(K4,2))**2+(COR(K2,3)-COR(K4,3))**2)
	DD=(DX1**2+DY1**2-DZ1**2)/(2.0*DX1*DY1)
	DD=ACOS(DD)
	S=S+0.5*DX1*DY1*SIN(DD)
	elseif(k1.eq.k2)then
!	print*,ie
	DX1=SQRT((COR(K1,1)-COR(K3,1))**2+(COR(K1,2)-COR(K3,2))**2+(COR(K1,3)-COR(K3,3))**2)
	DY1=SQRT((COR(K4,1)-COR(K1,1))**2+(COR(K4,2)-COR(K1,2))**2+(COR(K4,3)-COR(K1,3))**2)
	DZ1=SQRT((COR(K3,1)-COR(K4,1))**2+(COR(K3,2)-COR(K4,2))**2+(COR(K3,3)-COR(K4,3))**2)
	DD=(DX1**2+DY1**2-DZ1**2)/(2.0*DX1*DY1)
	DD=ACOS(DD)
	S=S+0.5*DX1*DY1*SIN(DD)
	end if
!	print*,ie
	RETURN
END SUBROUTINE MATRTRANI


SUBROUTINE MATRTRANJ(IE,A,S)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,K1,K5
	REAL(KIND=8)::DX,DY,DZ,S
	REAL(KIND=8),DIMENSION(NDM,NDM)::A
	K1=NOD(IE,1)
	K5=NOD(IE,5)
	DX=COR(K5,1)-COR(K1,1)
	DY=COR(K5,2)-COR(K1,2)
	DZ=COR(K5,3)-COR(K1,3)
	S=SQRT(DX**2+DY**2+DZ**2)
	A(1,1)=DX/S
	A(2,2)=DY/S
	A(3,3)=DZ/S
	RETURN
END SUBROUTINE MATRTRANJ


SUBROUTINE SDKSJ3(IE,KM,KS)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,J,K,L1,L2,LJ,LK,JJ,KK,JK,IH,IL,k1,k2,k3,k4
	REAL(KIND=8)::B,BL
	REAL(KIND=8),DIMENSION(NDS)::A,a1
	REAL(KIND=8),DIMENSION(NDM,NDM)::DD,TR,RT,RR
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS*NDM)::KS
	DATA A/4.0,2.0,1.0,2.0,4.0,2.0,1.0,2.0/
	data a1/2.0,1.0,1.0,1.0,2.0,1.0,1.0,1.0/
    k1=nod(ie,1)
    k2=nod(ie,2)
    k3=nod(ie,3)
    k4=nod(ie,4)
	IH=NFS*3
	IL=NFS*3
	DO J=1,NDM
		DO K=1,NDM
			DD(J,K)=0.0
			TR(J,K)=0.0
			RT(J,K)=0.0
		END DO
	END DO
	DD(1,1)=EST(IE,1)
	DD(2,2)=EST(IE,2)
	DD(3,3)=EST(IE,3)
	DO J=1,NDM
		DO K=1,NDM
			DMT(IE,J,K)=DD(J,K)
		END DO
	END DO
	DO J=1,NDM
		DO K=1,NDM
			TR(J,K)=DMT(IE,NDM+J,NDM+K)
			RT(K,J)=TR(J,K)
		END DO
	END DO
	CALL MATRIX_M3(RT,NDM,NDM,DD,NDM,NDM,TR,NDM,NDM,RR)
	DO J=1,NDS
		JJ=(J-1)*NDM
		JK=NDS-J+1
		DO K=1,JK
			KK=(K-1)*NDM
			DO L1=1,NDM
				LJ=JJ+KK+L1
				DO L2=1,NDM
					LK=KK+L2
					B=1.0
					IF(LJ>IH.AND.LK<=IL) B=-1.0
					if(k1/=k4.and.k2/=k3.and.k3/=k4)then
					KS(LJ,LK)=A(J)*RR(L1,L2)*B
					else
                    KS(LJ,LK)=A1(J)*RR(L1,L2)*B
					endif
					KS(LK,LJ)=KS(LJ,LK)
				END DO
			END DO
		END DO
	END DO
					if(k1/=k4.and.k2/=k3.and.k3/=k4)then
					    BL=EST(IE,4)/36.0
					else
                    	BL=EST(IE,4)/24.0
					endif

	DO J=1,NDS*NDM
		DO K=1,NDS*NDM
			KS(J,K)=KS(J,K)*BL
		END DO
	END DO
	RETURN
END SUBROUTINE SDKSJ3


SUBROUTINE SDKSJ4(IE,KM,KS)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,KM,J,K,L
	REAL(KIND=8)::DLT,AE,AA,AL,BL,RR
	REAL(KIND=8),DIMENSION(NDM)::XLMN
	REAL(KIND=8),DIMENSION(NDS*NDM)::BB
	REAL(KIND=8),DIMENSION(NDS*NDM,NDS*NDM)::KS
	XLMN(1)=DMT(IE,4,4)
	XLMN(2)=DMT(IE,5,5)
	XLMN(3)=DMT(IE,6,6)
	AE=EST(IE,1)
	AA=EST(IE,3)
	AL=EST(IE,4)
	BL=1.0/REAL(NFS*NFS)
	DO J=1,NDS
		K=(J-1)*NDM
		DLT=1.0
		IF(J<=NFS) DLT=-1.0
		DO L=1,NDM
			BB(K+L)=DLT*XLMN(L)
		END DO
	END DO
	DO J=1,NDS*NDM
		DO K=1,NDS*NDM
			KS(J,K)=BB(J)*BB(K)
		END DO
	END DO
	RR=AE*AA*BL/AL
	DO K=1,NDS*NDM
		DO L=1,NDS*NDM
			KS(K,L)=KS(K,L)*RR
		END DO
	END DO
	RETURN
END SUBROUTINE SDKSJ4


SUBROUTINE JNTLDING1(IE)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE
	REAL(KIND=8)::UU,VV,DU,DV,DW
	DU=DSTN(IE,6)
	DV=DSTN(IE,5)
	DW=DSTN(IE,3)
	UU=STN(IE,6)
	VV=STN(IE,5)
	LDG(IE)=1
	IF(DW>0.0) LDG(IE)=-1
	PLS(IE,1)=1.0
	PLS(IE,2)=1.0
	IF(DU*UU<0.0) PLS(IE,1)=-1.0
	IF(DV*VV<0.0) PLS(IE,2)=-1.0
	RETURN
END SUBROUTINE JNTLDING1

SUBROUTINE RLGDF(IE,LM,DT,FF)
	USE COMM
	IMPLICIT NONE
	INTEGER::IE,LM,LT,J,K,L
	REAL(KIND=8),DIMENSION(6)::SS,SN,DS
	REAL(KIND=8),DIMENSION(NDS*NDM)::FF
	REAL(KIND=8),DIMENSION(6,6)::DD
	REAL(KIND=8)::DT,SG3,PP,QQ,SL,RB,RC,RD,BN,CN,DN,ALF,VS,RS,DLT1,DLT2,SK
	DO J=1,NDS*NDM
		FF(J)=0.0
	END DO
	DO K=1,6
		SN(K)=0.0
	END DO
	DO L=1,6
		SS(L)=STS(IE,L)
	END DO
	RB=CR(LM,1)
	RC=CR(LM,2)
	RD=CR(LM,3)
	BN=CR(LM,4)
	CN=CR(LM,5)
	DN=CR(LM,6)
	ALF=CR(LM,7)
	SG3=STS(IE,9)
	PP=STS(IE,10)
	QQ=STS(IE,11)
	SL=STS(IE,12)
	IF(SG3<PAR) SG3=PAR
	IF(SL>0.90) SL=0.90
	SL=SL/(1.0-SL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!沈珠江双曲线模型
vs=rb*(sg3/par)
rs=rd*sl
vs=alf/vs*(vs-rhs(ie,1))**2*dt
rs=alf/rs*(rs-rhs(ie,2))**2*dt
RHS(IE,1)=RHS(IE,1)+VS
RHS(IE,2)=RHS(IE,2)+RS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO L=1,3
		SS(L)=SS(L)-PP
	END DO
	DO K=1,6
		SK=SS(K)
		DLT1=1.0
		DLT2=1.0
		IF(K>3)THEN
			DLT1=0.0
			DLT2=2.0
		END IF
		SN(K)=(1.0/3.0)*VS*DLT1+(SK/QQ)*RS*DLT2
		rhstn(ie,k)=sn(k)                         !流变产生的应变
	END DO
	LT=KGSS(LM,1)
	CALL EVSNCALC(IE)
	CALL DMAT(IE,LM,LT,DD)
	DO J=1,6
		DS(J)=0.0
		DO K=1,6
			DS(J)=DS(J)+DD(J,K)*SN(K)
		END DO
	END DO
	CALL NDFRCI(IE,DS,FF)
	RETURN
END SUBROUTINE RLGDF


SUBROUTINE ACCDIS(IP,ACX,ACY,ACZ,AM)
	USE COMM
	IMPLICIT NONE
	INTEGER::IP,IE,ME,LE,J,K
	REAL(KIND=8)::ACX,ACY,ACZ,AM,AN,HM,HN,DH,ZC,ZO,AME
	ME=LQ(IP,2)
	DO IE=1,NE
		DO J=1,NDM
			EQK(IE,J)=0.0
		END DO
	END DO
	LE=LUS(1)
	K=NOD(LE,1)
	HM=COR(K,NDM)
	HN=COR(K,NDM)
	DO IE=1,ME
		LE=LUS(IE)
		DO J=1,NDS
			K=NOD(LE,J)
			IF(COR(K,NDM)>HM) HM=COR(K,NDM)
			IF(COR(K,NDM)<HN) HN=COR(K,NDM)
		END DO
	END DO
	DH=HM-HN
	IF(DH<=40.0)THEN
		DO IE=1,ME
			LE=LUS(IE)
			ZC=0.0
			DO J=1,NDS
				K=NOD(LE,J)
				ZC=ZC+COR(K,NDM)
			END DO
			ZC=ZC/REAL(NDS)
			AME=1.0+((AM-1.0)/DH)*(ZC-HN)
			EQK(LE,1)=AME*ACX
			EQK(LE,2)=AME*ACY
			EQK(LE,3)=AME*ACZ
		END DO
	ELSE
		ZO=HN+0.6*DH
		AN=1.0+(AM-1.0)/3.0
		DO IE=1,ME
			LE=LUS(IE)
			ZC=0.0
			DO J=1,NDS
				K=NOD(LE,J)
				ZC=ZC+COR(K,NDM)
			END DO
			ZC=ZC/REAL(NDS)
			IF(ZC<=ZO)THEN
				AME=1.0+((AN-1.0)/(0.6*DH))*(ZC-HN)
			ELSE
				AME=AN+((AM-AN)/(0.4*DH))*(ZC-ZO)
			END IF
			EQK(LE,1)=AME*ACX
			EQK(LE,2)=AME*ACY
			EQK(LE,3)=AME*ACZ
		END DO
	END IF
	RETURN
END SUBROUTINE ACCDIS