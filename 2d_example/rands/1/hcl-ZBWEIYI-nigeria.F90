!------------------------------------------------------------
!-----------------------剖面数据整理-------------------------
!-------------------------------------------------------------
DIMENSION::X(150000),Y(150000),Z(150000),KF(150000,8),L1(150000),L2(150000),KF2(150000,8),X2(150000),Y2(150000),Z2(150000)
DIMENSION::DISPX(150000),DISPY(150000),DISPZ(150000),NND(150000),ma(100000),mb(100000),noe(200000)
DIMENSION::SIGMA1(150000),SIGMA2(150000),SIGMA3(150000),STRLEL(150000)
DIMENSION::MF(150000,4),MN(150000),NEM(150000),II(4),NB(30),MC(150000,4)
DIMENSION::NOD(150000),S1(150000),S3(150000),SL(150000),NTY(150000)
DIMENSION::FISPX(150000),FISPY(150000),FISPZ(150000),HEAD(150000)
dimension::eledis(150000,4)
DIMENSION::ID(150000,6)
character(len=60)::file1
real(kind=8),dimension(150000)::xt,yt,zt
!-------------------------------------------------------------THE Z-CORORDATE OF PLANE IS BETWEEN FOLLOWING TWO NUMBER
REAL(KIND=8),PARAMETER::ZP1=29        !min
REAL(KIND=8),PARAMETER::ZP2=31        !max
real(kind=8)::ra,rs,DX,DY,DZ
!-------------------------------------------------------------
DIMENSION::DPA(50,150000,5),DPB(50,150000,5),DPC(50,150000,5),DP1(150000,4),DP2(150000,4) !!存放各级节点位移计算结果
DIMENSION::ADP1(150000,4),BDP1(150000,4),DPD(50,150000,5),CDP1(150000,4),disold(100000,4)
DIMENSION::STR(50,150000,12),SS1(150000,12),SS2(150000,12)!!存放单元应力计算结果
DIMENSION::MME(50)       !THE NUMBER OF ELEMENTS OF EACH MATERICAL
real(kind=8):: xpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(123,FILE='ELENOD.dat',STATUS='OLD')
OPEN(2,FILE='DASAC1.DAT',STATUS='OLD')
OPEN(3,FILE='DASAC2.DAT',STATUS='OLD') !!!!!!!!!!!!!!!!!!!!!!!!!!_____________2-正常；3-缝———————3-张拉+张；5、6切向
OPEN(11,FILE='EL-STS-j1.PLT')
OPEN(12,FILE='EL-DIS-j1.PLT')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!单元数
!节点数

NM=14      !material 总数
NNK=1      !MATERIAL NEED TO ANALYSIS
!________________用于两级荷载之间的增量计算
PAR=1   !!!!!!!!!!!!!!!!!!___________________0-全量；1-增量
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NF=23       !荷载分级数
HZ1=22
HZ2=2      !地应力计算完毕
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!读入单元、节点、材料、荷载分级

read(123,*)ne !,nx2
DO I=1,NE
READ(123,*)KW,(KF2(I,J),J=1,8),l1(i) 
if(L1(I).le.13)then
l1(i)=0
!if(L1(I).eq.1.or.L1(I).eq.3)then
!L1(i)=0
!else if(L1(I).eq.12.or.L1(I).eq.13.or.L1(I).eq.14)then
!l1(i)=0
!else if(l1(i).eq.12)then
!l1(i)=1
else
l1(i)=1
end if
END DO


read(123,*)nx2
DO I=1,NX2
READ(123,*)KW,z2(I),x2(I),y2(I) 
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ma=1
do i=1,nx2-1
    if(ma(i).eq.0)cycle  
  do j=i+1,nx2
	  if(ma(j).eq.0)cycle
	  if(abs(z2(i)-z2(j)).lt.1e-3.and.abs(x2(i)-x2(j)).lt.1e-3.and.abs(y2(i)-y2(j)).lt.1e-3)ma(j)=0
  end do
end do

nx=0
do i=1,nx2
  if(ma(i).ne.0)then
     nx=nx+1
	 z(nx)=z2(i)
	 x(nx)=x2(i)
	 y(nx)=y2(i)
  end if
end do

do i=1,ne
 do j=1,8
    do jj=1,nx
       k=kf2(i,j)
	   kf(i,j)=kf2(i,j)
        	  if(abs(z(jj)-z2(k)).lt.1e-3.and.abs(x(jj)-x2(k)).lt.1e-3.and.abs(y(jj)-y2(k)).lt.1e-3)then
                kf(i,j)=jj
				exit
              end if
	end do
 end do
end do

print*,nx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!读入各级节点位移
DO I=1,NF
   READ(2,*)LD
   DO J=1,NX2
   READ(2,*)KW,(DPA(I,J,K),K=1,4) 
	 do k=1,3
      DPA(I,J,K)=DPA(I,J,K)
	 end do
   END DO
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!! 读入各级单元应力
DO I=1,NF
   READ(3,*)LD
   DO J=1,NE
    READ(3,*)KW,(STR(I,J,K),K=1,12) 
        do k=1,11
          STR(I,J,K)=STR(I,J,K)/1e+3
		end do
   END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!选择需要分析的分级荷载
print*,'file-read finished!'
DO I=1,NF
  IF(I==HZ1)THEN
    DO J=1,NX2
        DO K=1,4
  		  DP1(J,K)=DPA(I,J,K)
           ADP1(J,K)=DPB(I,J,K)
  		  BDP1(J,K)=DPC(I,J,K)
           CDP1(J,K)=DPD(I,J,K)
		END DO
	END DO
    DO J=1,NE
        DO K=1,12
		 SS1(J,K)=STR(I,J,K)
		END DO
	END DO
   
  ELSEIF(I==HZ2)THEN
  print*,111111111
    DO J=1,NX2
        DO K=1,4
		 DP2(J,K)=DPA(I,J,K)
		END DO
	END DO
    DO J=1,NE
        DO K=1,12
		 SS2(J,K)=STR(I,J,K)
		END DO
	END DO
  END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!赋值
DO I=1,NE
   SIGMA1(I)=SS1(I,1) ! Z-STRESS
   SIGMA2(I)=SS1(I,2) ! X-STRESS
   SIGMA3(I)=SS1(I,3) ! Y-STRESS
   STRLEL(I)=abs(SS1(I,5)) ! THE STRESS LEVEL
END DO

IF(PAR==0)THEN
DO I=1,NX2
DISPZ(I)=DP1(I,1)
DISPX(I)=DP1(I,2)
DISPY(I)=DP1(I,3)
HEAD(I)=DP1(I,4)
END DO
ELSE
DO I=1,NX2
DISPZ(I)=DP1(I,1)-DP2(I,1) !+disold(i,1)
DISPX(I)=DP1(I,2)-DP2(I,2) !+disold(i,2)
DISPY(I)=DP1(I,3)-DP2(I,3) !+disold(i,3)
HEAD(I)=DP1(I,4)
END DO
END IF

!!!!!!!!!!!!!!!!!!3D TECPLOT!!!!!!!!!!!!!!!!!!!!  NODE DISPLACEMENT
DO I=1,NM
  DO J=1,NE
   IF(L1(J)==I)MME(I)=MME(I)+1
  END DO
END DO

do i=1,ne
forall(j=1:4)eledis(i,j)=0.0
DZ=0
DX=0
do j=1,8
eledis(i,1)=eledis(i,1)+dispz(kf2(i,j))/8.0
eledis(i,2)=eledis(i,2)+dispx(kf2(i,j))/8.0
eledis(i,3)=eledis(i,3)+dispy(kf2(i,j))/8.0
eledis(i,4)=eledis(i,4)+(head(kf2(i,j))-y(kf2(i,j)))/8.0
end do
end do
! print*,DISPY(558)
!------------------------------------------------------------- ELEMEMTS STRESS
write(11,*)'TITLE="Example=FE-Volume Brick Data"'
write(11,*)'VARIABLES="x","y","z","sx","sy","sz","S1","S2","S3","SLEVEL"'

DO I=1,NM
if(mme(i).eq.0)cycle
WRITE(11,*)'ZONE N=',NX,',E=',MME(I),',DATAPACKING=BLOCK,ZONETYPE=FEBRICK,VARLOCATION=([4-10]=CELLCENTERED)'
DO J=1,NX
  WRITE(11,'(F12.3)') x(J)
END DO
DO J=1,NX
  WRITE(11,'(F12.3)') z(J)
END DO
DO J=1,NX
  WRITE(11,'(F12.3)') y(J)
END DO
DO J=1,NE
  IF(L1(J)==I) then
     WRITE(11,'(F12.3)') abs(SS1(J,3)) !*1e3     !!!!!!!!!!!!!!坝轴
  end if
END DO
DO J=1,NE
  IF(L1(J)==I) then
     WRITE(11,'(F12.3)') abs(SS1(J,5)) !*1e3        水平向
  end if
END DO
DO J=1,NE
  if(l1(j).eq.i)WRITE(11,'(F12.3)') abs(SS1(J,6))   !*1e3  竖向
END DO
DO J=1,NE
  IF(L1(J)==I)WRITE(11,'(F12.3)') SS1(J,7)
END DO
DO J=1,NE
  IF(L1(J)==I)WRITE(11,'(F12.3)') SS1(J,8)
END DO
DO J=1,NE
  IF(L1(J)==I)WRITE(11,'(F12.3)') SS1(J,9)
END DO
DO J=1,NE
ss1(j,12)=sqrt(ss1(j,5)**2+ss1(j,6)**2)/tan(32.5*3.14/180)*abs(ss1(j,3))
if(ss1(j,12).gt.0.95)ss1(j,12)=0.95
  IF(L1(J)==I)WRITE(11,'(F12.3)') SS1(J,12)
END DO
DO J=1,NE
  IF(L1(J)==I)THEN
    WRITE(11,'(8(I8,4X))')(KF(J,K),K=1,8)
  END IF
END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(12,*)'TITLE="Example=FE-Volume Brick Data"'
write(12,*)'VARIABLES="x","y","z","disx","disy","disz","head"'
DO I=1,NM
if(mme(i).eq.0)cycle
WRITE(12,*)'ZONE N=',NX,',E=',MME(I),',DATAPACKING=BLOCK,ZONETYPE=FEBRICK,VARLOCATION=([4-10]=CELLCENTERED)'
DO J=1,NX
  WRITE(12,'(F12.3)') x(J)
END DO
DO J=1,NX
  WRITE(12,'(F12.3)') z(J)
END DO
DO J=1,NX
  WRITE(12,'(F12.3)') y(J)
END DO
DO J=1,NE
  IF(L1(J)==I) then
     WRITE(12,'(F12.3)') eledis(j,1)*1e3     !!!!!!!!!!!!!!坝轴
  end if
END DO
DO J=1,NE
  IF(L1(J)==I) then
     WRITE(12,'(F12.3)') eledis(j,2)*1e3
  end if
END DO
DO J=1,NE
  if(l1(j).eq.i)WRITE(12,'(F12.3)') eledis(j,3)*1e3
END DO

DO J=1,NE
  if(l1(j).eq.i)WRITE(12,'(F12.3)') eledis(j,4)*0.0
END DO


DO J=1,NE
  IF(L1(J)==I)THEN
    WRITE(12,'(8(I8,4X))')(KF(J,K),K=1,8)
  END IF
END DO
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------
END
!------------------------------------------------------------
