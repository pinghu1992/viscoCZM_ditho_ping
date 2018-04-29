      PROGRAM miniABQ
!     ****************************************************************************
!   
!     PROGRAM: miniABQ
!   
!     PURPOSE:  For cohesive elements only 
!         for which all DOFs are known.
!   
!     ****************************************************************************
      INCLUDE 'ABA_PARAM.INC'
!     PROGRAM TO TEST A UMAT BEFORE USING IT INSIDE ABAQUS
!
!     PARAMETERS FOR THE TYPE OF UMAT
!     3D SOLID: NDI:3 NSHR:3, 2D SHELL: NDI:2, NSHR:1
      PARAMETER (NDI=3,NSHR=3,NTENS=NDI+NSHR, NDOFEL=24, NRHS=1)
      PARAMETER (NSVARS=12,NNODE=8,MCRD=3,MLVARX=NDOFEL) 
      PARAMETER (NPROPS=8,MDLOAD=3,NPREDF=24)
      PARAMETER (LAYER=1, CELENT=1.0)
      !MNDLOAD,NPREDF IS NOT USED IN THIS UEL
!
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),PROPS(NPROPS),
     1 SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(3),
     3 JDLTYP(MDLOAD,1),ADLMAG(MDLOAD,1),DDLMAG(MDLOAD,1),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(1),JPROPS(1)
!
      CHARACTER*80 CMNAME
      CHARACTER*72 STRING72
!      PARAMETER (CMNAME = 'WHATEVER MATERIAL YOU WANT TO DO')
!
!     YOUR TESTER EXECUTABLE STATEMENTS START HERE
      WRITE(6,"(/'       NDI,      NSHR,    NPROPS,   NSVARS')")
      WRITE(6,"(4I10)") NDI, NSHR, NPROPS, NSVARS
!
!     NDI: # OF DIRECT COMPONENTS (11,...) OF DDSDDE, DDSDDT, AND DRPLDE
!     NSHR: # OF ENGINEERING SHEAR COMPONENTS (12,...) OF DDSDDE, DDSDDT, AND DRPLDE
!     NTENS = NDI + NSHR: SIZE OF THE STRESS OR STRAIN COMPONENT ARRAY
!
!     IN ABAQUS
!     UNIT=6 WRITES TO THE .DAT FILE
!     UNIT=7 WRITES TO THE .MSG FILE
!     WHEN DOING STEP STATIC GENERAL, 1ST CALL TO UMAT IS WITH ZERO STRAIN AND DSTRAIN
!     SO YOU HAVE TO BE CAREFUL NOT TO DIVIDE BY ZERO
!     WHEN DOING STEP LINEAR PERTURBATION, ALL CALLS TO UMAT ARE WITH ZERO DSTRAIN

!     FILE PROPS CONTAINS LINES FORM THE ABAQUS INP FILE AS FOLLOWS
!     *USER MATERIAL, CONSTANTS=12
!
! -------------------- MATERIAL PARAMETERS START HERE ------------------------
      OPEN(UNIT=4,FILE='PROPS.TXT',STATUS='OLD',ACTION='READ')
      DO J=1,1
        READ(4,*) PROPS((J-1)*8+1), PROPS((J-1)*8+2), PROPS((J-1)*8+3),
     1        PROPS((J-1)*8+4), PROPS((J-1)*8+5), PROPS((J-1)*8+6),
     2 PROPS((J-1)*8+7), PROPS((J-1)*8+8)
      ENDDO
!
      CLOSE(UNIT=4)
! ----------------------- MATERIAL PARAMETER END HERE ----------------------
!
! ---- ZERO INITIALIZATION ----
      CALL k_vector_zero(SVARS,NSVARS)
      CALL k_matrix_zero(AMATRX,NDOFEL,NDOFEL)
      CALL k_vector_zero(U,NDOFEL)
      CALL k_matrix_zero(DU,NDOFEL,1)
      

!     OPEN FILE FOR THE OUTPUT FIELD VARIABLES
      OPEN (UNIT = 7, FILE = "FIELD_HIST.CSV")
      
!
!
!
!     HERE YOU DEFINE THE STRAIN AND DELTA STRAIN YOU WANT TO PASS TO THE UMAT
      KSTEP = 1 !EXPERIMENT WITH THIS. READ ABAQUS DOCUMENTATION
      KINC = 2  !EXPERIMENT WITH THIS. READ ABAQUS DOCUMENTATION
      TIME = 0.D0
      DTIME = 0.01D0
      MAXINC = 1001
      on_shear1 = 10
      on_shear2 = 1     
      on_normal = 10
!
      !original coordinates
      COORDS(1,1) = 0.d0
      COORDS(1,2) = 1.d0
      COORDS(1,3) = 1.d0
      COORDS(1,4) = 0.d0
      COORDS(1,5) = 0.d0
      COORDS(1,6) = 1.d0
      COORDS(1,7) = 1.d0
      COORDS(1,8) = 0.d0
      COORDS(2,1) = 0.d0
      COORDS(2,2) = 0.d0
      COORDS(2,3) = 1.d0
      COORDS(2,4) = 1.d0
      COORDS(2,5) = 0.d0
      COORDS(2,6) = 0.d0
      COORDS(2,7) = 1.d0
      COORDS(2,8) = 1.d0
      COORDS(3,1) = 0.d0
      COORDS(3,2) = 0.d0
      COORDS(3,3) = 0.d0
      COORDS(3,4) = 0.d0
      COORDS(3,5) = 0.001d0
      COORDS(3,6) = 0.001d0
      COORDS(3,7) = 0.001d0
      COORDS(3,8) = 0.001d0
!
!      WRITE(7,*) 'Time,Disp1,Stress1,Disp2,Stress2,Disp3,Stress3,
!     1 Dispm,Stressm' 
!      WRITE(7,15) TIME(1),',', 0,',', 0,',', 0,',', 0,',', 0,',', 0,
!     1 ',', 0,',', 0
      WRITE(7,*) 'displacement,ex,ey,exy,e1,e2,e12' 
      WRITE(7,15) TIME(1),',', 0,',', 0,',', 0,',', 0,',', 0,',', 0     
!
!     LOOPING FOR UMAT CALLER
      DO KINC=2,MAXINC
!
      if (on_normal == 1) then
        DU(15,1) = 0.00005D0
        DU(18,1) = 0.00005D0
        DU(21,1) = 0.00005D0
        DU(24,1) = 0.00005D0
      endif
      if (on_shear1 == 1) then
        DU(13,1) = 0.00006D0
        DU(16,1) = 0.00006D0
        DU(19,1) = 0.00006D0
        DU(22,1) = 0.00006D0
      end if
      if (on_shear2 == 1) then      
        DU(14,1) = 0.00000D0
        DU(17,1) = 0.00000D0
        DU(20,1) = 0.00001D0
        DU(23,1) = 0.00001D0
      end if        
        !
      do i=1,ndofel
          u(i) = u(i) + du(i,1);
      end do
!  
!  
! 
   
      CALL UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1 PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2 DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3 PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     4 NJPROP, PERIOD,sep_m1,Tm1,sep_m2,Tm2,sep_m3,Tm3,sep_m,Tm,ex,ey,
     5 exy,e1,e2,e12)
        
      TIME(1) = TIME(1) + DTIME
!  
!  
!        WRITE(7,15) TIME(1),',', sep_m1,',',Tm1,',', sep_m2,',',Tm2,','
!     1   , sep_m3,',',Tm3,',', sep_m,',',Tm
!15    FORMAT (F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6
!     1 ,A,F12.6,A,F12.6)
      displacement = 0.00001D0 * (KINC-1)
        WRITE(7,15) displacement,',', ex,',',ey,',', exy,',',e1,','
     1   , e2,',',e12
15    FORMAT (F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6)
!
      END DO
!
!      WRITE(6,"(/'   ---------DDSDDE--------   ')")
!      DO I=1,NTENS
!        WRITE(6,2) (DDSDDE(I,J),J=1,NTENS)
!      ENDDO
      
!      WRITE(6,"(/'    -------- STRESS --------   ')")
!      WRITE(6,2) STRESS
!2     FORMAT(6F12.3)
!
      CLOSE(7)
!
!
      END PROGRAM miniABQ

