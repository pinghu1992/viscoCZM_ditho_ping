c 3D Bilinear CZM with coupling UEL for ABAQUS Code (Compatible with 8 node Brick Elements)
c By: Ping Hu & Ditho Ardiansyah Pulungan
c March 15th, 2018
c
c =====================================================================
      SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1 PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2 DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3 PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     4 NJPROP, PERIOD,xstress)
c
      INCLUDE 'ABA_PARAM.INC'
c     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
c
      DIMENSION ds1(4),ds2(4),dn(4),Trac(MCRD,NRHS),
     1 Trac_Jacob(MCRD,MCRD),R(MCRD,MCRD),coord_l(MCRD,NNODE),
     2 GP_coord(2),sf(4),B(MCRD,NDOFEL),co_de_m(3,4),
     3 B_t(NDOFEL,MCRD), Transformation_M(NDOFEL,NDOFEL),
     4 Transformation_M_T(NDOFEL,NDOFEL),temp1(MCRD,NDOFEL)
c
      DIMENSION stiff_l(NDOFEL,NDOFEL),temp2(NDOFEL,NDOFEL),
     1 stiff_g(NDOFEL,NDOFEL),residual_l(NDOFEL,NRHS),
     2 residual_g(NDOFEL,NRHS),aJacob_M(2,3),delu_loc_gp(mcrd),
     3 co_de(mcrd,nnode), disp(mcrd,nnode), disp_l(mcrd,nnode),
     4 dxdxi(2,2),dxidx(2,2), dNdxi(4,2), dNdx(4,2), uv(2,4), 
     5 uv_xy(2,2)
c
      DOUBLE PRECISION G_fn, G_ft, f_tn, f_tt, K_n, K_t,
     1 sep_n0, sep_nf, sep_t0, sep_tf, tmax1, tmax2, pmax, opn, 
     2 opt1,opt2,tmax,c_r,c_s,detdxdxi,eps_x,eps_y,eps_xy,up_eps_vm,
     3 up_sigma,bt_sigma,eps_11,eps_12,eps_21,eps_22,bt_eps_vm
c
c Define Inputs========================================================
c 
      G_fn=PROPS(1)
      G_ft=PROPS(2)
      f_tn=PROPS(3)
      f_tt=PROPS(4)
      K_n=PROPS(5)
      K_t=PROPS(6)
      up_sigma=PROPS(7)
      bt_sigma=PROPS(8)
      GP_n=4d0
c write out results for debugging
c      OPEN(unit=16,file='C:\Temp\writeout.txt',status='unknown')
c      OPEN(unit=20,file='C:\Temp\writeout1.txt',status='unknown')
      OPEN(unit=18,file='C:\Temp\writeout2.txt',status='unknown')
c     WRITE(16,*) G_fn,G_ft,f_tn,f_tt,K_n,K_t
c
c Initialize Matrices and Vectors======================================
c 
      call k_vector_zero(ds1,4)
      call k_vector_zero(ds2,4)
      call k_vector_zero(dn,4)
      call k_matrix_zero(Trac,mcrd,nrhs)
      call k_matrix_zero(Trac_Jacob,mcrd,mcrd)
      call k_matrix_zero(R,mcrd,mcrd)
      call k_matrix_zero(coord_l,mcrd,nnode)
      call k_vector_zero(GP_coord,2)
      call k_vector_zero(sf,4)
      call k_matrix_zero(Transformation_M,ndofel,ndofel)
      call k_matrix_zero(Transformation_M_T,ndofel,ndofel)
      call k_matrix_zero(B,mcrd,ndofel)
      call k_matrix_zero(B_t,ndofel,mcrd)
      call k_matrix_zero(temp1,mcrd,ndofel)
      call k_matrix_zero(stiff_l,ndofel,ndofel)
      call k_matrix_zero(temp2,ndofel,ndofel)
      call k_matrix_zero(stiff_g,ndofel,ndofel)
      call k_matrix_zero(residual_l,ndofel,nrhs)
      call k_matrix_zero(residual_g,ndofel,nrhs)
      call k_matrix_zero(aJacob_M,2,3)
      call k_matrix_zero(rhs,ndofel,nrhs)
      call k_matrix_zero(amatrx,ndofel,ndofel)
      call k_matrix_zero(co_de,mcrd,nnode)
      call k_matrix_zero(disp,mcrd,nnode)
      call k_matrix_zero(disp_l,mcrd,nnode)      
      a_Jacob=0.d0
c
c Determine Inputs to Cohesive Model===================================
c sep_n0 is damage initial separation in the normal direction==========
c
      sep_n0=f_tn/K_n
      sep_nf=2.d0*G_fn/f_tn
      sep_t0=f_tt/K_t
      sep_tf=2.d0*G_ft/f_tt
!     WRITE(16,*) sep_n0,sep_nf,sep_t0,sep_tf
c
c Do local computations================================================
c 
      do i = 1, mcrd
         do j = 1, nnode
            co_de(i,j)=coords(i,j)+U(3.d0*(j-1.d0)+i)
            disp(i,j)=U(3.d0*(j-1.d0)+i)
         end do
      end do
c
c Do Calculations at Gauss Points======================================
c
      do i = 1, GP_n
c
      call k_matrix_zero(aJacob_M,2,3)
c
      gpt = i
c
c get the local coordinate and local displacement======================
c
      call k_local_coordinates(gpt,co_de,R,coord_l,Transformation_M,
     & Transformation_M_T,a_Jacob,aJacob_M,coords,u,ndofel,nnode,
     & mcrd, SVARS,disp,disp_l)   
c
c Compute shear and normal local separation(opening displacments)======
c 
      do j = 1, 4
!         ds1(j)=coord_l(1,j+4)-coord_l(1,j)
!         ds2(j)=coord_l(2,j+4)-coord_l(2,j)
!         dn(j) =coord_l(3,j+4)-coord_l(3,j)
         ds1(j)=disp_l(1,j+4)-disp_l(1,j)
         ds2(j)=disp_l(2,j+4)-disp_l(2,j)
         dn(j) =disp_l(3,j+4)-disp_l(3,j)         
      end do
!     WRITE(16,*) coord_l
c
c Determine the values of the shape function at each Gauss Point/Newton integration points
c
         call k_shape_fun(i,sf)
c
         call k_vector_zero(delu_loc_gp,mcrd)
c
c Determine shear and normal separation (opening displamenets) at integration points==============================
c         
         do j = 1, 4
            delu_loc_gp(1)=delu_loc_gp(1)+ds1(j)*sf(j)
            delu_loc_gp(2)=delu_loc_gp(2)+ds2(j)*sf(j)
            delu_loc_gp(3)=delu_loc_gp(3)+ dn(j)*sf(j)
         end do
         
!        if (i .EQ. 1) then
!            WRITE(16,*) 'separation', delu_loc_gp
!        end if

         opt1=ABS(delu_loc_gp(1))
         opt2=ABS(delu_loc_gp(2))
         opn=delu_loc_gp(3)
!        WRITE(16,*) opt1,opt2,opn
c
c update the max separation of the current step to calculate interface damage==========
c
         if ((Svars(GP_n*(i-1.0)+1.0) .LT. opt1) .AND.
     & (opt1 .GE. 0.0)) then
            Svars(GP_n*(i-1.0)+1.0)=opt1
         end if
         if ((Svars(GP_n*(i-1.0)+2.0) .LT. opt2) .AND.
     & (opt2 .GE. 0.0)) then
            Svars(GP_n*(i-1.0)+2.0)=opt2
         end if
         if ((Svars(GP_n*(i-1.0)+3.0) .LT. opn) .AND.
     & (opn .GE. 0.0)) then
            Svars(GP_n*(i-1.0)+3.0)=opn
         end if
         tmax1=Svars(GP_n*(i-1.0)+1.0)
         tmax2=Svars(GP_n*(i-1.0)+2.0)
         pmax=Svars(GP_n*(i-1.0)+3.0)
!        tmax=sqrt(tmax1**2.0d0+tmax2**2.0d0)         
! ---------------------------------------------------------
c
c Determine Traction vector and tangent modulus matrix
c
c
c========================Do Strain Calculations at integration Points==================
c
      if (gpt .eq. 1) then
         c_r=-1.0d0
         c_s=-1.0d0
      elseif (gpt .eq. 2) then
         c_r= 1.0d0
         c_s=-1.0d0
      elseif (gpt .eq. 3) then
         c_r= 1.0d0
         c_s= 1.0d0
      elseif (gpt .eq. 4) then
         c_r=-1.0d0
         c_s= 1.0d0
      end if
c
c shape function derivative to natural coordinate=================
c
      dNdxi(1,1) =-0.25d0*(1-c_s)
      dNdxi(2,1) = 0.25d0*(1-c_s)
      dNdxi(3,1) = 0.25d0*(1+c_s)
      dNdxi(4,1) =-0.25d0*(1+c_s)
      dNdxi(1,2) =-0.25d0*(1-c_r)
      dNdxi(2,2) =-0.25d0*(1+c_r)
      dNdxi(3,2) = 0.25d0*(1+c_r)
      dNdxi(4,2) = 0.25d0*(1-c_r)
c
c get Jacobian matrix(local-natural coordinate)=========================
c
      do ii = 1,2
         do j = 1,2
            do k =1, 4
               dxdxi(ii,j) = dxdxi(ii,j) + dNdxi(k,j)*coord_l(ii,k)
            end do
         end do
      end do
      dxdxi(1,1) = dNdxi(1,1)*coord_l(1,1)+dNdxi(2,1)*coord_l(1,2)+
     & dNdxi(3,1)*coord_l(1,3)+dNdxi(4,1)*coord_l(1,4)
      dxdxi(1,2) = dNdxi(1,1)*coord_l(2,1)+dNdxi(2,1)*coord_l(2,2)+
     & dNdxi(3,1)*coord_l(2,3)+dNdxi(4,1)*coord_l(2,4)
      dxdxi(2,1) = dNdxi(1,2)*coord_l(1,1)+dNdxi(2,2)*coord_l(1,2)+
     & dNdxi(3,2)*coord_l(1,3)+dNdxi(4,2)*coord_l(1,4)
      dxdxi(2,2) = dNdxi(1,2)*coord_l(2,1)+dNdxi(2,2)*coord_l(2,2)+
     & dNdxi(3,2)*coord_l(2,3)+dNdxi(4,2)*coord_l(2,4)       
c
c inverse dxdxi========================================================     
c
      detdxdxi=dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
      dxidx(1,1)=1.d0/detdxdxi*dxdxi(2,2)
      dxidx(2,2)=1.d0/detdxdxi*dxdxi(1,1)
      dxidx(1,2)=-1.d0/detdxdxi*dxdxi(1,2)
      dxidx(2,1)=-1.d0/detdxdxi*dxdxi(2,1)
c
c shape function derivative to local coordinate========================
c
      do ii = 1,4
         do j = 1,2
             do k = 1,2
                 dNdx(ii,j) = dNdx(ii,j) + dNdxi(ii,k) * dxidx(k,j)
             end do
         end do
      end do
      dNdx(1,1) = dNdxi(1,1)*dxidx(1,1)+dNdxi(1,2)*dxidx(2,1)
      dNdx(1,2) = dNdxi(1,1)*dxidx(1,2)+dNdxi(1,2)*dxidx(2,2)
      dNdx(2,1) = dNdxi(2,1)*dxidx(1,1)+dNdxi(2,2)*dxidx(2,1)
      dNdx(2,2) = dNdxi(2,1)*dxidx(1,2)+dNdxi(2,2)*dxidx(2,2)
      dNdx(3,1) = dNdxi(3,1)*dxidx(1,1)+dNdxi(3,2)*dxidx(2,1)
      dNdx(3,2) = dNdxi(3,1)*dxidx(1,2)+dNdxi(3,2)*dxidx(2,2)
      dNdx(4,1) = dNdxi(4,1)*dxidx(1,1)+dNdxi(4,2)*dxidx(2,1)
      dNdx(4,2) = dNdxi(4,1)*dxidx(1,2)+dNdxi(4,2)*dxidx(2,2)
c
c get strains
c
      do ii = 1,2
         do j = 1,4
             uv(ii,j) = disp_l(ii,j)
         end do
      end do
c
      call k_matrix_multiply(uv,dNdx,uv_xy,2,4,2)
c
      uv_xy(1,1) = uv(1,1)*dNdx(1,1)+uv(1,2)*dNdx(2,1)+uv(1,3)*dNdx(3,1)
     & +uv(1,4)*dNdx(4,1)
      uv_xy(1,2) = uv(1,1)*dNdx(1,2)+uv(1,2)*dNdx(2,2)+uv(1,3)*dNdx(3,2)
     & +uv(1,4)*dNdx(4,2)
      uv_xy(2,1) = uv(2,1)*dNdx(1,1)+uv(2,2)*dNdx(2,1)+uv(2,3)*dNdx(3,1)
     & +uv(2,4)*dNdx(4,1)
      uv_xy(2,2) = uv(2,1)*dNdx(1,2)+uv(2,2)*dNdx(2,2)+uv(2,3)*dNdx(3,2)
     & +uv(2,4)*dNdx(4,2)
c      
      eps_11 = uv_xy(1,1)*cos(up_sigma)+uv_xy(1,2)*sin(up_sigma)
      eps_12 = -uv_xy(1,1)*sin(up_sigma)+uv_xy(1,2)*cos(up_sigma)
      eps_21 = uv_xy(2,1)*cos(up_sigma)+uv_xy(2,2)*sin(up_sigma)
      eps_22 = -uv_xy(2,1)*sin(up_sigma)+uv_xy(2,2)*cos(up_sigma)
      
      up_eps_vm = abs(eps_22)
      if (i.EQ.1) then
         WRITE(20,*) uv_xy(1,1),',', uv_xy(2,1),',', uv_xy(1,2)
         WRITE(20,*) uv_xy(2,2)         
      end if      
c===========do it again to get bottom equivalent strain=======================
      dxdxi(1,1) = dNdxi(1,1)*coord_l(1,5)+dNdxi(2,1)*coord_l(1,6)+
     & dNdxi(3,1)*coord_l(1,7)+dNdxi(4,1)*coord_l(1,8)
      dxdxi(1,2) = dNdxi(1,1)*coord_l(2,5)+dNdxi(2,1)*coord_l(2,6)+
     & dNdxi(3,1)*coord_l(2,7)+dNdxi(4,1)*coord_l(2,8)
      dxdxi(2,1) = dNdxi(1,2)*coord_l(1,5)+dNdxi(2,2)*coord_l(1,6)+
     & dNdxi(3,2)*coord_l(1,7)+dNdxi(4,2)*coord_l(1,8)
      dxdxi(2,2) = dNdxi(1,2)*coord_l(2,5)+dNdxi(2,2)*coord_l(2,6)+
     & dNdxi(3,2)*coord_l(2,7)+dNdxi(4,2)*coord_l(2,8)

      detdxdxi=dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
      dxidx(1,1)=1.d0/detdxdxi*dxdxi(2,2)
      dxidx(2,2)=1.d0/detdxdxi*dxdxi(1,1)
      dxidx(1,2)=-1.d0/detdxdxi*dxdxi(1,2)
      dxidx(2,1)=-1.d0/detdxdxi*dxdxi(2,1)   

      dNdx(1,1) = dNdxi(1,1)*dxidx(1,1)+dNdxi(1,2)*dxidx(2,1)
      dNdx(1,2) = dNdxi(1,1)*dxidx(1,2)+dNdxi(1,2)*dxidx(2,2)
      dNdx(2,1) = dNdxi(2,1)*dxidx(1,1)+dNdxi(2,2)*dxidx(2,1)
      dNdx(2,2) = dNdxi(2,1)*dxidx(1,2)+dNdxi(2,2)*dxidx(2,2)
      dNdx(3,1) = dNdxi(3,1)*dxidx(1,1)+dNdxi(3,2)*dxidx(2,1)
      dNdx(3,2) = dNdxi(3,1)*dxidx(1,2)+dNdxi(3,2)*dxidx(2,2)
      dNdx(4,1) = dNdxi(4,1)*dxidx(1,1)+dNdxi(4,2)*dxidx(2,1)
      dNdx(4,2) = dNdxi(4,1)*dxidx(1,2)+dNdxi(4,2)*dxidx(2,2) 

      do ii = 1,2
         do j = 1,4
             uv(ii,j) = disp_l(ii,j)
         end do
      end do
     
      uv_xy(1,1) = uv(1,1)*dNdx(1,1)+uv(1,2)*dNdx(2,1)+uv(1,3)*dNdx(3,1)
     & +uv(1,4)*dNdx(4,1)
      uv_xy(1,2) = uv(1,1)*dNdx(1,2)+uv(1,2)*dNdx(2,2)+uv(1,3)*dNdx(3,2)
     & +uv(1,4)*dNdx(4,2)
      uv_xy(2,1) = uv(2,1)*dNdx(1,1)+uv(2,2)*dNdx(2,1)+uv(2,3)*dNdx(3,1)
     & +uv(2,4)*dNdx(4,1)
      uv_xy(2,2) = uv(2,1)*dNdx(1,2)+uv(2,2)*dNdx(2,2)+uv(2,3)*dNdx(3,2)
     & +uv(2,4)*dNdx(4,2) 

      eps_11 = uv_xy(1,1)*cos(up_sigma)+uv_xy(1,2)*sin(up_sigma)
      eps_12 = -uv_xy(1,1)*sin(up_sigma)+uv_xy(1,2)*cos(up_sigma)
      eps_21 = uv_xy(2,1)*cos(up_sigma)+uv_xy(2,2)*sin(up_sigma)
      eps_22 = -uv_xy(2,1)*sin(up_sigma)+uv_xy(2,2)*cos(up_sigma)
      
      bt_eps_vm = abs(eps_22)     
c===========================================================================      
      
      if (i.EQ.1) then
         WRITE(20,*) eps_y,',', bt_eps_vm,',', up_eps_vm   
      end if
c
c====================viscosity after coupling================================
c=========================give the old values================================
         Tt1_old = Svars(16.0+6.0*(i-1.0)+1.0)
         Tt2_old = Svars(16.0+6.0*(i-1.0)+2.0)
         Tn_old = Svars(16.0+6.0*(i-1.0)+3.0)   
         Et1_old = Svars(16.0+6.0*(i-1.0)+4.0)
         Et2_old = Svars(16.0+6.0*(i-1.0)+5.0)
         En_old = Svars(16.0+6.0*(i-1.0)+6.0)
c call cohesive law
      call k_cohesive_law(Trac,Trac_Jacob,G_fn,G_ft,f_tn,f_tt,K_n,K_t,
     & pmax,tmax,tmax1,tmax2,delu_loc_gp,sep_n0,sep_nf,sep_t0,sep_tf,
     & mcrd, nrhs, SVARS,Tt1_old,Tt2_old,Tn_old,Et1_old,Et2_old,En_old,
     & DTIME,up_eps_vm,bt_eps_vm)
!
      if (gpt.eq.1) then
      xstress = Trac(3,1)
      end if
     
     
c=============================store the new one==========================================
         Svars(16.0+6.0*(i-1.0)+1.0) = Trac(1,1)
         Svars(16.0+6.0*(i-1.0)+2.0) = Trac(2,1)
         Svars(16.0+6.0*(i-1.0)+3.0) = Trac(3,1)
         Svars(16.0+6.0*(i-1.0)+4.0) = Trac_Jacob(1,1)
         Svars(16.0+6.0*(i-1.0)+5.0) = Trac_Jacob(2,2)
         Svars(16.0+6.0*(i-1.0)+6.0) = Trac_Jacob(3,3)
      if ((Trac(1,1).eq.0.d0).and.(Trac(2,1).eq.0.d0).and.
     & (Trac(3,1).eq.0.d0).and.(Trac_Jacob(1,1).eq.0.d0).and.
     & (Trac_Jacob(2,2).eq.0.d0).and.(Trac_Jacob(3,3).eq.0.d0)) then 
          write(18,*) kinc,',',JELEM    
      end if         
          if (i.eq.1) then
              write(16,*) kinc,',', delu_loc_gp(3),',',trac(3,1)      
          end if
!      if ((KINC.EQ.1).AND.(i.EQ.1)) then
!         WRITE(16,*) delu_loc_gp(1), delu_loc_gp(2), delu_loc_gp(3)
!         WRITE(16,*) opt1, opt2, opn
!         WRITE(16,*) ds1(1), ds2(1), dn(1)
!         WRITE(16,*) ds1(2), ds2(2), dn(2)
!         WRITE(16,*) ds1(3), ds2(3), dn(3)
!         WRITE(16,*) ds1(4), ds2(4), dn(4)         
!         WRITE(16,*) disp_l(1,1), disp_l(2,1), disp_l(3,1)
!         WRITE(16,*) disp_l(1,5), disp_l(2,5), disp_l(3,5)
!         WRITE(16,*) disp(1,1), disp(2,1), disp(3,1)
!         WRITE(16,*) disp(1,5), disp(2,5), disp(3,5)   
!         WRITE(16,*) disp(1,2), disp(2,2), disp(3,2)
!         WRITE(16,*) disp(1,6), disp(2,6), disp(3,6) 
!         WRITE(16,*) disp(1,4), disp(2,4), disp(3,4)
!         WRITE(16,*) disp(1,8), disp(2,8), disp(3,8)          
!         WRITE(16,*) Trac(1,1), Trac(2,1), Trac(3,1)
!      end if
c
c Determine B matrix and its transpose
c 
         call k_Bmatrix(sf,B,mcrd,ndofel)
c
         call k_matrix_transpose(B,B_t,mcrd,ndofel)
c
c Compute the stiffness matrix
c Local Stiffness = B_t * Trac_Jacob * B
c
         call k_matrix_multiply(Trac_Jacob,B,temp1,mcrd,mcrd,
     & ndofel)
         call k_matrix_multiply(B_t,temp1,stiff_l,ndofel,
     & mcrd,ndofel)
c
c Compute Global stiffness matrix 
c Global_K = T' * K * T
c
         call k_matrix_multiply(Transformation_M_T,stiff_l,
     & temp2,ndofel,ndofel,ndofel)
         call k_matrix_multiply(temp2,Transformation_M,stiff_g,
     & ndofel,ndofel,ndofel)
c
c Multiply Jacobian with the Global stiffness and add contribution
c from each Gauss Point
c
         call k_matrix_plus_scalar(amatrx,stiff_g,a_Jacob,
     & ndofel,ndofel)
c
c Compute the global residual vector
c Local_residual = B_t * Trac
c Global_residual = T' * Local_residual
c 
         call k_matrix_multiply(B_t,Trac,residual_l,ndofel,
     & mcrd,nrhs)

         call k_matrix_multiply(Transformation_M_T,residual_l,
     & residual_g,ndofel,ndofel,nrhs)

c
c Multiply the Global residual by the Jacobian and add the 
c contribution from each point
c   
         call k_matrix_plus_scalar(rhs,residual_g,a_Jacob,
     & ndofel,nrhs)
!        if ((KINC.EQ.1) .AND. (gpt.EQ.1)) then 
!            WRITE(16,*) delu_loc_gp(3), Trac(3,1), r(1,1)
!             WRITE(16,*) r(1,1), r(1,2), r(1,3)
!             WRITE(16,*) r(2,1), r(2,2), r(2,3)
!             WRITE(16,*) r(3,1), r(3,2), r(3,3)
!             WRITE(16,*) rhs(8,1),rhs(7,1),rhs(6,1)
!             WRITE(16,*) rhs(5,1),rhs(4,1),rhs(3,1)
!             WRITE(16,*) rhs(2,1),rhs(1,1)
!             WRITE(16,*) amatrx(1,1), amatrx(1,2), amatrx(1,3)
!             WRITE(16,*) amatrx(2,1), amatrx(2,2), amatrx(2,3)
!             WRITE(16,*) amatrx(3,1), amatrx(3,2), amatrx(3,3)
!             WRITE(16,*) amatrx(4,4), amatrx(4,5), amatrx(4,6)
!             WRITE(16,*) amatrx(5,4), amatrx(5,5), amatrx(5,6)
!             WRITE(16,*) amatrx(6,4), amatrx(6,5), amatrx(6,6)
!             WRITE(16,*) amatrx(7,7), amatrx(7,8), amatrx(7,9)
!             WRITE(16,*) amatrx(8,7), amatrx(8,8), amatrx(8,9)
!             WRITE(16,*) amatrx(9,7), amatrx(9,8), amatrx(9,9)             
!         end if
      end do
c
      return
      end
c======================================================================
c=============================SUBROUTINES==============================
c======================================================================
c
c Determine the strain-displacement (B) matrix
c
      subroutine k_Bmatrix(sf,B,mcrd,ndofel)
      INCLUDE 'ABA_PARAM.INC'
      dimension sf(4),B(mcrd,ndofel)
      B(1,1) =  sf(1)
      B(1,4) =  sf(2)
      B(1,7) =  sf(3)
      B(1,10)=  sf(4)
      B(1,13)= -sf(1)
      B(1,16)= -sf(2)
      B(1,19)= -sf(3)
      B(1,22)= -sf(4)
      B(2,2) =  sf(1)
      B(2,5) =  sf(2)
      B(2,8) =  sf(3)
      B(2,11)=  sf(4)
      B(2,14)= -sf(1)
      B(2,17)= -sf(2)
      B(2,20)= -sf(3)
      B(2,23)= -sf(4)
      B(3,3) =  sf(1)
      B(3,6) =  sf(2)
      B(3,9) =  sf(3)
      B(3,12)=  sf(4)
      B(3,15)= -sf(1)
      B(3,18)= -sf(2)
      B(3,21)= -sf(3)
      B(3,24)= -sf(4)
c
      return
      end
c======================================================================
      subroutine k_cohesive_law(T,T_d,G_fn,G_ft,f_tn,f_tt,K_n,K_t,
     & pmax,tmax,tmax1,tmax2,delu,sep_n0,sep_nf,sep_t0,sep_tf,
     & mcrd, nrhs, SVARS,Tt1_old,Tt2_old,Tn_old,Et1_old,Et2_old,En_old,
     & DTIME,up_eps_vm,bt_eps_vm)
     
      INCLUDE 'ABA_PARAM.INC'
      dimension T(mcrd,nrhs),T_d(mcrd,mcrd),delu(mcrd),SVARS(40)
       DOUBLE PRECISION G_fn, G_ft, f_tn, f_tt, K_n, K_t,
     & tmax1, tmax2, pmax, tmax, popn, popt1, popt2, D_n, D_t1, D_t2,
     & Tn, Tt1, Tt2, Kn, Kt1, Kt2, T_d, T,Tt1_old,Tt2_old,Tn_old,
     & Et1_old,Et2_old,En_old,eta,
     & delu, sep_n0, sep_nf, sep_t0, sep_tf, up_eps_vm, bt_eps_vm
     
       PARAMETER (DMGMAX = 1.0D0, ONE = 1.d0, ZERO = 0.d0)
c     
      popt1=delu(1)
      popt2=delu(2)
      popn=delu(3)
      eps_max = 0.02d0
      eta = 0.00000d0
c
c 
c======================check coupling condition=======================================
c
      if ((up_eps_vm .LE. eps_max).or.(bt_eps_vm .LE. eps_max)) then
     
         if (tmax1 .LE. sep_t0) then
             Tt1=K_t*popt1
             Kt1=K_t
         else if (tmax1 .GT. sep_t0) Then
             D_t1=sep_tf*(tmax1-sep_t0)/tmax1/(sep_tf-sep_t0)
             !to prevent damage larger than 1
             if (D_t1.ge.DMGMAX) then
                D_t1 = DMGMAX
                Tt1 = ZERO
                Tt2 = ZERO
                Tn = ZERO
                Kt1=ZERO
                Kt2=ZERO
                Kn=ZERO
             end if
             if ((popt1.lt.tmax1).AND.(D_t1.lt.DMGMAX)) then
                Tt1= (ONE-D_t1)*K_t*popt1
                Kt1 = (ONE-D_t1)*K_t
             else if ((popt1.eq.tmax1).AND.(D_t1.lt.DMGMAX)) then
                Tt1=(ONE-D_t1)*K_t*popt1
                Kt1=(D_t1-sep_tf/(sep_tf-sep_t0))*K_t+(ONE-D_t1)*K_t   
             end if 
         end if
c
         if (tmax2 .LE. sep_t0) then
             Tt2=K_t*popt2
             Kt2=K_t
         else if (tmax2 .GT. sep_t0) Then
             D_t2=sep_tf*(tmax2-sep_t0)/tmax2/(sep_tf-sep_t0)
             !to prevent damage larger than 1
             if (D_t2.ge.DMGMAX) then
                D_t2 = DMGMAX
                Tt2 = ZERO
                Tt1 = ZERO
                Tn = ZERO
                Kt2=ZERO
                Kt1=ZERO
                Kn=ZERO
             end if
             if ((popt2.lt.tmax2).AND.(D_t2.lt.DMGMAX)) then
                Tt2= (ONE-D_t2)*K_t*popt2
                Kt2 = (ONE-D_t2)*K_t
             else if ((popt2.eq.tmax2).AND.(D_t2.lt.DMGMAX)) then
                Tt2=(ONE-D_t2)*K_t*popt2
                Kt2=(D_t2-sep_tf/(sep_tf-sep_t0))*K_t+(ONE-D_t2)*K_t   
             end if 
         end if
c
        if (popn .LE. ZERO) then
             Kn=K_n
             Tn=Kn*popn
        end if
        if (popn .GE. ZERO) then
         if (pmax .LE. sep_n0) then
             Tn=K_n*popn
             Kn=K_n
         else if (pmax .GT. sep_n0) Then
             D_n=sep_nf*(pmax-sep_n0)/pmax/(sep_nf-sep_n0)
             !to prevent damage larger than 1
             if (D_n.ge.DMGMAX) then
                D_n = DMGMAX
                Tn = ZERO
                Tt1 = ZERO
                Tt2 = ZERO
                Kn=ZERO
                Kt1=ZERO
                Kt2=ZERO
             end if
             if ((popn.lt.pmax).AND.(D_n.lt.DMGMAX)) then
                Tn= (ONE-D_n)*K_n*popn
                Kn = (ONE-D_n)*K_n
             else if ((popn.eq.pmax).AND.(D_n.lt.DMGMAX)) then
                Tn=(ONE-D_n)*K_n*popn
                Kn=(D_n-sep_nf/(sep_nf-sep_n0))*K_n+(ONE-D_n)*K_n   
             end if 
         end if
        end if
      else
         Tt1 = 0.d0
         Tt2 = 0.d0
         Tn = 0.d0
         Kt1 = 0.d0
         Kt2 = 0.d0
         Kn = 0.d0
c===========================calculate viscosity values=========================
         Tt1 = eta / (eta + DTIME) * Tt1_old + DTIME / (eta + DTIME)*Tt1
         Tt2 = eta / (eta + DTIME) * Tt2_old + DTIME / (eta + DTIME)*Tt2
         Tn = eta / (eta + DTIME) * Tn_old + DTIME / (eta + DTIME) * Tn
         Kt1 = eta / (eta + DTIME) * Et1_old + DTIME / (eta + DTIME)*Kt1
         Kt2 = eta / (eta + DTIME) * Et2_old + DTIME / (eta + DTIME)*Kt2
         Kn = eta / (eta + DTIME) * En_old + DTIME / (eta + DTIME) * Kn
      end if

c
c
      T(1,1) = Tt1
      T(2,1) = Tt2
      T(3,1) = Tn
c
      T_d(1,1)=Kt1
      T_d(1,2)=0.0d00
      T_d(1,3)=0.0d00
      T_d(2,1)=0.0d00
      T_d(2,2)=Kt2
      T_d(2,3)=0.0d00
      T_d(3,1)=0.0d00
      T_d(3,2)=0.0d00
      T_d(3,3)=Kn
c
      return
      end
c=====================================================================
      subroutine k_local_coordinates(gpt,co_de,R,coord_l,
     1 Transformation_M,
     & Transformation_M_T,a_Jacob,aJacob_M,coords,u,ndofel,nnode,
     & mcrd, SVARS,disp,disp_l)
      INCLUDE 'ABA_PARAM.INC'
      dimension R(mcrd,mcrd),coord_l(mcrd,nnode),aJacob_M(2,3),
     & Transformation_M(ndofel,ndofel),coords(mcrd,nnode),
     & Transformation_M_T(ndofel,ndofel),u(ndofel),
     & co_de(mcrd,nnode), co_de_m(3,4),SFD(2,4),SVARS(40),
     & disp(mcrd,nnode), disp_l(mcrd,nnode)
       DOUBLE PRECISION x1, x2, x3, x4, y1, y2, y3, y4, y5, y6, z1, z2,
     & z3, z4
c
      call k_matrix_zero(co_de_m,3,4)
c
      do i = 1, 3
         co_de_m(i,1)=(co_de(i,1)+co_de(i,5))*0.5d00
         co_de_m(i,2)=(co_de(i,2)+co_de(i,6))*0.5d00
         co_de_m(i,3)=(co_de(i,3)+co_de(i,7))*0.5d00
         co_de_m(i,4)=(co_de(i,4)+co_de(i,8))*0.5d00
      end do
c
      x1=co_de_m(1,1)
      x2=co_de_m(1,2)
      x3=co_de_m(1,3)
      x4=co_de_m(1,4)
c
      y1=co_de_m(2,1)
      y2=co_de_m(2,2)
      y3=co_de_m(2,3)
      y4=co_de_m(2,4)
c
      z1=co_de_m(3,1)
      z2=co_de_m(3,2)
      z3=co_de_m(3,3)
      z4=co_de_m(3,4)
c
      if (gpt .eq. 1) then
         c_r=-1.0d0
         c_s=-1.0d0
      elseif (gpt .eq. 2) then
         c_r= 1.0d0
         c_s=-1.0d0
      elseif (gpt .eq. 3) then
         c_r= 1.0d0
         c_s= 1.0d0
      elseif (gpt .eq. 4) then
         c_r=-1.0d0
         c_s= 1.0d0
      end if
c
      SFD(1,1) =-0.25d0*(1-c_s)
      SFD(1,2) = 0.25d0*(1-c_s)
      SFD(1,3) = 0.25d0*(1+c_s)
      SFD(1,4) =-0.25d0*(1+c_s)
      SFD(2,1) =-0.25d0*(1-c_r)
      SFD(2,2) =-0.25d0*(1+c_r)
      SFD(2,3) = 0.25d0*(1+c_r)
      SFD(2,4) = 0.25d0*(1-c_r)
c
      do i = 1,2
         do j = 1,3
            do k =1, 4
               aJacob_M(i,j) = aJacob_M(i,j) + SFD(i,k)*co_de_m(j,k)
            end do
         end do
      end do
c
      dum1 = aJacob_M(1,2)*aJacob_M(2,3) - aJacob_M(1,3)*aJacob_M(2,2)
      dum2 = aJacob_M(1,3)*aJacob_M(2,1) - aJacob_M(1,1)*aJacob_M(2,3)
      dum3 = aJacob_M(1,1)*aJacob_M(2,2) - aJacob_M(1,2)*aJacob_M(2,1)
c
      a_Jacob = sqrt(dum1**2 + dum2**2 + dum3**2)
c
      R(3,1) = dum1/a_Jacob
      R(3,2) = dum2/a_Jacob
      R(3,3) = dum3/a_Jacob
c
      aLen=sqrt(aJacob_M(1,1)**2.0d00+aJacob_M(1,2)**2.0d00
     1 +aJacob_M(1,3)**2.0d00)
      R(1,1)=aJacob_M(1,1)/aLen
      R(1,2)=aJacob_M(1,2)/aLen
      R(1,3)=aJacob_M(1,3)/aLen
c
      aLen2a=R(3,2)*R(1,3)-R(3,3)*R(1,2)
      aLen2b=R(3,3)*R(1,1)-R(3,1)*R(1,3)
      aLen2c=R(3,1)*R(1,2)-R(3,2)*R(1,1)
c
      aLen2 = sqrt(aLen2a**2.0d00 + aLen2b**2.0d00 + aLen2c**2.0d00)
c
      R(2,1) = aLen2a/aLen2
      R(2,2) = aLen2b/aLen2
      R(2,3) = aLen2c/aLen2
c
      a_J11 = (-0.25d0*(1.0d0-c_s))*x1 + (0.25d0*(1.0d0-c_s))*x2 + 
     & (0.25d0*(1.0d0+c_s))*x3 + (-0.25d0*(1.0d0+c_s))*x4
      a_J12 = (-0.25d0*(1.0d0-c_s))*z1 + (0.25d0*(1.0d0-c_s))*z2 + 
     & (0.25d0*(1.0d0+c_s))*z3 + (-0.25d0*(1.0d0+c_s))*z4
      a_J21 = (-0.25d0*(1.0d0-c_r))*x1 + (-0.25d0*(1.0d0+c_r))*x2 + 
     & (0.25d0*(1.0d0+c_r))*x3 + (0.25d0*(1.0d0-c_r))*x4
      a_J22 = (-0.25d0*(1.0d0-c_r))*z1 + (-0.25d0*(1.0d0+c_r))*z2 + 
     & (0.25d0*(1.0d0+c_r))*z3 + (0.25d0*(1.0d0-c_r))*z4
c
      b_J11 = (-0.25d0*(1.0d0-c_s))*x1 + (0.25d0*(1.0d0-c_s))*x2 + 
     & (0.25d0*(1.0d0+c_s))*x3 + (-0.25d0*(1.0d0+c_s))*x4
      b_J12 = (-0.25d0*(1.0d0-c_s))*y1 + (0.25d0*(1.0d0-c_s))*y2 + 
     & (0.25d0*(1.0d0+c_s))*y3 + (-0.25d0*(1.0d0+c_s))*y4
      b_J21 = (-0.25d0*(1.0d0-c_r))*x1 + (-0.25d0*(1.0d0+c_r))*x2 + 
     & (0.25d0*(1.0d0+c_r))*x3 + (0.25d0*(1.0d0-c_r))*x4
      b_J22 = (-0.25d0*(1.0d0-c_r))*y1 + (-0.25d0*(1.0d0+c_r))*y2 + 
     & (0.25d0*(1.0d0+c_r))*y3 + (0.25d0*(1.0d0-c_r))*y4
c
      c_J11 = (-0.25d0*(1.0d0-c_s))*y1 + (0.25d0*(1.0d0-c_s))*y2 + 
     & (0.25d0*(1.0d0+c_s))*y3 + (-0.25d0*(1.0d0+c_s))*y4
      c_J12 = (-0.25d0*(1.0d0-c_s))*z1 + (0.25d0*(1.0d0-c_s))*z2 + 
     & (0.25d0*(1.0d0+c_s))*z3 + (-0.25d0*(1.0d0+c_s))*z4
      c_J21 = (-0.25d0*(1.0d0-c_r))*y1 + (-0.25d0*(1.0d0+c_r))*y2 + 
     & (0.25d0*(1.0d0+c_r))*y3 + (0.25d0*(1.0d0-c_r))*y4
      c_J22 = (-0.25d0*(1.0d0-c_r))*z1 + (-0.25d0*(1.0d0+c_r))*z2 + 
     & (0.25d0*(1.0d0+c_r))*z3 + (0.25d0*(1.0d0-c_r))*z4
c
c
      a_Jacob1 = (a_J11*a_J22 - a_J12*a_J21)
      a_Jacob2 = (b_J11*b_J22 - b_J12*b_J21)
      a_Jacob3 = (c_J11*c_J22 - c_J12*c_J21)
c
      a_Jacob = sqrt(a_Jacob1**2.0d00+a_Jacob2**2.0d00+a_Jacob3**2.0d00)
c
c======================================================================
      num=nnode
c
      do i = 1, num
         dum=3.0*(i-1.0)
         Transformation_M(dum+1,dum+1)=R(1,1)
         Transformation_M(dum+1,dum+2)=R(1,2) 
         Transformation_M(dum+1,dum+3)=R(1,3)
         Transformation_M(dum+2,dum+1)=R(2,1)
         Transformation_M(dum+2,dum+2)=R(2,2)
         Transformation_M(dum+2,dum+3)=R(2,3)
         Transformation_M(dum+3,dum+1)=R(3,1)
         Transformation_M(dum+3,dum+2)=R(3,2)
         Transformation_M(dum+3,dum+3)=R(3,3)
      end do
c
      call k_matrix_transpose(Transformation_M,Transformation_M_T,
     1 ndofel,ndofel)
c
c
      do i = 1, nnode
         coord_l(1,i)=(R(1,1)*co_de(1,i)+R(1,2)*co_de(2,i)
     & +R(1,3)*co_de(3,i))
         coord_l(2,i)=(R(2,1)*co_de(1,i)+R(2,2)*co_de(2,i)
     & +R(2,3)*co_de(3,i))
         coord_l(3,i)=(R(3,1)*co_de(1,i)+R(3,2)*co_de(2,i)
     & +R(3,3)*co_de(3,i))
      end do
c
      do i = 1, nnode
         disp_l(1,i)=(R(1,1)*disp(1,i)+R(1,2)*disp(2,i)
     & +R(1,3)*disp(3,i))
         disp_l(2,i)=(R(2,1)*disp(1,i)+R(2,2)*disp(2,i)
     & +R(2,3)*disp(3,i))
         disp_l(3,i)=(R(3,1)*disp(1,i)+R(3,2)*disp(2,i)
     & +R(3,3)*disp(3,i))
      end do      
c    
      return
      end
c======================================================================
      subroutine k_shape_fun(i,sf)
      INCLUDE 'ABA_PARAM.INC'
      dimension sf(4), GP_coord(2)
c
      if (i .eq. 1) then
         GP_coord(1)=-1.0d0
         GP_coord(2)=-1.0d0
      elseif (i .eq. 2) then
         GP_coord(1)= 1.0d0
         GP_coord(2)=-1.0d0
      elseif (i .eq. 3) then
         GP_coord(1)= 1.0d0
         GP_coord(2)= 1.0d0
      elseif (i .eq. 4) then
         GP_coord(1)=-1.0d0
         GP_coord(2)= 1.0d0
      end if
c
      sf(1)=(1-GP_coord(1))*(1-GP_coord(2))*0.25d0
      sf(2)=(1+GP_coord(1))*(1-GP_coord(2))*0.25d0
      sf(3)=(1+GP_coord(1))*(1+GP_coord(2))*0.25d0
      sf(4)=(1-GP_coord(1))*(1+GP_coord(2))*0.25d0
c
      return
      end
c======================================================================
      subroutine k_matrix_multiply(A,B,C,l,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(l,n),B(n,m),C(l,m)
c
      call k_matrix_zero(C,l,m)
c
      do i = 1, l
         do j = 1, m
            do k = 1, n
               C(i,j)=C(i,j)+A(i,k)*B(k,j)
            end do
         end do
      end do
c
      return
      end
c======================================================================
      subroutine k_matrix_plus_scalar(A,B,c,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m),B(n,m)
c
      do i = 1, n
         do j = 1, m
            A(i,j)=A(i,j)+c*B(i,j)
         end do
      end do
c
      return
      end
c======================================================================
      subroutine k_matrix_transpose(A,B,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m),B(m,n)
c
      do i = 1, n
         do j = 1, m
            B(j,i)=A(i,j)
         end do
      end do
c
      return
      end
c======================================================================
      subroutine k_matrix_zero(A,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m)
c
      do i = 1, n
         do j = 1, m
            A(i,j)=0.d0
         end do
      end do
c
      return
      end
c======================================================================
      subroutine k_vector_zero(A,n)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n)
c
      do i = 1, n
         A(i)=0.d0
      end do
c
      return
      end
c======================================================================
      subroutine k_Mac(pM,a,b)
      INCLUDE 'ABA_PARAM.INC'
c
      if ((a-b) .GE. 0.0) then
         pM=a-b
      elseif ((a-b) .LT. 0.0) then
         pM=0.d0
      end if
c
      return
      end
c======================================================================
c=================================END==================================
c======================================================================









      
