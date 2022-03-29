 Program scf

      Implicit None
      Integer :: natm
      Integer :: nelctrn
      Integer :: nocc
      Integer :: unocc
      Integer :: nbasis
      Integer :: nbasfn
      Integer :: ldim
      Integer :: i,j,l,k,m,n
      Integer :: ia,jb
      Integer :: itr

      Integer :: p,q,r
      Integer :: a,b
      Integer :: tmp_1
      Integer :: mu,nu,lam,sig
      !Real :: sum_val
      Double precision :: tmp_2
      Double precision :: val_1, val_2
      Double precision :: E_nuc,sumn, tmp
      Double precision :: sum_E, E_old, E_new !HF_SCF Energy
      Double precision :: diff_E
      Double precision :: E_scf,E_mp2, sum_mp2
      Integer :: compare

      Double precision, DIMENSION(:), ALLOCATABLE :: oneh_tmp
      Double precision, DIMENSION(600) :: buf
      Double precision, DIMENSION(600) :: ibuf
      Double precision, DIMENSION(600) :: buf2
      Double precision, DIMENSION(600) :: ibuf2
      Double precision, DIMENSION(:,:), ALLOCATABLE :: oneh

      Double precision, DIMENSION(:), ALLOCATABLE :: ovrlp_tmp
      Double precision, DIMENSION(:,:), ALLOCATABLE :: ovrlp

      !ovrlp_mat is for check c'sc = 1
      Double precision, DIMENSION(:,:), ALLOCATABLE :: ovrlp_mat

      !Eigen vector matrix of overlap mat  
      Double precision, DIMENSION(:,:), ALLOCATABLE :: ovrlp_evec

      !Ovrlp's Eigen value matrix's square root matrix == s
      Double precision, DIMENSION(:,:), ALLOCATABLE :: s
      Double precision, DIMENSION(:,:), ALLOCATABLE :: x_tmp
      Double precision, DIMENSION(:,:), ALLOCATABLE :: x_untry

      Double precision, DIMENSION(:,:), ALLOCATABLE :: G
      !Double precision, DIMENSION(:,:), ALLOCATABLE :: sumn

      Double precision, DIMENSION(:,:), ALLOCATABLE :: coeff_mat
      Double precision, DIMENSION(:,:), ALLOCATABLE :: C_dash
      Double precision, DIMENSION(:,:), ALLOCATABLE :: C

      Double precision, DIMENSION(:,:), ALLOCATABLE :: cdash_sc
      Double precision, DIMENSION(:,:), ALLOCATABLE :: tmp_sc

      Double precision, DIMENSION(:,:), ALLOCATABLE :: dens_mat
      !P is the density  matrix
      Double precision, DIMENSION(:,:), ALLOCATABLE :: old_P
      Double precision, DIMENSION(:,:), ALLOCATABLE :: new_P
      Double precision, DIMENSION(:,:), ALLOCATABLE :: P_mat
      Double precision, DIMENSION(:,:), ALLOCATABLE :: nw_G
      Double precision, DIMENSION(:,:), ALLOCATABLE :: nwfock_mat
      Double precision, DIMENSION(:,:), ALLOCATABLE :: nwfock_tmp
      Double precision, DIMENSION(:,:), ALLOCATABLE :: nwfock_dash
      Double precision, DIMENSION(:), ALLOCATABLE :: orb_E



      Double precision, DIMENSION(:,:), ALLOCATABLE :: fock_mat
      Double precision, DIMENSION(:,:), ALLOCATABLE :: fock_tmp
      Double precision, DIMENSION(:,:), ALLOCATABLE :: fock_dash

      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: twoeint
      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: mo_2eint
      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: i_mo2eint
      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: j_mo2eint
      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: k_mo2eint
      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: l_mo2eint
      Double precision, DIMENSION(:,:,:,:), ALLOCATABLE :: spin_int
      Double precision, DIMENSION(:,:), ALLOCATABLE :: H_cis
      Double precision, DIMENSION(:,:), ALLOCATABLE :: fs


      Double precision, DIMENSION(:,:), ALLOCATABLE :: fs_vec
      Double precision, DIMENSION(:,:), ALLOCATABLE :: e_vec



      Call aces_init_rte
      Call aces_ja_init

C Your scf logic can start here. Needs to read from ACES II files
C One and two electron integrals. Number of basis functions, Number
C of electrons (charge of the molecule), nuclear repulsion energy


      call getrec(1, 'JOBARC', 'NATOMS', 1, natm)
      call getrec(1, 'JOBARC', 'NBASTOT', 1, nbasis)
      call getrec(1, 'JOBARC', 'NAOBASFN', 1, nbasfn)
      call getrec(1, 'JOBARC', 'NMPROTON', 1, nelctrn)
      call getrec(1, 'JOBARC', 'NUCREP', 1, E_nuc)
      print *,  "********************* HF-SCF CODE*********************"
      print *, ""
      print *, "#atoms= ", natm
      print *, "#basis= ", nbasis
      print *, "#basis func = ", nbasfn
      print *, "#electrons= #protons= ",nelctrn
      print *, "E_nuc=", E_nuc
      print *, ""

      ldim = nbasis*(nbasis+1)/2
      nocc = nelctrn/2
      unocc = nbasis - nocc

      print *, "#nocc= ",nocc
      print *, "#unocc=", unocc
      print *, ""


C     _____________________  ONE-ELECTRON INTEGRALS
      allocate ( oneh_tmp(ldim) )
      allocate ( ovrlp_tmp(ldim) )
      allocate ( oneh(nbasis, nbasis) )
      allocate ( ovrlp(nbasis, nbasis) )
      allocate ( ovrlp_mat(nbasis, nbasis) )

      call Get1EInt(oneh_tmp,ovrlp_tmp,buf,ibuf,ldim)
      call EXPND2(oneh_tmp,oneh,nbasis)
      call EXPND2(ovrlp_tmp,ovrlp,nbasis)
      ovrlp_mat = ovrlp

      print *, ""
      print *, "Read 1-e Integrals"
      print *, ""



C     ____________________ TWO-ELECTRON INTEGRALS 
      allocate ( twoeint(nbasfn, nbasfn, nbasfn, nbasfn) )
      twoeint = 0.0d0

      call Get2Ints(twoeint,buf2,ibuf2,nbasfn,ldim)

      print *, "Read 2-e Integrals"
      print *, ""

C     ______________Diagonalize the OVerlap MAtrix 
C     (eq: 3.166, Szabo &Ostlund) 

      allocate ( ovrlp_evec(nbasis, nbasis) )
      Call EIG(ovrlp, ovrlp_evec, i, nbasis, 1)
      print *, "Diagonalized the ovrlp matrix"
      print *, ""
      ! Now ovrlp is a matrix of eigenvalues
      ! Now ovrlp_evec is a matrix of eigenvectors

C     ____________________GET X-UNITARY MATRIX_SYMM ORTHOGONALIZATION
C     (eq: 3.167, Szabo &Ostlund) 

      allocate ( s(nbasis, nbasis) )
      allocate ( x_tmp(nbasis, nbasis) )
      allocate ( x_untry(nbasis, nbasis) )
      x_tmp = 0.0
      x_untry = 0.0


      s = 0.0

      do i=1,nbasis
         s(i,i) = 1.0/sqrt(ovrlp(i,i))
      end do


      call DGEMM("n","n",nbasis,nbasis,nbasis,1.0D0,ovrlp_evec,nbasis,
     &    s,nbasis,0.0D0,x_tmp,nbasis)
      call DGEMM("n","t",nbasis,nbasis,nbasis,1.0D0,x_tmp,nbasis,
     &   ovrlp_evec,nbasis,0.0D0,x_untry,nbasis)

      print *, "Done Symmetric Orthogonalization, got X-unitary matrix"
      print *, ""
      !print *, x_untry



C      ________________INITIAL GUESS DENSITY MATRIX__________________
      allocate ( dens_mat(nbasis, nbasis) )
      dens_mat = 0.0d0
      print *, "Obtained Initial Guess Density Matrix"
      print *, ""
      !print *, dens_mat
      print *, ""





C_________________________________________________________________
C_________________________________________________________________
C_________________________________________________________________
C_________________________________________________________________
C     ____________ITERATIONS___________      


      allocate (nw_G(nbasis, nbasis))
      allocate (nwfock_mat(nbasis, nbasis))
      allocate (nwfock_tmp(nbasis, nbasis))
     allocate (nwfock_dash(nbasis, nbasis))
      allocate (C_dash(nbasis, nbasis))
      allocate (C(nbasis, nbasis))
      allocate (new_P(nbasis, nbasis))
      allocate (P_mat(nbasis, nbasis))
      allocate (tmp_sc(nbasis, nbasis))
      allocate (cdash_sc(nbasis, nbasis))
      allocate (orb_E( nbasis))


      allocate (i_mo2eint(nbasis,nbasis,nbasis,nbasis))
      allocate (j_mo2eint(nbasis,nbasis,nbasis,nbasis))
      allocate (k_mo2eint(nbasis,nbasis,nbasis,nbasis))
      allocate (l_mo2eint(nbasis,nbasis,nbasis,nbasis))
      allocate (mo_2eint(nbasis, nbasis, nbasis,nbasis))

      new_P = dens_mat
      itr =  1
      E_old = 0.0
      diff_E = 1.2
C     ______________Moving new_P to old _P
      allocate (old_P(nbasis,nbasis))

C________LOOP____________________
      do while (diff_E .ge. 0.00000001)
      print *, ""
      print *,"#",itr,"TH ITERATION ####################"
      print *, ""
          old_P = 0.0d0
          old_P = new_P
          !print *, "In itr OLD MAT will be New MAt of 0th itr"
          !print *, old_P


      nwfock_mat = 0.0d0
      nwfock_tmp = 0.d0
      nwfock_dash= 0.0d0
      !Calculate new G from density mat old_P and twoeint
C     (eq: 3.154, Szabo &Ostlund) 

      do i = 1,nbasis
         do j = 1,nbasis
            nw_G(i,j) = 0.0d0
            do k = 1,nbasis
               do l = 1,nbasis
      nw_G(i,j) = nw_G(i,j)
     &  +old_P(k,l)*(twoeint(i,j,l,k)-(0.5*twoeint(i,k,l,j)))
               end do
            end do
          end do
      end do



      nwfock_mat = oneh + nw_G
      print *, "Obtained **New Fock** Matrix"
      !print *, nwfock_mat

      !_______________ENERGY CALCULATION__________
C     (eq: 3.184, Szabo &Ostlund) 

      print *, ""
      print *, "****Calculating E_SCF******** "
      print *, ""
      sum_E = 0.0d0
      E_new = 0.0d0
      do i = 1, nbasis
         do j = 1, nbasis
           sum_E = sum_E + old_P(i,j)*(oneh(i,j) + nwfock_mat(i,j))
         end do
      end do
      E_new = (0.5*sum_E) + E_nuc
      print *, ""
      print *, "E_new =",E_new
      print *, ""
      print *, "E_old=", E_old
      diff_E = abs(E_new - E_old)

      !________MOVING NEW_E TO OLD_E 
      E_old = E_new
      !print *, ""
      !print *, "E_old=", E_old

      !_________ TRANSFORMED FOCK MATRIX_______
C     (eq: 3.177, Szabo &Ostlund) 

      nwfock_tmp = 0.0d0

      call DGEMM('t','n',nbasis,nbasis,nbasis,1.0d0,x_untry,
     & nbasis,nwfock_mat,nbasis,0.0d0,nwfock_tmp,nbasis)

      print *, "Obtained New Intermediate Fock Matrix"
      print *, ""
      !print *, nwfock_tmp

      call DGEMM('n','n',nbasis,nbasis,nbasis,1.0d0,nwfock_tmp,
     &       nbasis,x_untry,nbasis,0.0d0,nwfock_dash,nbasis)
      print *, "Obtained New Transformed Fock Matrix"
      print *, ""
      !print *, nwfock_dash

      !___________DIAGONALIZE THE FOCK MATRIX_____
C     (eq: 3.178, Szabo &Ostlund) 

      C_dash = 0.0d0
      orb_E = 0.0d0
      Call EIG(nwfock_dash, C_dash, i, nbasis, 1)
      do i = 1, nbasis
         orb_E(i) = nwfock_dash(i,i)
      end do
      print *, "*****Orbital Energies"
      !print *, orb_E

        !Now fock_dash contains the eigenvalues i.e. \varepsilon mat
        !Now C_dash is a matrix of eigenvectors of transformed fock mat 
        !print *, "Diagonalized Transformed Fock Matrix"



C      ______________NEW COEFFICIENT MARIX_______
C     (step 9, page146, Szabo &Ostlund) 


      C = 0.0d0
      call DGEMM('n','n',nbasis,nbasis,nbasis,1.0d0,x_untry,
     &       nbasis,C_dash,nbasis,0.0d0,C,nbasis)
      print *, "Obtained New Coefficient Matrix"
      print *, ""
      !print *, C
      print *, ""

C     ______________ CHECK C'SC SHOULD BE 1
      tmp_sc = 0.0d0
      cdash_sc = 0.0d0
      call DGEMM('t','n',nbasis,nbasis,nbasis,1.0d0,C,
     &       nbasis,ovrlp_mat,nbasis,0.0d0,tmp_sc,nbasis)
      call DGEMM('n','n',nbasis,nbasis,nbasis,1.0d0,tmp_sc,
     &       nbasis,C,nbasis,0.0d0,cdash_sc,nbasis)
      !print *, ""
      !print *, "CHECK C'SC = 1"
      !print *, ""
      !print *, cdash_sc
      !print *, ""

C     ______________OBTAIN NEW DENSITY MATRIX__
C     (eq: 3.145, Szabo &Ostlund) 

      new_P = 0.0d0
      P_mat = 0.0d0

      do i=1, nbasis
        do j=1, nbasis
           do k=1, nelctrn/2
           P_mat(i,j) = P_mat(i,j) + C(i,k)*C(j,k)
           end do
        end do


      new_P = P_mat*2.0d0

      print *, "Obtained New Density Matrix"
      print *, ""
      !print *, new_ P
      print *, ""
      !print *, "Printing Old Dens Mat"
      !print *, old_P

C     ______________compare new_P to old _P
      !m = compare(old_P,new_P,nbasis) 
      !print *, ""
      !print *, "m =",m
      !print *, ""
      !   if (m == 1) then
      !       print *, "New Density Matrix != Old Density Matrix"
      !   else 
      !       print *, "New Density Matrix = Old Density Matrix"
      !   end if
      print *, ""
        itr = itr +1
      print *,"***************************************************"
      print *,"***************************************************"

      end do  !ITERATION loop ends here
      print *, "#Total Iterations =",itr-1


C____________________________________________________
C     _______________MP2_____________
      print *, "___________  MP2 CORRECTION ______________"
      E_scf = E_old
      do mu = 1, nbasis
        do nu = 1, nbasis
          do lam = 1, nbasis
            do l = 1, nbasis
              do sig = 1, nbasis
               l_mo2eint(mu,nu,lam,l)= l_mo2eint(mu,nu,lam,l)
     &           + C(sig,l)*twoeint(mu,nu,lam,sig)
              end do
            end do
          end do
        end do
      end do

      do mu = 1, nbasis
        do nu = 1, nbasis
          do k = 1, nbasis
            do lam = 1, nbasis
              do l = 1, nbasis
               k_mo2eint(mu,nu,k,l)= k_mo2eint(mu,nu,k,l)
     &           + C(lam,k)*l_mo2eint(mu,nu,lam,l)
              end do
            end do
          end do
        end do
      end do
      !print *, "check2"


      do mu = 1, nbasis
        do j = 1, nbasis
          do nu = 1, nbasis
            do k = 1, nbasis
              do l = 1, nbasis
               j_mo2eint(mu,j,k,l)= j_mo2eint(mu,j,k,l)
     &           + C(nu,j)*k_mo2eint(mu,nu,k,l)
              end do
            end do
          end do
        end do
      end do


      do i = 1, nbasis
        do mu = 1, nbasis
          do j = 1, nbasis
            do k = 1, nbasis
              do l = 1, nbasis

               i_mo2eint(i,j,k,l)= i_mo2eint(i,j,k,l)
     &           + C(mu,i)*j_mo2eint(mu,j,k,l)

              end do
            end do
          end do
        end do
      end do

      !FINAL MO BASIS 2-e Integrals
      mo_2eint = i_mo2eint
      print *, "Obtained MO basis 2-e Integrals"
      E_mp2 = 0.0d0

      ! MP2 Energy Correction Calculation
      do i = 1, nocc
        do j = 1, nocc
          do a = nocc+1,nbasis
            do b = nocc+1, nbasis
       sum_mp2 = 0.0d0
       E_mp2=(mo_2eint(i,a,j,b)*(2*mo_2eint(i,a,j,b)-mo_2eint(i,b,j,a))/
     & (orb_E(i)+orb_E(j)-orb_E(a)-orb_E(b)))
       !print *, "sum_mp2=", sum_mp2
     & + E_mp2

       !sum_mp2=mo_2eint(i,a,j,b)*(2*mo_2eint(i,a,j,b)-mo_2eint(i,b,j,a))
       !E_mp2 = E_mp2 +
      !& (sum_mp2/(orb_E(i)+orb_E(j)-orb_E(a)-orb_E(b)))
       !print *, "E_mp2=", E_mp2
            end do
          end do
        end do
      end do

      print *, ""
      print *, "E_mp2=", E_mp2
      print *, ""
      print *, "E_tot=",E_mp2+E_scf

      print *, ""
      print *, "orb E"
      print *, orb_E




      call aces_ja_fin
      Call aces_fin

      Stop
      End





      function compare(a_mat,b_mat,nbas) result(d)
      Implicit None
      Integer :: nbas, d
      Double precision, DIMENSION(nbas,nbas):: a_mat
      Double precision, DIMENSION(nbas,nbas):: b_mat
      REAL :: diff = 0.0
      Integer ::p,q
      d = 0
      ilp: do p = 1,nbas
         jlp: do q = 1,nbas
            diff = abs(a_mat(p,q)-b_mat(p,q))
            do while (diff > 0.0000001)
               d= 1
                EXIT
            end do
            !else return 0
         end do jlp
      end do ilp
!------------------------------------------------------------------------------------------------------------------------------------------- 

