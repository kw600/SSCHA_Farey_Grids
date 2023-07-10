 MODULE constants
!--------------------------------------------------------------!
! Numerical constants and constants for variable declarations. !
!--------------------------------------------------------------!
 IMPLICIT NONE

 INTEGER,PARAMETER :: dp=kind(1.d0)
 REAL(dp),PARAMETER :: two_pi=8.d0*atan(1.d0)

 END MODULE constants


 MODULE utils
!-------------------!
! Various utilites. !
!-------------------!
 USE constants

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE errstop(sub,message)
!-------------------------------------------!
! Report an error in a subroutine and stop. !
!-------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: sub,message

 write(*,*)
 write(*,*)'Error in subroutine '//trim(adjustl(sub))//'.'
 write(*,*)
 call wordwrap(trim(adjustl(message)))
 write(*,*)
 stop

 END SUBROUTINE errstop


 SUBROUTINE erralloc(arg)
!--------------------------------------!
! Report an allocation error and stop. !
!--------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: arg
 
 write(*,*)
 write(*,*)'Problem allocating '//trim(adjustl(arg))//' array.'
 write(*,*)

 END SUBROUTINE erralloc


 CHARACTER(12) FUNCTION i2s(n)
!------------------------------------------------------------------------------!
! Convert integers to left justified strings that can be printed in the middle !
! of a sentence without introducing large amounts of white space.              !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 INTEGER :: i,j
 INTEGER,PARAMETER :: ichar0=ichar('0')

 i2s=''
 i=abs(n)

 do j=len(i2s),1,-1
  i2s(j:j)=achar(ichar0+mod(i,10))
  i=i/10
  if(i==0)exit
 enddo ! j

 if(n<0)then
  i2s='-'//adjustl(i2s)
 else
  i2s=adjustl(i2s)
 endif ! n<0

 END FUNCTION i2s


 SUBROUTINE wordwrap(text)
!-------------------------------------------------------------------------!
! This subroutine prints out the contents of the character string 'text', !
! ensuring that line breaks only occur at space characters.               !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: text
 CHARACTER(260) :: temp
 INTEGER :: i,unit,lentext,startpoint,stoppoint,lastpos,linelength

 unit=6

 lentext=len(trim(text))
 if(lentext<1)then
  write(unit,*)
  return
 endif

 linelength=79
 startpoint=1

 do i=1,huge(1)
  stoppoint=startpoint+linelength-1
  if(stoppoint<=lentext)then
   lastpos=index(trim(text(startpoint:stoppoint))," ",.true.)
   if(lastpos>0)stoppoint=startpoint+lastpos-1
  else
   stoppoint=lentext
  endif ! stoppoint<=lentext
  if(i==1)then
   temp=text(startpoint:stoppoint)
   write(unit,*)trim(temp)
  else
   temp=text(startpoint:stoppoint)
   write(unit,*)trim(adjustl(temp))
  endif ! i==1
  if(stoppoint==lentext)then
   exit
  else
   startpoint=stoppoint+1
  endif ! stoppoint==lentext
 enddo ! i
 
 END SUBROUTINE wordwrap


 SUBROUTINE read_input_files(no_symm_ops,basis,grid,prim_latt_vecs,symm_ops)
!-------------------------!
! Read basic input files. !
!-------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_symm_ops
 INTEGER,INTENT(out) :: basis,grid(3)
 REAL(dp),INTENT(out) :: prim_latt_vecs(3,3),symm_ops(4,3,no_symm_ops)
 INTEGER :: ierr,i_symm,i_row

! Read basis.dat file
 open(unit=11,file='equilibrium.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening basis.dat file.')
 read(11,*,iostat=ierr)basis
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading basis.dat file.')
 close(11)

! Read grid.dat file
 open(unit=11,file='grid.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening grid.dat file.')
 read(11,*,iostat=ierr)grid(1:3)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading grid.dat file.')
 close(11)

! Read prim.dat file
 open(unit=11,file='lattice.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening prim.dat file.')
 read(11,*,iostat=ierr)prim_latt_vecs(1,1:3)
 if(ierr==0)read(11,*,iostat=ierr)prim_latt_vecs(2,1:3)
 if(ierr==0)read(11,*,iostat=ierr)prim_latt_vecs(3,1:3)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading prim.dat file.')
 close(11)

! Read symmetry.dat file
 open(unit=11,file='symmetry.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening symmetry.dat &
  &file.')
 read(11,*,iostat=ierr)
 do i_symm=1,no_symm_ops
  do i_row=1,4
   read(11,*,iostat=ierr)symm_ops(i_row,1:3,i_symm)
   if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading symmetry.dat &
    &file.')
  enddo ! i_row
 enddo ! i_symm
 close(11)

 END SUBROUTINE read_input_files


 SUBROUTINE read_kpoints(no_kpoints,kpoints,multiplicity,kpoint_to_supercell)
!--------------------------------------------------!
! Read input files related to k-points in the IBZ. !
!--------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-8
 INTEGER,INTENT(in) :: no_kpoints
 INTEGER,INTENT(out) :: multiplicity(no_kpoints),&
  &kpoint_to_supercell(no_kpoints)
 REAL(dp),INTENT(out) :: kpoints(3,no_kpoints)
 INTEGER :: ierr,i_point
 REAL(dp) :: kpoint(3)

! Read ibz.dat file
 open(unit=11,file='ibz.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_KPOINTS','Problem opening ibz.dat file.')
 do i_point=1,no_kpoints
  read(11,*,iostat=ierr)kpoints(1:3,i_point),multiplicity(i_point)
  if(ierr/=0)call errstop('READ_KPOINTS','Problem reading ibz.dat file.')
 enddo ! i_point
 close(11)

! Read kpoint_to_supercell.dat file
 open(unit=11,file='kpoint_to_supercell.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_KPOINTS','Problem opening &
  &kpoint_to_supercell.dat file.')
 do i_point=1,no_kpoints
  read(11,*,iostat=ierr)kpoint(1:3),kpoint_to_supercell(i_point)
  if(ierr/=0)call errstop('READ_KPOINTS','Problem reading &
   &kpoint_to_supercell.dat file.')
  if(any(abs(kpoint(1:3)-kpoints(1:3,i_point))>tol))then
   call errstop('READ_KPOINTS','k-points in ibz.dat file and &
    &kpoints_to_supercell.dat file disagree.')
  endif ! tol
 enddo ! i_point
 close(11)
 
! Check k-points are in expected order
! do i_point=2,no_kpoints
!  if(.not.(kpoint_to_supercell(i_point)==kpoint_to_supercell(i_point-1)).and.&
!   &.not.(kpoint_to_supercell(i_point)==kpoint_to_supercell(i_point-1)+1))then
!   call errstop('READ_KPOINTS','k-points are not in expected order.')
!  endif ! kpoint_to_supercell
! enddo ! i_point

 do i_point=1,no_kpoints
  kpoints(1:3,i_point)=modulo(kpoints(1:3,i_point)+0.5d0+tol,1.d0)-0.5d0-tol
 enddo ! i_points

 END SUBROUTINE read_kpoints


 SUBROUTINE read_dyn_mats(basis,mass,species,atom_prim_frac,no_kpoints,dyn_mats,&
  &kpoint_to_supercell)
!--------------------------------------------------------!
! Read in dynamical matrices at each k-point in the IBZ. !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_kpoints,kpoint_to_supercell(no_kpoints)
 REAL(dp),INTENT(out) :: mass(basis),atom_prim_frac(3,basis)
 CHARACTER(2) :: species(basis),dummy_species
 COMPLEX(dp),INTENT(out) :: dyn_mats(basis,3,basis,3,no_kpoints)
 REAL(dp),PARAMETER :: mass_tol=1.d-4,frac_tol=1.d-8
 INTEGER :: ierr,i_atom,i_cart,j_atom,j_cart,atom1,cart1,atom2,cart2,ibz_point,&
  &supercell,atom_map(basis)
 REAL(dp) :: real_part,imag_part,temp_mass,temp_frac(3)
 LOGICAL :: found_atom(basis)

 open(unit=11,file='atoms_in_primitive_cell.1.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening &
  &atoms_in_primitive_cell.1.dat file.')
 do i_atom=1,basis
  read(11,*,iostat=ierr)species(i_atom),mass(i_atom),atom_prim_frac(1:3,i_atom)
  if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading &
  &atoms_in_primitive_cell.1.dat file.')
 enddo ! i_atom
 if(any(atom_prim_frac(1:3,1:basis)<0.d0.or.&
  &atom_prim_frac(1:3,1:basis)>=1.d0))call errstop('READ_DYN_MATS',&
   &'Fractional atomic coordinates are not in range [0.0,1.0)')
 close(11)

 open(unit=11,file='dyn_mat.1.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening dyn_mat.1.dat file.')
 do i_atom=1,basis
  do i_cart=1,3
   do j_atom=1,basis
    do j_cart=1,3
     read(11,*,iostat=ierr)atom1,cart1,atom2,cart2,real_part,imag_part
     if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading dyn_mat.1.dat &
      &file.')
     if(atom1/=i_atom.or.cart1/=i_cart.or.atom2/=j_atom.or.cart2/=j_cart)call &
      errstop('READ_DYN_MATS','dyn_mat.1.dat file does not seem to be in the &
       &expected order.')
      dyn_mats(atom1,cart1,atom2,cart2,1)=cmplx(real_part,imag_part,dp)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom
 close(11)

 do ibz_point=2,no_kpoints
  supercell=kpoint_to_supercell(ibz_point)
  open(unit=11,file='atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat',&
   &status='old',iostat=ierr)
  if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening &
   &atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat file.')
  found_atom(1:basis)=.false.
  atom_map(1:basis)=0
  do i_atom=1,basis
   read(11,*,iostat=ierr)dummy_species,temp_mass,temp_frac(1:3)
   if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading &
    &atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat file.')
   do j_atom=1,basis
    if(abs(temp_mass-mass(j_atom))<mass_tol)then
     if(all(abs(temp_frac(1:3)-atom_prim_frac(1:3,j_atom))<frac_tol))then
      found_atom(j_atom)=.true.
      atom_map(i_atom)=j_atom
     endif ! frac_tol
    endif ! mass_tol
   enddo ! j_atom
  enddo ! i_atom
  if(.not.any(found_atom(1:basis)))call errstop('READ_DYN_MATS','Unable to &
   &find all atoms in supercell '//trim(i2s(supercell))//'.')
  close(11)
  open(unit=11,file='dyn_mat.'//trim(i2s(ibz_point))//'.dat',status='old',&
   &iostat=ierr)
  if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening dyn_mat.'//&
   &trim(i2s(ibz_point))//'.dat file.')
  do i_atom=1,basis
   do i_cart=1,3
    do j_atom=1,basis
     do j_cart=1,3
      read(11,*,iostat=ierr)atom1,cart1,atom2,cart2,real_part,imag_part
      if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading dyn_mat.'//&
       &trim(i2s(ibz_point))//'.dat file.')
      if(atom1/=i_atom.or.cart1/=i_cart.or.atom2/=j_atom.or.cart2/=j_cart)call &
       errstop('READ_DYN_MATS','dyn_mat.'//trim(i2s(ibz_point))//'.dat file &
        &does not seem to be in the expected order.')
       dyn_mats(atom_map(atom1),cart1,atom_map(atom2),cart2,ibz_point)=&
        &cmplx(real_part,imag_part,dp)
     enddo ! j_cart
    enddo ! j_atom
   enddo ! i_cart
  enddo ! i_atom
  close(11)
 enddo ! ibz_point

 END SUBROUTINE read_dyn_mats


 SUBROUTINE read_path(no_points,path)
!-----------------------------------------------------------------!
! Read in the high symmetry points on the phonon dispersion path. !
!-----------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_points
 REAL(dp),INTENT(out) :: path(3,no_points)
 INTEGER :: ierr,i_point

 open(unit=11,file='path.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_PATH','Problem opening path.dat file.')
 do i_point=1,no_points
  read(11,*,iostat=ierr)path(1:3,i_point)
  if(ierr/=0)call errstop('READ_PATH','Problem reading path.dat file.')
 enddo ! i_point
 close(11)

 END SUBROUTINE read_path

 END MODULE utils


 MODULE linear_algebra
!----------------!
! LINEAR_ALGEBRA !
!----------------!
 USE constants
 USE utils

 IMPLICIT NONE

 INTERFACE

 SUBROUTINE zcopy(N,ZX,INCX,ZY,INCY)
  INTEGER,INTENT(in) :: INCX,INCY,N
  COMPLEX(KIND(1.d0)),INTENT(in) :: ZX(*)
  COMPLEX(KIND(1.d0)),INTENT(out) :: ZY(*)
 END SUBROUTINE zcopy

 SUBROUTINE zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
  CHARACTER(1),INTENT(in) :: JOBZ,UPLO
  INTEGER,INTENT(out) :: INFO
  INTEGER,INTENT(in) :: LDA,LWORK,N
  REAL(KIND(1.d0)),INTENT(out) :: W(*)
  REAL(KIND(1.d0)),INTENT(inout) :: RWORK(*)
  COMPLEX(KIND(1.d0)),INTENT(inout) :: A(LDA,*),WORK(*)
 END SUBROUTINE zheev

 END INTERFACE

 CONTAINS

 REAL(dp) FUNCTION det_33(A)
!-----------------------------------------------------!
! Given a 3x3 matrix A, this function returns det(A). !
!-----------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: A(3,3)

 det_33=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))+&
  &A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

 END FUNCTION det_33


 SUBROUTINE inv_33(A,B)
!-------------------------------------------------------------------------!
! This subroutine calculates the inverse B of matrix A, where A and B are !
! real 3x3 matrices.                                                      !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: A(3,3)
 REAL(dp),INTENT(out) :: B(3,3)
 REAL(dp) :: d
 
 d=det_33(A)
 if(d==0.d0)call errstop('INV_33','Trying to invert a singular matrix.')  
 d=1.d0/d

 B(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))*d
 B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))*d
 B(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))*d
 B(2,1)=(A(3,1)*A(2,3)-A(2,1)*A(3,3))*d
 B(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))*d
 B(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))*d
 B(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))*d
 B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))*d
 B(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))*d

 END SUBROUTINE inv_33

 END MODULE linear_algebra


 MODULE minimum_image
!---------------!
! MINIMUM_IMAGE !
!---------------!
 USE constants
 USE utils
 USE linear_algebra

 IMPLICIT NONE

 CONTAINS 

 SUBROUTINE min_images_brute_force(a,lat_vec,b,nim)
!---------------------------------------------------------------------!
! Compute the minimum image vectors b of vector a with respect to the !
! lattice specified by the rows of lat_vec.                           !
!---------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: a(3),lat_vec(3,3)
 REAL(dp),INTENT(out) :: b(3,8)
 INTEGER,INTENT(out) :: nim
 REAL(dp) :: rec_vec(3,3),delta1(3),delta2(3),delta3(3),mag_b_sq,dist2,tol_L2
 INTEGER :: n(3),i,j,k
 INTEGER,PARAMETER :: check_shell=3
 REAL(dp),PARAMETER :: tol=1.d-8

 call inv_33(lat_vec,rec_vec)
 rec_vec=transpose(rec_vec)

 tol_L2=tol*dot_product(lat_vec(1,1:3),lat_vec(1,1:3))
 n(1)=floor(dot_product(a(1:3),rec_vec(1,1:3)))
 n(2)=floor(dot_product(a(1:3),rec_vec(2,1:3)))
 n(3)=floor(dot_product(a(1:3),rec_vec(3,1:3)))

 mag_b_sq=-1.d0
 nim=-1

 do i=n(1)-check_shell,n(1)+check_shell+1
  delta1=a-dble(i)*lat_vec(1,1:3)
  do j=n(2)-check_shell,n(2)+check_shell+1
   delta2=delta1-dble(j)*lat_vec(2,1:3)
   do k=n(3)-check_shell,n(3)+check_shell+1
    delta3=delta2-dble(k)*lat_vec(3,1:3)
    dist2=dot_product(delta3,delta3)
    if(abs(dist2-mag_b_sq)<=tol_L2)then
     nim=nim+1
     if(nim>8)call errstop('MIN_IMAGES_BRUTE_FORCE','Need to increase &
      &maxim parameter.')
     b(1:3,nim)=delta3(1:3)
    elseif(dist2<mag_b_sq.or.nim==-1)then
     mag_b_sq=dist2
     nim=1
     b(1:3,1)=delta3(1:3)
    endif
   enddo ! k
  enddo ! j
 enddo ! i
 
 if(nim<=0)call errstop('MIN_IMAGES_BRUTE_FORCE','Bug.')
 
 END SUBROUTINE min_images_brute_force

 END MODULE minimum_image


 MODULE symmetry
!----------!
! SYMMETRY !
!----------!
 USE constants
 USE utils
 USE linear_algebra

 IMPLICIT NONE

 CONTAINS

 LOGICAL FUNCTION lattice_point(cart,latt_vecs)
!------------------------------------------------------------------------------!
! Given a vector with Cartesian coordinates cart, this function is true if the !
! vector is a point on the lattice specified by latt_vecs.                     !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-4
 REAL(dp),INTENT(in) :: cart(3),latt_vecs(3,3)
 REAL(dp) :: frac(3),rec_vecs(3,3)

 call inv_33(latt_vecs,rec_vecs)
 rec_vecs=transpose(rec_vecs)

 frac(1:3)=matmul(rec_vecs(1:3,1:3),cart(1:3))

 frac(1:3)=modulo(frac(1:3)+tol,1.d0)-tol

 lattice_point=all(abs(frac(1:3))<tol)

 END FUNCTION lattice_point


 SUBROUTINE kpoint_symmetry_maps(rec_vecs,no_points,points_cart,&
  &no_symms,point_symms,forwards,backwards)
!-----------------------------------------------------------------------------!
! Determine the mapping under each symmetry operation for each k point on the !
! grid.                                                                       !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-8
 INTEGER,INTENT(in) :: no_points,no_symms
 REAL(dp),INTENT(in) :: rec_vecs(3,3),points_cart(3,no_points),&
  &point_symms(3,3,no_symms)
 INTEGER,INTENT(out) :: forwards(no_points,no_symms),&
  &backwards(no_points,no_symms)
 INTEGER :: i_grid,i_symm,j_grid
 REAL(dp) :: rot_point_cart(3)

 forwards(1:no_points,1:no_symms)=0
 backwards(1:no_points,1:no_symms)=0

 do i_symm=1,no_symms
  do i_grid=1,no_points
   rot_point_cart(1:3)=matmul(point_symms(1:3,1:3,i_symm),&
    &points_cart(1:3,i_grid))
   do j_grid=1,no_points
    if(lattice_point(rot_point_cart(1:3)-points_cart(1:3,j_grid),&
     &rec_vecs(1:3,1:3)))then
     if(forwards(i_grid,i_symm)/=0)then
      call errstop('KPOINT_SYMMETRY_GROUPS','Grid point '//trim(i2s(i_grid))//&
       &' is transformed to more than one grid point by symmetry operation '&
       &//trim(i2s(i_symm))//'.')
     else
      forwards(i_grid,i_symm)=j_grid
     endif ! forwards
     if(backwards(j_grid,i_symm)/=0)then
      call errstop('KPOINT_SYMMETRY_GROUPS','More than one grid point is &
       &transformed to grid point '//trim(i2s(i_grid))//' by symmetry &
       &operation '//trim(i2s(i_symm))//'.')
     else
      backwards(j_grid,i_symm)=i_grid
     endif ! backwards
    endif ! lattice_point
   enddo ! j_grid
  enddo ! i_grid
 enddo ! i_symm

 END SUBROUTINE kpoint_symmetry_maps


 SUBROUTINE atom_symmetry_maps(latt_vecs,basis,atom_cart,no_symms,point_symms,&
  &trans_symms,forwards,backwards)
!---------------------------------------------------------------------!
! Construct the mapping of each atom in the primitive cell under each !
! symmetry operation.                                                 !
!---------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_symms
 REAL(dp),INTENT(in) :: latt_vecs(3,3),atom_cart(3,basis),&
  &point_symms(3,3,no_symms),trans_symms(3,no_symms)
 INTEGER,INTENT(out) :: forwards(basis,no_symms),backwards(basis,no_symms)
 INTEGER :: i_symm,i_atom,j_atom
 REAL(dp) :: symm_pos(3)
 LOGICAL :: found_atom(basis)

 forwards=0
 backwards=0

 do i_symm=1,no_symms
  found_atom=.false.
  do i_atom=1,basis
   symm_pos(1:3)=matmul(point_symms(1:3,1:3,i_symm),atom_cart(1:3,i_atom))+&
    &trans_symms(1:3,i_symm)
   do j_atom=1,basis
    if(lattice_point(symm_pos-atom_cart(1:3,j_atom),latt_vecs))then
     found_atom(i_atom)=.true.
     if(forwards(i_atom,i_symm)/=0)then
      call errstop('ATOM_SYMMETRY_MAPS','Atom '//trim(i2s(i_atom))//' is &
       &transformed to more than one atom by symmetry operation '&
       &//trim(i2s(i_symm))//'.')
     else
      forwards(i_atom,i_symm)=j_atom
     endif ! forwards
     if(backwards(j_atom,i_symm)/=0)then
      call errstop('ATOM_SYMMETRY_MAPS','More than one atom is mapped to atom '&
       &//trim(i2s(j_atom))//' by symmetry operation '//trim(i2s(i_symm))//'.')
     else
      backwards(j_atom,i_symm)=i_atom
     endif ! backwards
    endif ! lattice_point
   enddo ! j_atom
  enddo ! i_atom
  if(any(.not.found_atom))call errstop('ATOM_SYMMETRY_MAPS','Unable to &
   &map all atoms under symmetry operation '//trim(i2s(i_symm))//'.')
 enddo ! i_symm

 END SUBROUTINE atom_symmetry_maps


 SUBROUTINE match_kpoints(rec_latt_vecs,no_grid_points,grid_points_cart,&
  &no_ibz_points,ibz_points_cart,no_symms,point_symms,map_ibz,map_symm,&
  &time_reversal)
!-------------------------------------------------!
! Map each k point on the grid to one in the IBZ. !
!-------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_grid_points,no_ibz_points,no_symms
 REAL(dp),INTENT(in) :: rec_latt_vecs(3,3),grid_points_cart(3,no_grid_points),&
  &ibz_points_cart(3,no_ibz_points),point_symms(3,3,no_symms)
 INTEGER,INTENT(out) :: map_ibz(no_grid_points),map_symm(no_grid_points)
 LOGICAL,INTENT(out) :: time_reversal(no_grid_points)
 INTEGER :: i_grid,i_point,i_symm
 REAL(dp) :: rot_ibz_point(3)
 LOGICAL :: found_grid_point(no_grid_points)

 found_grid_point=.false.

 do i_point=1,no_ibz_points
  do i_symm=1,no_symms
   rot_ibz_point(1:3)=matmul(point_symms(1:3,1:3,i_symm),&
    &ibz_points_cart(1:3,i_point))
   do i_grid=1,no_grid_points
    if(found_grid_point(i_grid))cycle
    if(lattice_point(rot_ibz_point(1:3)-grid_points_cart(1:3,i_grid),&
     &rec_latt_vecs(1:3,1:3)))then
     found_grid_point(i_grid)=.true.
     map_ibz(i_grid)=i_point
     map_symm(i_grid)=i_symm
     time_reversal(i_grid)=.false.
    elseif(lattice_point(rot_ibz_point(1:3)+grid_points_cart(1:3,i_grid),&
     &rec_latt_vecs(1:3,1:3)))then
     found_grid_point(i_grid)=.true.
     map_ibz(i_grid)=i_point
     map_symm(i_grid)=i_symm
     time_reversal(i_grid)=.true.      
    endif ! lattice_point
   enddo ! i_grid
  enddo ! i_symm
 enddo ! i_point

 if(any(.not.found_grid_point))call errstop('MATCH_KPOINTS','Unable to map &
  &all k points on grid to the IBZ.')

 END SUBROUTINE match_kpoints


 SUBROUTINE g_matrix_phases(kpoint_cart,point_symm,trans_symm,basis,atoms_cart,&
  &phase)
!------------------------------------------!
! Calculate phases used in gamma matrices. !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis
 REAL(dp),INTENT(in) :: kpoint_cart(3),point_symm(3,3),trans_symm(3),&
  &atoms_cart(3,basis)
 COMPLEX(dp),INTENT(out) :: phase(basis,basis)
 INTEGER :: i_atom,j_atom
 REAL(dp) :: rot_kpoint(3),symm_pos(3),arg

 rot_kpoint(1:3)=matmul(point_symm(1:3,1:3),kpoint_cart(1:3))

 do i_atom=1,basis
  do j_atom=1,basis
   symm_pos(1:3)=matmul(point_symm(1:3,1:3),atoms_cart(1:3,j_atom))+&
    &trans_symm(1:3)
   arg=dot_product(rot_kpoint(1:3),atoms_cart(1:3,i_atom)-symm_pos(1:3))
   phase(i_atom,j_atom)=cmplx(cos(arg),sin(arg),dp)
  enddo ! j_atom
 enddo ! i_atom

 END SUBROUTINE g_matrix_phases


 SUBROUTINE apply_symmetry_to_dyn_mat(basis,forwards,phase,dyn_mat_in,&
  &point_symm,dyn_mat_out)
!-----------------------------------------------------!
! Apply a symmetry operation to the dynamical matrix. !
!-----------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,forwards(basis)
 REAL(dp),INTENT(in) :: point_symm(3,3)
 COMPLEX(dp),INTENT(in) :: phase(basis,basis),dyn_mat_in(basis,3,basis,3)
 COMPLEX(dp),INTENT(out) :: dyn_mat_out(basis,3,basis,3)
 INTEGER :: i_atom,j_atom,symm_i_atom,symm_j_atom
 REAL(dp) :: trans_symm(3,3)
 COMPLEX(dp) :: temp_mat(3,3)

 trans_symm=transpose(point_symm)

 dyn_mat_out=cmplx(0.d0,0.d0,dp)

 do i_atom=1,basis
  symm_i_atom=forwards(i_atom)
  do j_atom=1,basis
  symm_j_atom=forwards(j_atom)
  temp_mat=dyn_mat_in(i_atom,1:3,j_atom,1:3)
  temp_mat=matmul(matmul(point_symm,temp_mat),trans_symm)
  dyn_mat_out(symm_i_atom,1:3,symm_j_atom,1:3)=phase(symm_j_atom,j_atom)*&
   &temp_mat(1:3,1:3)*conjg(phase(symm_i_atom,i_atom))
  enddo ! j_atom
 enddo ! i_atom

 END SUBROUTINE apply_symmetry_to_dyn_mat

 END MODULE symmetry


 MODULE phonon
!--------!
! PHONON !
!--------!
 USE constants
 USE utils
 USE linear_algebra
 USE symmetry

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,dyn_mat)
!--------------------------------------------------------------------------!
! Construct the dynamical matrix at an arbitrary wave vector. We take into !
! account all minimum image primitive cells in the supercell.              !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_cells,no_ims(basis,basis,no_cells)
 REAL(dp),INTENT(in) :: kpoint(3),mass(basis),cell_vecs(3,8,basis,basis,&
  &no_cells),force_consts(basis,3,basis,3,no_cells)
 COMPLEX(dp),INTENT(out) :: dyn_mat(basis,3,basis,3)
 INTEGER :: i_atom,j_atom,i_cart,j_cart,i_cell,i_im
 REAL(dp) :: k_dot_r,prefactor
 COMPLEX(dp) :: exp_i_k_dot_r
 
 dyn_mat=cmplx(0.d0,0.d0,dp)

 do i_atom=1,basis
  do i_cart=1,3
   do j_atom=1,basis
    do j_cart=1,3
     prefactor=1.d0/sqrt(mass(i_atom)*mass(j_atom))
     do i_cell=1,no_cells
      exp_i_k_dot_r=cmplx(0.d0,0.d0,dp)
      do i_im=1,no_ims(i_atom,j_atom,i_cell)
       k_dot_r=dot_product(kpoint(1:3),cell_vecs(1:3,i_im,i_atom,j_atom,i_cell))
       exp_i_k_dot_r=exp_i_k_dot_r+cmplx(cos(k_dot_r),-sin(k_dot_r),dp)
      enddo ! i_im
      exp_i_k_dot_r=exp_i_k_dot_r/dble(no_ims(i_atom,j_atom,i_cell))
      dyn_mat(i_atom,i_cart,j_atom,j_cart)=&
       &dyn_mat(i_atom,i_cart,j_atom,j_cart)+&
       &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)*exp_i_k_dot_r
     enddo ! i_cell
     dyn_mat(i_atom,i_cart,j_atom,j_cart)=prefactor*&
      &dyn_mat(i_atom,i_cart,j_atom,j_cart)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom

 END SUBROUTINE construct_dyn_mat


 SUBROUTINE calculate_frequencies(basis,dyn_mat,freqs)
!-----------------------------------------------------------------!
! Diagonalise the dynamical matrix and calculate its eigenvalues. !
!-----------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis
 COMPLEX(dp),INTENT(in) :: dyn_mat(basis,3,basis,3)
 REAL(dp),INTENT(out) :: freqs(3*basis)
 INTEGER :: i_atom,j_atom,i_cart,j_cart,i_index,j_index,info
 REAL(dp) :: rwork(9*basis-2),minus_freqs_sq(3*basis)
 COMPLEX(dp) :: temp,temp_mat(3*basis,3*basis),work(6*basis-1)

 i_index=0
 do i_atom=1,basis
  do i_cart=1,3
   i_index=i_index+1
   j_index=0
   do j_atom=1,basis
    do j_cart=1,3
     j_index=j_index+1
     temp_mat(i_index,j_index)=dyn_mat(i_atom,i_cart,j_atom,j_cart)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom

 do i_index=1,3*basis
  temp_mat(i_index,i_index)=cmplx(real(temp_mat(i_index,i_index)),0.d0,dp)
  do j_index=i_index+1,3*basis
   temp=0.5d0*(temp_mat(i_index,j_index)+conjg(temp_mat(j_index,i_index)))
   temp_mat(i_index,j_index)=temp
   temp_mat(j_index,i_index)=conjg(temp)
  enddo ! j_index
 enddo ! i_index

 call zheev('N','U',3*basis,temp_mat(1,1),3*basis,minus_freqs_sq(1),work(1),&
  &6*basis-1,rwork(1),info)
 if(info/=0)call errstop('CALCULATE_FREQUENCIES','ZHEEV failed.')

 i_index=3*basis
 do j_index=1,3*basis
  if(minus_freqs_sq(i_index)>=0.d0)then
   freqs(j_index)=-sqrt(minus_freqs_sq(i_index))
  else
   freqs(j_index)=sqrt(-minus_freqs_sq(i_index))
  endif ! minus_freqs_sq
  i_index=i_index-1
 enddo ! j_index

 ENDSUBROUTINE calculate_frequencies


 SUBROUTINE calculate_polarisation(basis,dyn_mat,freqs,vecs,reference)
!----------------------------------------------------------------------------------!
! Diagonalise the dynamical matrix and calculate its eigenvalues and eigenvectors. !
!----------------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,reference
 COMPLEX(dp),INTENT(in) :: dyn_mat(basis,3,basis,3)
 REAL(dp),INTENT(out) :: freqs(3*basis)
 COMPLEX(dp),INTENT(out) :: vecs(3*basis,3*basis)
 INTEGER :: i_atom,j_atom,i_cart,j_cart,i_index,j_index,info
 REAL(dp) :: rwork(9*basis-2),minus_freqs_sq(3*basis)
 COMPLEX(dp) :: temp,temp_mat(3*basis,3*basis),work(6*basis-1)

 i_index=0
 do i_atom=1,basis
  do i_cart=1,3
   i_index=i_index+1
   j_index=0
   do j_atom=1,basis
    do j_cart=1,3
     j_index=j_index+1
     temp_mat(i_index,j_index)=dyn_mat(i_atom,i_cart,j_atom,j_cart)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom

 do i_index=1,3*basis
  temp_mat(i_index,i_index)=cmplx(real(temp_mat(i_index,i_index)),0.d0,dp)
  do j_index=i_index+1,3*basis
   temp=0.5d0*(temp_mat(i_index,j_index)+conjg(temp_mat(j_index,i_index)))
   temp_mat(i_index,j_index)=temp
   temp_mat(j_index,i_index)=conjg(temp)
  enddo ! j_index
 enddo ! i_index

 if(reference==0)then
  do i_index=1,3*basis
   do j_index=i_index+1,3*basis
    temp_mat(i_index,j_index)=cmplx(real(temp_mat(i_index,j_index)),0.d0,dp)
    temp_mat(j_index,i_index)=cmplx(real(temp_mat(j_index,i_index)),0.d0,dp)
   enddo ! j_index
  enddo ! i_index
 endif ! reference

 call zheev('V','U',3*basis,temp_mat(1,1),3*basis,minus_freqs_sq(1),work(1),&
  &6*basis-1,rwork(1),info)
 if(info/=0)call errstop('CALCULATE_POLARISATION','ZHEEV failed.')

 i_index=3*basis
 do j_index=1,3*basis
  if(minus_freqs_sq(i_index)>=0.d0)then
   freqs(j_index)=-sqrt(minus_freqs_sq(i_index))
  else
   freqs(j_index)=sqrt(-minus_freqs_sq(i_index))
  endif ! minus_freqs_sq
  call zcopy(3*basis,temp_mat(1,i_index),1,vecs(1,j_index),1)
  i_index=i_index-1
 enddo ! j_index

 ENDSUBROUTINE calculate_polarisation


 SUBROUTINE generate_dispersion(rec_vecs,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,no_points,path)
!---------------------!
! GENERATE_DISPERSION !
!---------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_cells,no_ims(basis,basis,no_cells),no_points
 REAL(dp),INTENT(in) :: rec_vecs(3,3),mass(basis),cell_vecs(3,8,basis,basis,&
  &no_cells),force_consts(basis,3,basis,3,no_cells),path(3,no_points)
 INTEGER :: ialloc,ierr,i_path,path_length,i_point,j_point,i_dof
 REAL(dp) :: k_start(3),k_stop(3),k_diff(3),k_dist,total_k_dist,delta_k,&
  &kpoint(3),omega(3*basis)
 COMPLEX(dp) :: dyn_mat(basis,3,basis,3)

 total_k_dist=0.d0
 do i_point=1,no_points-1
  k_start(1:3)=matmul(path(1:3,i_point),rec_vecs(1:3,1:3))
  k_stop(1:3)=matmul(path(1:3,i_point+1),rec_vecs(1:3,1:3))
  k_diff(1:3)=k_stop(1:3)-k_start(1:3)
  total_k_dist=total_k_dist+sqrt(dot_product(k_diff(1:3),k_diff(1:3)))
 enddo ! i_point

 delta_k=total_k_dist/1000.d0

 open(unit=14,file='phonon_dispersion_curve.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISPERSION','Problem opening &
  &phonon_dispersion_curve.dat file')

 open(unit=15,file='high_symmetry_points.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISPERSION','Problem opening &
  &high_symmetry_points.dat file.')

 total_k_dist=0.d0
 do i_point=1,no_points-1
  write(15,*)i_point,total_k_dist
  k_start(1:3)=matmul(path(1:3,i_point),rec_vecs(1:3,1:3))
  k_stop(1:3)=matmul(path(1:3,i_point+1),rec_vecs(1:3,1:3))
  k_diff(1:3)=k_stop(1:3)-k_start(1:3)
  k_dist=sqrt(dot_product(k_diff(1:3),k_diff(1:3)))
  path_length=int(k_dist/delta_k)
  do i_path=0,path_length-1
   kpoint(1:3)=k_start(1:3)+dble(i_path)*k_diff(1:3)/dble(path_length)
   call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
    &force_consts,dyn_mat)
   call calculate_frequencies(basis,dyn_mat,omega)
   write(14,*)total_k_dist,omega
   total_k_dist=total_k_dist+k_dist/dble(path_length)
  enddo ! i_path
 enddo ! i_point
 call construct_dyn_mat(k_stop,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,dyn_mat)
 call calculate_frequencies(basis,dyn_mat,omega)
 write(14,*)total_k_dist,omega
 write(15,*)no_points,total_k_dist

 close(14)
 close(15)

 ENDSUBROUTINE generate_dispersion


 SUBROUTINE generate_dos(rec_vecs,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,temperature)
!--------------!
! GENERATE_DOS !
!--------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: no_bins=1000,no_prelims=10000,no_samples=10000,no_sets=100
 REAL(dp),PARAMETER :: freq_tol=1.d-6,safety_factor=1.1d0
 INTEGER,INTENT(in) :: basis,no_cells,no_ims(basis,basis,no_cells)
 REAL(dp),INTENT(in) :: rec_vecs(3,3),mass(basis),cell_vecs(3,8,basis,basis,&
  &no_cells),force_consts(basis,3,basis,3,no_cells),temperature
 INTEGER :: ialloc,ierr,i_sample,i_freq,i_bin,i_set,no_dof
 REAL(dp) :: max_freq,min_freq,frac(3),kpoint(3),freqs(3*basis),bin_width,&
  &rec_bin_width,freq_dos(no_bins,no_sets),free,error_free,square_free,omega,&
  &energy,square_energy,error_energy,bin_energy(no_bins),bin_free(no_bins),&
  &value_energy,value_free
 COMPLEX(dp) :: dyn_mat(basis,3,basis,3)
 LOGICAL :: soft_modes

 no_dof=3*basis
 soft_modes=.false.

! Initialise the random number generator
 call random_seed()

 max_freq=-1.d0
 min_freq=huge(1.d0)

 do i_sample=1,no_prelims
  call random_number(frac)
  kpoint(1:3)=matmul(frac(1:3),rec_vecs(1:3,1:3))
  call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
    &force_consts,dyn_mat)
  call calculate_frequencies(basis,dyn_mat,freqs)
  if(freqs(1)<min_freq)min_freq=freqs(1)
  if(freqs(no_dof)>max_freq)max_freq=freqs(no_dof)
 enddo ! i_sample
 soft_modes=(min_freq<-freq_tol)
 if(max_freq<=0.d0)call errstop('GENERATE_DOS','The system is pathologically &
  &unstable.')

 bin_width=safety_factor*max_freq/dble(no_bins)
 rec_bin_width=1.d0/bin_width
 freq_dos(1:no_bins,1:no_sets)=0.d0

 do i_set=1,no_sets
  write(*,*)'Calculating DOS set '//trim(i2s(i_set))//' of '//&
   &trim(i2s(no_sets))//'.'
  do i_sample=1,no_samples
   call random_number(frac(1:3))
   kpoint(1:3)=matmul(frac(1:3),rec_vecs(1:3,1:3))
   call construct_dyn_mat(kpoint(1:3),basis,mass(1:basis),no_cells,&
    &no_ims(1:basis,1:basis,1:no_cells),&
    &cell_vecs(1:3,1:8,1:basis,1:basis,1:no_cells),&
    &force_consts(1:basis,1:3,1:basis,1:3,1:no_cells),&
    &dyn_mat(1:basis,1:3,1:basis,1:3))
   call calculate_frequencies(basis,dyn_mat(1:basis,1:3,1:basis,1:3),&
    &freqs(1:no_dof))
   if(freqs(1)<-freq_tol)soft_modes=.true.
   do i_freq=1,no_dof
    if(freqs(i_freq)>freq_tol)then
     i_bin=ceiling(rec_bin_width*freqs(i_freq))
     if(i_bin>no_bins)call errstop('GENERATE_DOS','Frequency too high to be &
      &binned.')
     freq_dos(i_bin,i_set)=freq_dos(i_bin,i_set)+1.d0
    endif ! freqs
   enddo ! i_freq
  enddo ! i_sample
 enddo ! i_set

 if(soft_modes)write(*,*)'Soft modes present.'

 do i_bin=1,no_bins
  omega=bin_width*(dble(i_bin)-0.5d0)
  bin_energy(i_bin)=harmonic_energy(temperature,omega)
  bin_free(i_bin)=harmonic_free_energy(temperature,omega)
 enddo ! i_bin

 energy=0.d0 ; square_energy=0.d0 ; free=0.d0 ; square_free=0.d0

 do i_set=1,no_sets
  value_energy=dot_product(bin_energy(1:no_bins),freq_dos(1:no_bins,i_set))
  value_energy=value_energy/dble(no_samples)
  energy=energy+value_energy
  square_energy=square_energy+value_energy*value_energy
  value_free=dot_product(bin_free(1:no_bins),freq_dos(1:no_bins,i_set))
  value_free=value_free/dble(no_samples)
  free=free+value_free
  square_free=square_free+value_free*value_free
 enddo ! i_set

 energy=energy/dble(no_sets)
 free=free/dble(no_sets)
 square_energy=square_energy/dble(no_sets)
 square_free=square_free/dble(no_sets)
 error_energy=sqrt((square_energy-energy*energy)/dble(no_sets-1))
 error_free=sqrt((square_free-free*free)/dble(no_sets-1))

 open(unit=16,file='interpolated_energy.dat',position='append',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening &
  &interpolated_energy.dat file.')
 write(16,*)temperature,energy,error_energy
 close(16)

 open(unit=16,file='interpolated_free_energy.dat',position='append',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening &
  &interpolated_free_energy.dat file.')
 write(16,*)temperature,free,error_free
 close(16)

 freq_dos(1:no_bins,1:no_sets)=freq_dos(1:no_bins,1:no_sets)*rec_bin_width/&
  &dble(no_samples*no_sets)

 open(unit=16,file='freq_dos.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening freq_dos.dat file.')
 do i_bin=1,no_bins
  write(16,*)bin_width*(dble(i_bin)-0.5d0),sum(freq_dos(i_bin,1:no_sets))
 enddo ! i_bin
 close(16)

 END SUBROUTINE generate_dos


 SUBROUTINE generate_disp_patterns(no_prim_cells,cell_vecs,no_grid_points,&
  &grid_vecs,reference,basis,mass,dyn_mats)
!-------------------------!
! GENERATE_DISP_PATTERNS !
!-------------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol_omega=1.d-6
 INTEGER,INTENT(in) :: basis,no_grid_points,no_prim_cells,&
  &reference(no_grid_points)
 REAL(dp),INTENT(in) :: cell_vecs(3,no_prim_cells),grid_vecs(3,no_grid_points)
 COMPLEX(dp),INTENT(in) :: dyn_mats(basis,3,basis,3,no_grid_points)
 INTEGER :: ierr,i_atom,i_cart,i_cell,i_index,i_point,j_index
 REAL(dp) :: freqs(3*basis),disp_pattern(3),k_dot_r,mass(basis),prefactor
 COMPLEX(dp) :: exp_i_k_dot_r(no_prim_cells),pol_vec(3,basis),&
  &vecs(3*basis,3*basis)

 open(unit=101,file='disp_patterns.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISP_PATTERNS','Problem opening &
  &disp_patterns.dat file.')

 open(unit=102,file='kdisp_patterns.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISP_PATTERNS','Problem opening &
  &kdisp_patterns.dat file.')

 open(unit=103,file='pol_vec.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISP_PATTERNS','Problem opening pol_vec.dat &
  &file.')

 do i_point=1,no_grid_points
  call calculate_polarisation(basis,dyn_mats(1:basis,1:3,1:basis,1:3,i_point),&
   &freqs(1:3*basis),vecs(1:3*basis,1:3*basis),reference(i_point))
  do i_cell=1,no_prim_cells
   if(reference(i_point)==0)then
    k_dot_r=-dot_product(grid_vecs(1:3,i_point),cell_vecs(1:3,i_cell))
   elseif(reference(i_point)>i_point)then
    k_dot_r=-dot_product(grid_vecs(1:3,i_point),cell_vecs(1:3,i_cell))
   else
    k_dot_r=-dot_product(grid_vecs(1:3,reference(i_point)),&
     &cell_vecs(1:3,i_cell))
   endif ! reference
   exp_i_k_dot_r(i_cell)=cmplx(cos(k_dot_r),sin(k_dot_r),dp)
  enddo ! i_cell
  do j_index=1,3*basis
   i_index=0
   do i_atom=1,basis
    do i_cart=1,3
     i_index=i_index+1
     pol_vec(i_cart,i_atom)=vecs(i_index,j_index)
    enddo ! i_cart
   enddo ! i_atom
   if(freqs(j_index)<-tol_omega)then
    write(101,*)'Frequency : ',freqs(j_index),'(SOFT)' 
    write(102,*)'Frequency : ',freqs(j_index),'(SOFT)'
    write(103,*)'Mode number : ',j_index,'Frequency : ',freqs(j_index),'(SOFT)'
   else
    write(101,*)'Frequency : ',freqs(j_index)
    write(102,*)'Frequency : ',freqs(j_index)
    write(103,*)'Mode number : ',j_index,'Frequency : ',freqs(j_index)
   endif ! freqs
   write(101,*)grid_vecs(1:3,i_point)
   write(102,*)grid_vecs(1:3,i_point)
   write(103,*)grid_vecs(1:3,i_point)
   write(101,*)'Displacement pattern for each atom:'
   write(102,*)'Displacement pattern for each atom:'
   write(103,*)'Polarisation vector:'
   do i_atom=1,basis
     do i_cell=1,no_prim_cells
     if(reference(i_point)==0)then
      disp_pattern=real(pol_vec(1:3,i_atom)*exp_i_k_dot_r(i_cell))
      prefactor=1.d0
     elseif(reference(i_point)>i_point)then
      disp_pattern=real(pol_vec(1:3,i_atom)*exp_i_k_dot_r(i_cell))
      prefactor=sqrt(2.d0)
     else
      disp_pattern=aimag(pol_vec(1:3,i_atom)*exp_i_k_dot_r(i_cell))
      prefactor=sqrt(2.d0)     
     endif ! reference
     write(101,*)disp_pattern(1:3)/sqrt(mass(i_atom)),prefactor
     write(102,*)disp_pattern(1:3),prefactor
     write(103,*)real(pol_vec(1:3,i_atom))/sqrt(mass(i_atom))
     write(103,*)aimag(pol_vec(1:3,i_atom))/sqrt(mass(i_atom))
    enddo ! i_atom
   enddo ! i_cell
   write(101,*)
   write(102,*)
   write(103,*)
  enddo ! j_index
 enddo ! i_point

 close(101)
 close(102)
 close(103)

 END SUBROUTINE generate_disp_patterns


 SUBROUTINE eval_freqs_on_grid(basis,no_grid_points,dyn_mats,temperature)
!--------------------!
! EVAL_FREQS_ON_GRID !
!--------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_grid_points
 REAL(dp),INTENT(in) :: temperature
 COMPLEX(dp),INTENT(in) :: dyn_mats(basis,3,basis,3,no_grid_points)
 REAL(dp),PARAMETER :: tol_omega=1.d-8
 INTEGER :: ierr,i_dof,i_grid,no_dof
 REAL(dp) :: energy,free_energy,freqs(3*basis)

 no_dof=3*basis

 energy=0.d0
 free_energy=0.d0

 do i_grid=1,no_grid_points
  call calculate_frequencies(basis,dyn_mats(1:basis,1:3,1:basis,1:3,i_grid),&
   &freqs(1:no_dof))
  do i_dof=1,no_dof
   if(freqs(i_dof)>tol_omega)then
    energy=energy+harmonic_energy(temperature,freqs(i_dof))
    free_energy=free_energy+harmonic_free_energy(temperature,freqs(i_dof))
   endif ! freqs
  enddo ! i_dof
 enddo ! i_grid

 energy=energy/dble(no_grid_points)
 free_energy=free_energy/dble(no_grid_points)

 open(unit=16,file='grid_energy.dat',position='append',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening grid_energy.dat file.')
 write(16,*)temperature,energy
 close(16)

 open(unit=16,file='grid_free_energy.dat',position='append',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening grid_free_energy.dat &
  &file.')
 write(16,*)temperature,free_energy
 close(16)

 END SUBROUTINE eval_freqs_on_grid


 REAL(dp) FUNCTION harmonic_energy(temperature,omega)
!-----------------!
! HARMONIC_ENERGY !
!-----------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: omega,temperature
 REAL(dp),PARAMETER :: kB_au_per_K=3.16679002948702d-6,tol=1.d-8
 REAL(dp) :: denominator

 if(temperature<tol)then
  harmonic_energy=0.5d0*omega
 else
  denominator=exp(omega/(kB_au_per_K*temperature))-1.d0
  if(denominator>0.d0)then
   harmonic_energy=(1.d0/denominator+0.5d0)*omega
  else
   harmonic_energy=kB_au_per_K*temperature
  endif ! denominator
 endif ! temperature

 END FUNCTION harmonic_energy


 REAL(dp) FUNCTION harmonic_free_energy(temperature,omega)
!----------------------!
! HARMONIC_FREE_ENERGY !
!----------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: omega,temperature
 REAL(dp),PARAMETER :: kB_au_per_K=3.16679002948702d-6,tol=1.d-8
 REAL(dp) :: difference

 if(temperature<tol)then
  harmonic_free_energy=0.5d0*omega
 else
  difference=1.d0-exp(-omega/(kB_au_per_K*temperature))
  if(difference>0.d0)then
   harmonic_free_energy=0.5d0*omega+kB_au_per_K*temperature*log(difference)
  else
   harmonic_free_energy=-huge(0.d0)
  endif ! difference
 endif ! temperature

 END FUNCTION harmonic_free_energy

 END MODULE phonon


 PROGRAM fourier_interpolation
!-----------------------!
! FOURIER_INTERPOLATION !
!-----------------------!
 USE constants
 USE utils
 USE linear_algebra
 USE minimum_image
 USE symmetry
 USE phonon

 IMPLICIT NONE

 REAL(dp),PARAMETER :: tol=1.d-8 

 INTEGER,ALLOCATABLE :: atom_map_symm_backwards(:,:)
 INTEGER,ALLOCATABLE :: atom_map_symm_forwards(:,:)
 INTEGER,ALLOCATABLE :: grid_map_symm_backwards(:,:)
 INTEGER,ALLOCATABLE :: grid_map_symm_forwards(:,:)
 INTEGER,ALLOCATABLE :: grid_to_ibz_map(:)
 INTEGER,ALLOCATABLE :: ibz_to_grid_symm(:)
 INTEGER,ALLOCATABLE :: ibz_to_supercell_map(:)
 INTEGER,ALLOCATABLE :: identity_map(:)
 INTEGER,ALLOCATABLE :: multiplicity(:)
 INTEGER,ALLOCATABLE :: no_im_cells(:,:,:)
 INTEGER,ALLOCATABLE :: reference(:)

 REAL(dp),ALLOCATABLE :: atom_prim_cart(:,:)
 REAL(dp),ALLOCATABLE :: atom_prim_frac(:,:)
 REAL(dp),ALLOCATABLE :: atom_super_cart(:,:,:)
 REAL(dp),ALLOCATABLE :: cell_pos_cart(:,:)
 REAL(dp),ALLOCATABLE :: delta_prim(:,:,:,:,:)
 REAL(dp),ALLOCATABLE :: force_consts(:,:,:,:,:)
 REAL(dp),ALLOCATABLE :: grid_points_cart(:,:)
 REAL(dp),ALLOCATABLE :: grid_points_frac(:,:)
 REAL(dp),ALLOCATABLE :: ibz_points_cart(:,:)
 REAL(dp),ALLOCATABLE :: ibz_points_frac(:,:)
 REAL(dp),ALLOCATABLE :: mass(:)
 REAL(dp),ALLOCATABLE :: path(:,:)
 REAL(dp),ALLOCATABLE :: point_symms(:,:,:)
 REAL(dp),ALLOCATABLE :: symm_ops(:,:,:)
 REAL(dp),ALLOCATABLE :: trans_symms(:,:)

 CHARACTER(2),ALLOCATABLE :: species(:)

 COMPLEX(dp),ALLOCATABLE :: dyn_mats_grid(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: dyn_mats_ibz(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: phase(:,:)
 COMPLEX(dp),ALLOCATABLE :: temp_dyn_mat(:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: dyn_mats_symm(:,:,:,:,:)

 LOGICAL,ALLOCATABLE :: time_reversal(:)
 LOGICAL,ALLOCATABLE :: two_k_is_a_G(:)

 INTEGER :: task

 INTEGER :: basis
 INTEGER :: grid(3)
 INTEGER :: dof_prim
 INTEGER :: no_grid_points
 INTEGER :: no_ibz_points
 INTEGER :: no_kpoints_path
 INTEGER :: no_prim_cells
 INTEGER :: no_symm_ops

 INTEGER :: ialloc
 INTEGER :: ierr
 INTEGER :: istat

 INTEGER :: counter
 INTEGER :: i_atom,j_atom
 INTEGER :: i_back
 INTEGER :: i_cart,j_cart
 INTEGER :: i_cell
 INTEGER :: i_grid,j_grid
 INTEGER :: i_im
 INTEGER :: i_point
 INTEGER :: i_symm
 INTEGER :: m1,m2,m3

 REAL(dp) :: delta_r_corr(3)
 REAL(dp) :: delta_r_ims(3,8)
 REAL(dp) :: identity(3,3)
 REAL(dp) :: k_dot_r
 REAL(dp) :: kpoint(3)
 REAL(dp) :: prefactor
 REAL(dp) :: prim_latt_vecs(3,3)
 REAL(dp) :: prim_rec_vecs(3,3)
 REAL(dp) :: super_latt_vecs(3,3)
 REAL(dp) :: super_rec_vecs(3,3)
 REAL(dp) :: temperature

 COMPLEX(dp) :: exp_i_k_dot_r

! Get total number of symmetry operations and allocate corresponding arrays
 open(unit=10,file='symmetry.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening symmetry.dat file in main program.'
  stop
 endif ! ierr
 read(10,*,iostat=ierr)no_symm_ops
 if(ierr/=0)then
  write(*,*)'Problem reading symmetry.dat file in main program.'
  stop
 endif ! ierr
 close(10)
 
 allocate(symm_ops(4,3,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('SYMM_OPS')
 allocate(point_symms(3,3,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('POINT_SYMMS')
 allocate(trans_symms(3,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('TRANS_SYMMS')

! Read basic input files and allocate corresponding arrays
 call read_input_files(no_symm_ops,basis,grid,prim_latt_vecs,symm_ops)

 dof_prim=3*basis
 no_grid_points=product(grid(1:3))
 no_prim_cells=no_grid_points

 do i_symm=1,no_symm_ops
  point_symms(1:3,1:3,i_symm)=symm_ops(1:3,1:3,i_symm)
  trans_symms(1:3,i_symm)=matmul(symm_ops(4,1:3,i_symm),prim_latt_vecs(1:3,1:3))
 enddo ! i_symm

 allocate(identity_map(basis),stat=ialloc)
 if(ialloc/=0)call erralloc('IDENTITY_MAP')

 allocate(grid_to_ibz_map(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_TO_IBZ_MAP')

 allocate(ibz_to_grid_symm(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_TO_GRID_SYMM')

 allocate(reference(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('REFERENCE')

 allocate(time_reversal(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('TIME_REVERSAL')

 allocate(two_k_is_a_G(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('TWO_K_IS_A_G')

 allocate(grid_map_symm_backwards(no_grid_points,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_MAP_SYMM_BACKWARDS')

 allocate(grid_map_symm_forwards(no_grid_points,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_MAP_SYMM_FORWARDS')

 allocate(grid_points_cart(3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_POINTS_CART')

 allocate(grid_points_frac(3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_POINTS_FRAC')

 allocate(dyn_mats_grid(basis,3,basis,3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('DYN_MATS_GRID')

 allocate(temp_dyn_mat(basis,3,basis,3),stat=ialloc)
 if(ialloc/=0)call erralloc('TEMP_DYN_MAT')

 allocate(dyn_mats_symm(basis,3,basis,3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('DYN_MATS_SYMM')

 allocate(force_consts(basis,3,basis,3,no_prim_cells),stat=ialloc)
 if(ialloc/=0)call erralloc('FORCE_CONSTS')

 identity=0.d0
 identity(1,1)=1.d0
 identity(2,2)=1.d0
 identity(3,3)=1.d0

 do i_atom=1,basis
  identity_map(i_atom)=i_atom
 enddo ! i_atom

 call inv_33(prim_latt_vecs,prim_rec_vecs)
 prim_rec_vecs=transpose(prim_rec_vecs)
 prim_rec_vecs=two_pi*prim_rec_vecs

 super_latt_vecs(1,1:3)=dble(grid(1))*prim_latt_vecs(1,1:3)
 super_latt_vecs(2,1:3)=dble(grid(2))*prim_latt_vecs(2,1:3)
 super_latt_vecs(3,1:3)=dble(grid(3))*prim_latt_vecs(3,1:3)

 call inv_33(super_latt_vecs,super_rec_vecs)
 super_rec_vecs=transpose(super_rec_vecs)
 super_rec_vecs=two_pi*super_rec_vecs

 i_grid=0
 do m1=0,grid(1)-1
  do m2=0,grid(2)-1
   do m3=0,grid(3)-1
    i_grid=i_grid+1
    if(i_grid>no_grid_points)then
     write(*,*)'Found more k-points than on grid.'
     stop
    endif ! i_grid
    grid_points_frac(1,i_grid)=dble(m1)/dble(grid(1))
    grid_points_frac(2,i_grid)=dble(m2)/dble(grid(2))
    grid_points_frac(3,i_grid)=dble(m3)/dble(grid(3))
    grid_points_frac(1:3,i_grid)=&
     &modulo(0.5d0+grid_points_frac(1:3,i_grid)+tol,1.d0)-0.5d0-tol
   enddo ! m3
  enddo ! m2
 enddo ! m1
 if(i_grid<no_grid_points)then
  write(*,*)'Not found all k-points on grid.'
  stop
 endif ! i_grid

! Convert grid points from fractional to Cartesian coodinates
 do i_grid=1,no_grid_points
  grid_points_cart(1:3,i_grid)=matmul(grid_points_frac(1:3,i_grid),&
   &prim_rec_vecs(1:3,1:3))
  two_k_is_a_G=lattice_point(2.d0*grid_points_cart(1:3,i_grid),&
   &prim_rec_vecs(1:3,1:3))
 enddo ! i_grid

! Determine +/- k-point pairs
 reference(1:no_grid_points)=0
 do i_grid=1,no_grid_points
  do j_grid=1,i_grid-1
   if(reference(j_grid)>0)cycle
   if(lattice_point(grid_points_cart(1:3,i_grid)+grid_points_cart(1:3,j_grid),&
    &prim_rec_vecs(1:3,1:3)))then
    reference(i_grid)=j_grid
    reference(j_grid)=i_grid
    exit
   endif ! lattice_point
  enddo ! j_grid
 enddo ! i_grid

! Determine which symmetry operations are in the point group and inverse group
! for each wave vector
 call kpoint_symmetry_maps(prim_rec_vecs,no_grid_points,grid_points_cart,&
  &no_symm_ops,point_symms,grid_map_symm_forwards,grid_map_symm_backwards)

 allocate(atom_prim_frac(3,basis),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_PRIM_FRAC')

 allocate(atom_prim_cart(3,basis),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_PRIM_CART')

 allocate(atom_super_cart(3,basis,no_prim_cells),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_SUPER_CART')

 allocate(mass(basis),stat=ialloc)
 if(ialloc/=0)call erralloc('MASS')

 allocate(species(basis),stat=ialloc)
 if(ialloc/=0)call erralloc('SPECIES')

 allocate(atom_map_symm_forwards(basis,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_MAP_SYMM_FORWARDS')

 allocate(atom_map_symm_backwards(basis,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_MAP_SYMM_BACKWARDS')

 allocate(phase(basis,basis),stat=ialloc)
 if(ialloc/=0)call erralloc('PHASE')

 allocate(cell_pos_cart(3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('CELL_POS_CART')

 i_cell=0
 do m1=0,grid(1)-1
  do m2=0,grid(2)-1
   do m3=0,grid(3)-1
    i_cell=i_cell+1
    if(i_cell>no_grid_points)then
     write(*,*)'Found more primitive cells than in supercell.'
     stop
    endif ! i_cell
    cell_pos_cart(1:3,i_cell)=dble(m1)*prim_latt_vecs(1,1:3)+&
     &dble(m2)*prim_latt_vecs(2,1:3)+dble(m3)*prim_latt_vecs(3,1:3)
   enddo ! m3
  enddo ! m2
 enddo ! m1
 if(i_cell<no_prim_cells)then
  write(*,*)'Not found all primitive cells in supercell.'
  stop
 endif ! i_cell

! Get the number of k-points in the IBZ and allocate corresponding arrays
 call system("echo $(wc -l ibz.dat | awk '{print $1}') > tempfile.dat",istat)
 if(istat/=0)then
  write(*,*)'Problem counting the number of lines in ibz.dat.'
  stop
 endif ! istat
 open(unit=10,file='tempfile.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening tempfile.dat file.'
  stop
 endif ! ierr
 read(10,*,iostat=ierr)no_ibz_points
 if(ierr/=0)then
  write(*,*)'Problem reading tempfile.dat file.'
  stop
 endif ! ierr
 close(10,status='delete')

 allocate(ibz_points_cart(3,no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_POINTS_CART')

 allocate(ibz_points_frac(3,no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_POINTS_FRAC')

 allocate(ibz_to_supercell_map(no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_TO_SUPERCELL_MAP')

 allocate(multiplicity(no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('MULTIPLICITY')

 allocate(dyn_mats_ibz(basis,3,basis,3,no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('DYN_MATS_IBZ')

! Read input files related to k-points in the IBZ
 call read_kpoints(no_ibz_points,ibz_points_frac,multiplicity,&
  &ibz_to_supercell_map)

! Convert IBZ points from fractional to Cartesian coodinates
 do i_point=1,no_ibz_points
  ibz_points_cart(1:3,i_point)=matmul(ibz_points_frac(1:3,i_point),&
   &prim_rec_vecs(1:3,1:3))
 enddo ! i_point

! Read in the dynamical matrix at each k-point in the IBZ
 call read_dyn_mats(basis,mass,species,atom_prim_frac,no_ibz_points,dyn_mats_ibz,&
  &ibz_to_supercell_map)

! Determine Cartesian coordinates of atoms in primitive cell
 do i_atom=1,basis
  atom_prim_cart(1:3,i_atom)=matmul(atom_prim_frac(1:3,i_atom),&
   &prim_latt_vecs(1:3,1:3))
 enddo ! i_atom

! Determine Cartesian coordinates of atoms in supercell
 do i_cell=1,no_prim_cells
  do i_atom=1,basis
   atom_super_cart(1:3,i_atom,i_cell)=atom_prim_cart(1:3,i_atom)+&
    &cell_pos_cart(1:3,i_cell)
  enddo ! i_atoms
 enddo ! i_cell

 open(1,file='super_equilibrium.dat')
 write(1,*)basis*no_prim_cells
 do i_atom=1,basis
  do i_cell=1,no_prim_cells
    write(1,*)species(i_atom),mass(i_atom),atom_super_cart(1:3,i_atom,i_cell)
  enddo ! i_cell
 enddo ! i_atom
 close(1)

 

 allocate(no_im_cells(basis,basis,no_prim_cells),stat=ialloc)
 if(ialloc/=0)call erralloc('NO_IM_CELLS')

 allocate(delta_prim(3,8,basis,basis,no_prim_cells),stat=ialloc)
 if(ialloc/=0)call erralloc('DELTA_PRIM')

 no_im_cells=0
 delta_prim=0.d0
 do i_atom=1,basis
  do j_atom=1,basis
   delta_r_corr(1:3)=atom_prim_cart(1:3,i_atom)-atom_prim_cart(1:3,j_atom)
   do i_cell=1,no_prim_cells
    call min_images_brute_force(atom_super_cart(1:3,j_atom,i_cell)-&
     &atom_prim_cart(1:3,i_atom),super_latt_vecs(1:3,1:3),delta_r_ims(1:3,1:8),&
     &no_im_cells(i_atom,j_atom,i_cell))
    do i_im=1,no_im_cells(i_atom,j_atom,i_cell)
     delta_prim(1:3,i_im,i_atom,j_atom,i_cell)=delta_r_ims(1:3,i_im)+&
      &delta_r_corr(1:3)
    !write(*,*) i_cell,i_atom,j_atom,i_im,delta_prim(1:3,i_im,i_atom,j_atom,i_cell)
    enddo ! i_im
   enddo ! i_cell
  enddo ! j_atom
 enddo ! i_atom

 open(unit=10,file='delta_prim.dat',status='replace')
 do i_cell=1,no_prim_cells
  do i_atom=1,basis
   do j_atom=1,basis
    do i_im=1,no_im_cells(i_atom,j_atom,i_cell)
     write(10,*) i_cell,i_atom,j_atom,i_im,delta_prim(1:3,i_im,i_atom,j_atom,i_cell)
    enddo ! i_im
   enddo ! j_atom
  enddo ! i_atom
 enddo ! i_cell
 close(10)


! Determine the mapping of the atoms in the primitive cell under each symmetry
! operation
  call atom_symmetry_maps(prim_latt_vecs,basis,atom_prim_cart,no_symm_ops,&
  &point_symms,trans_symms,atom_map_symm_forwards,atom_map_symm_backwards)

! Map each k-point on the grid to one in the IBZ
 call match_kpoints(prim_rec_vecs,no_grid_points,grid_points_cart,&
  &no_ibz_points,ibz_points_cart,no_symm_ops,point_symms,grid_to_ibz_map,&
  &ibz_to_grid_symm,time_reversal)

! Determine the dynamical matrix at each point on the k-point grid by applying 
! the appropriate unitary trasformation to the dynamical matrix at a point in 
! the IBZ
 do i_grid=1,no_grid_points
  i_point=grid_to_ibz_map(i_grid)
  i_symm=ibz_to_grid_symm(i_grid)
  call g_matrix_phases(ibz_points_cart(1:3,i_point),&
   point_symms(1:3,1:3,i_symm),trans_symms(1:3,i_symm),basis,&
   &atom_prim_cart(1:3,1:basis),phase(1:basis,1:basis))
  call apply_symmetry_to_dyn_mat(basis,atom_map_symm_forwards(1:basis,i_symm),&
   &phase(1:basis,1:basis),dyn_mats_ibz(1:basis,1:3,1:basis,1:3,i_point),&
   &point_symms(1:3,1:3,i_symm),dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_grid))
  if(time_reversal(i_grid))dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_grid)=&
   conjg(dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_grid))
  if(two_k_is_a_G(i_grid))dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_grid)=&
   &cmplx(real(dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_grid)),0.d0,dp)
 enddo ! i_grid

! Symmetrize the dynamical matrix at each wave vector with respect to the
! translation operations of the supercell
 dyn_mats_symm(1:basis,1:3,1:basis,1:3,1:no_grid_points)=cmplx(0.d0,0.d0,dp)
 do i_grid=1,no_grid_points
  do i_cell=1,no_grid_points
   call g_matrix_phases(grid_points_cart(1:3,i_grid),identity(1:3,1:3),&
    &cell_pos_cart(1:3,i_cell),basis,atom_prim_cart(1:3,1:basis),&
    &phase(1:basis,1:basis))
   call apply_symmetry_to_dyn_mat(basis,identity_map(1:basis),&
    &phase(1:basis,1:basis),dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_grid),&
    &identity(1:3,1:3),temp_dyn_mat(1:basis,1:3,1:basis,1:3))
   dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)=&
    &dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)+&
    &temp_dyn_mat(1:basis,1:3,1:basis,1:3)
  enddo ! i_cell
  dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)=&
   &dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)/dble(no_grid_points)
  if(two_k_is_a_G(i_grid))dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)=&
   &cmplx(real(dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)),0.d0,dp)
 enddo ! i_grid
 dyn_mats_grid(1:basis,1:3,1:basis,1:3,1:no_grid_points)=&
  &dyn_mats_symm(1:basis,1:3,1:basis,1:3,1:no_grid_points)

! Symmetrize the dynamical matrix at each wave vector with respect to its
! point and inverse groups
 dyn_mats_symm(1:basis,1:3,1:basis,1:3,1:no_grid_points)=cmplx(0.d0,0.d0,dp)
 do i_grid=1,no_grid_points
  counter=0
  do i_symm=1,no_symm_ops
   if(grid_map_symm_backwards(i_grid,i_symm)>0)then
    counter=counter+1
    i_back=grid_map_symm_backwards(i_grid,i_symm)
    call g_matrix_phases(grid_points_cart(1:3,i_back),&
     &point_symms(1:3,1:3,i_symm),trans_symms(1:3,i_symm),basis,&
     &atom_prim_cart(1:3,1:basis),phase(1:basis,1:basis))
    call apply_symmetry_to_dyn_mat(basis,&
     &atom_map_symm_forwards(1:basis,i_symm),phase(1:basis,1:basis),&
     &dyn_mats_grid(1:basis,1:3,1:basis,1:3,i_back),&
     &point_symms(1:3,1:3,i_symm),temp_dyn_mat(1:basis,1:3,1:basis,1:3))
    dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)=&
     &dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)+&
     &temp_dyn_mat(1:basis,1:3,1:basis,1:3)
   endif ! grid_map_symm_backwards
  enddo ! i_symm
  dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)=&
   &dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)/dble(counter)
  if(two_k_is_a_G(i_grid))dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)=&
   &cmplx(real(dyn_mats_symm(1:basis,1:3,1:basis,1:3,i_grid)),0.d0,dp)
 enddo ! i_grid
 dyn_mats_grid(1:basis,1:3,1:basis,1:3,1:no_grid_points)=&
  &dyn_mats_symm(1:basis,1:3,1:basis,1:3,1:no_grid_points)

! do i_grid=1,no_grid_points
!  open(unit=10,file='dyn_mat_symm.'//trim(i2s(i_grid))//'.dat',status='replace')
!  do i_atom=1,basis
!   do i_cart=1,3
!    do j_atom=1,basis
!     do j_cart=1,3
!      write(10,*)i_atom,i_cart,j_atom,j_cart,&
!       &real(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid)),&
!       &aimag(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid))
!     enddo ! j_cart
!    enddo ! j_atom
!   enddo ! i_cart
!  enddo ! i_atom
!  close(10)
! enddo ! i_grid

! Construct the matrix of force constants
 force_consts=0.d0
 do i_cell=1,no_prim_cells
  do i_atom=1,basis
   do i_cart=1,3
    do j_atom=1,basis
     do j_cart=1,3
      prefactor=sqrt(mass(i_atom)*mass(j_atom))/&
       &dble(no_grid_points)
      do i_grid=1,no_grid_points
       exp_i_k_dot_r=cmplx(0.d0,0.d0,dp)
       do i_im=1,no_im_cells(i_atom,j_atom,i_cell)
        k_dot_r=dot_product(grid_points_cart(1:3,i_grid),&
         &delta_prim(1:3,i_im,i_atom,j_atom,i_cell))
        exp_i_k_dot_r=exp_i_k_dot_r+cmplx(cos(k_dot_r),sin(k_dot_r),dp)
       enddo ! i_im
       exp_i_k_dot_r=exp_i_k_dot_r/dble(no_im_cells(i_atom,j_atom,i_cell))
       force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)=&
        &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)+&
        &real(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid)*exp_i_k_dot_r)
      enddo ! i_grid
      force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)=prefactor*&
       &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)
     enddo ! j_cart
    enddo ! j_atom
   enddo ! i_cart
  enddo ! i_atom
 enddo ! i_cell

 open(unit=10,file='force.dat',status='replace')
 do i_cell=1,no_prim_cells
  do i_atom=1,basis
   do i_cart=1,3
    do j_atom=1,basis
     do j_cart=1,3
      write(10,*) i_cell,i_atom,i_cart,j_atom,j_cart,force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)
     enddo ! j_cart
    enddo ! j_atom
   enddo ! i_cart
  enddo ! i_atom
 enddo ! i_cell
 close(10)

! Get the number of high symmetry points on the dispersion path and allocate 
! corresponding arrays
 call system("echo $(wc -l path.dat | awk '{print $1}') > tempfile.dat",istat)
 if(istat/=0)then
  write(*,*)'Problem counting the number of lines in path.dat file.'
  stop
 endif ! istat
 open(unit=10,file='tempfile.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening tempfile.dat file.'
  stop
 endif ! ierr
 read(10,*,iostat=ierr)no_kpoints_path
 if(ierr/=0)then
  write(*,*)'Problem reading tempfile.dat file.'
  stop
 endif ! ierr
 close(10,status='delete')

 allocate(path(3,no_kpoints_path),stat=ialloc)
 if(ialloc/=0)call erralloc('PATH')

 call read_path(no_kpoints_path,path)


! open(unit=101,file='temperature.dat',status='old',iostat=ierr)
! if(ierr/=0)then
!  write(*,*)'Problem opening temperature.dat file.'
!  stop
! endif ! ierr
! read(101,*,iostat=ierr)temperature
! if(ierr/=0)then
!  write(*,*)'Problem reading temperature.dat file.'
!  stop
! endif ! ierr
! close(101)

 call generate_dispersion(prim_rec_vecs,basis,mass,no_grid_points,&
  &no_im_cells,delta_prim,force_consts,no_kpoints_path,path)

 call generate_dos(prim_rec_vecs,basis,mass,no_grid_points,no_im_cells,&
  &delta_prim,force_consts,temperature)

 call generate_disp_patterns(no_prim_cells,cell_pos_cart,no_grid_points,&
  &grid_points_cart,reference,basis,mass,dyn_mats_grid)

 call eval_freqs_on_grid(basis,no_grid_points,&
  &dyn_mats_grid(1:basis,1:3,1:basis,1:3,1:no_grid_points),temperature)

 END PROGRAM fourier_interpolation
