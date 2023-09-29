!! This file is part of BSA Library.
!! Copyright (C) 2023  Michele Esposito Marzino 
!!
!! BSA Library is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BSA Library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BSA Library.  If not, see <https://www.gnu.org/licenses/>.
submodule(Logging) LoggingImpl

   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG, NOTEMSG, int32
   use BsaLib_Data, only: bsa_Abort
   implicit none
   character(len = 256) :: fmt
   character(len=0),  parameter :: lev0 = ''
   character(len=3),  parameter :: lev1 = '   '


contains


   module subroutine init(this, iun, fname)
      class(logger_t),  intent(inout) :: this
      integer(int32),   intent(in)    :: iun
      character(len=*), intent(in), optional :: fname

      this%iun_ = iun
      if (present(fname)) this%fileName_ = fname
   end subroutine



   module function name(this) result(nam)
      class(logger_t), intent(in)   :: this
      character(len=:), allocatable :: nam

      nam = this%fileName_
   end function

   module subroutine setName(this, fname)
      class(logger_t), intent(inout) :: this
      character(len=*), intent(in)   :: fname

      this%fileName_ = fname
   end subroutine setName



   module function unit(this) result(iun)
      class(logger_t), intent(in) :: this
      integer :: iun

      iun = this%iun_
   end function


   module subroutine setUnit(this, iun)
      class(logger_t), intent(inout) :: this
      integer, intent(in), target    :: iun

      call this%init(iun)
   end subroutine
   



   module subroutine logZonePremeshingTotTime(this, zname, rtime, npts, print2console)
      class(logger_t), intent(in)  :: this
      character(len=*), intent(in) :: zname
      real(real64), intent(in) :: rtime
      integer, intent(in), optional :: npts
      logical, intent(in), optional :: print2console

      logical :: do_print = .false.

      if (present(print2console) .and. print2console) do_print = .true.

      write(unit=fmt, fmt='(a, "Done zone  ", a)') INFOMSG, zname
      write(this%iun_, '(1x, a)') fmt(1 : len_trim(fmt))
      if (do_print) print '(1x, a)', fmt(1 : len_trim(fmt))

      if (present(npts)) then

         write(unit=fmt, fmt='( a, "- n. of zone points:  ", i0 )') &
            MSGCONT, npts
         write(this%iun_, '(1x, a)') fmt(1 : len_trim(fmt))
         if (do_print) print '(1x, a)', fmt(1 : len_trim(fmt))
      endif

      write(unit=fmt, fmt='(a,     "- TOT. TIME:          ", g0, " [s]")') &
         MSGCONT, rtime
      write(this%iun_, '(1x, a, /)') fmt(1 : len_trim(fmt))
      if (do_print) print '(1x, a, /)', fmt(1 : len_trim(fmt))
   end subroutine logZonePremeshingTotTime














!*****************************************************************************************
!     ALLOCATION
!*****************************************************************************************


! #ifdef _BSA_ALLOC_DEBUG
!    module subroutine allocOKMsg_scalar_(name_, iloc, nbytes)
!       character(len = *), intent(in) :: name_
!       integer(kind = 8), intent(in), optional  :: iloc, nbytes

!       write(unit_debug_, fmt='(3a)', advance='no') &
!          'variable  "', name_, '"  allocated.'

!       if (present(iloc)) then
!          write(unit_debug_, '(a, i0)', advance='no') &
!             'Location in memory:  ', iloc
!       endif

!       if (present(nbytes)) then
!          write(unit_debug_, fmt='(a, i0, ".")', advance='no') &
!             'Occupancy (bytes):  ', nbytes
!          write(unit_debug_, *) ''
!       endif
!    end subroutine


!    module subroutine allocOKMsg_array_(name_, dims, iloc, nbytes)
!       character(len = *), intent(in) :: name_
!       integer, intent(in)            :: dims(..)
!       integer(kind = 8), intent(in), optional  :: iloc, nbytes
!       integer :: dim, ndims

!       write(unit_debug_, fmt='(3a)', advance='no') &
!          'variable  "', name_, '"  allocated.'

!       select rank (dims)
!          rank (0)
!             dim = dims
!             write(unit_debug_, '(a, i0)') &
!                'Dimension:  ', dim
!          rank (1)
!             ndims = size(dims)
!             write(unit_debug_, '(a, i0, *(" - ", i0) )') &
!                'Dimension:  ', (dims(dim), dim = 1, ndims)

!          rank default

!             print '(1x, a, a)', &
!                ERRMSG, &  
!                ' Dimensions for a NDrank array allocation must be at most 1D-rank array.'
!             call bsa_Abort()
!       end select

!       if (present(iloc)) then
!          write(unit_debug_, '(a, i0)', advance='no') &
!             'Location in memory:  ', iloc
!       endif

!       if (present(nbytes)) then
!          write(unit_debug_, fmt='(a, i0, ".")', advance='no') &
!             'Occupancy (bytes):  ', nbytes
!          write(unit_debug_, *) ''
!       endif
!    end subroutine


!    module subroutine deallocOKMsg(name_)
!       character(len = *), intent(in) :: name_

!       write(unit_debug_, fmt='(3a)') &
!          'variable  "', name_, '"  de-allocated.'
!    end subroutine
! #endif



   module subroutine allocKOMsg(name_, istat, emsg)
      character(len = *), intent(in) :: name_, emsg
      integer, intent(in) :: istat

      write(unit_debug_, fmt='(3a)') &
         '[ERROR] variable  "', name_, '"  could not be allocated.'
      write(unit_debug_, fmt='(15x, a, i0, 2a)') &
         'Exit code  ', istat, '. Error message:  ', emsg(1 : len_trim(emsg))

      call bsa_Abort()
   end subroutine


   module subroutine deallocKOMsg(name_, istat, emsg)
      character(len = *), intent(in) :: name_, emsg
      integer, intent(in) :: istat

      write(unit_debug_, fmt='(3a)') &
         '[ERROR] variable  "', name_, '"  could not be de-allocated.'
      write(unit_debug_, fmt='(15x, a, i0, 2a)') &
         'Exit code  ', istat, '. Error message:  ', emsg(1 : len_trim(emsg))

      call bsa_Abort()
   end subroutine



end submodule


























! subroutine LogWindSpectralAnalysisHeader(logger, nz, ivaru, isu)
!    ! in
!    class(logger_t), intent(in) :: logger
!    integer, intent(in) :: nz, ivaru
!    integer, intent(in), optional :: isu

!    10 format( //,&
!       ' ANALYSE SPECTRALE VENT',/, &
!       ' **********************',/, &
!       ' NOMBRE DE ZONES    :  ',I4,/, &
!       ' TYPE DE VARIATION  DU VENT AVEC L ALTITUDE : ',A10,/ &
!       ' DENSITE SPECTRALE DE VENT : VON KARMAN' / )

!    11 format( //,&
!       ' ANALYSE SPECTRALE VENT',/, &
!       ' **********************',/, &
!       ' NOMBRE DE ZONES    :  ',I4,/, &
!       ' TYPE DE VARIATION  DU VENT AVEC L ALTITUDE : ',A12,/ &
!       ' DENSITE SPECTRALE DE VENT : ', A12 / )


!    if (present(isu)) then
!       write(logger%iun_, fmt=11), nz, CST_WIND_V_PROFILES(ivaru), CST_PSD_TYPES(isu)
!    else
!       write(logger%iun_, fmt=10), nz, CST_WIND_V_PROFILES(ivaru)
!    endif
! end subroutine LogWindSpectralAnalysisHeader






! subroutine LogWindSPeedProfileType(logger, ivaru)
!    class(logger_t), intent(in) :: logger
!    integer, intent(in) :: ivaru

!    10 format( & 
!       '  !!!!!!!!!!!LOI SPECIALE MILLAU ZREF= 260.07  !!!!!!!!!!', /& 
!       '  -------------------------------------------------------' )
!    11 format( & 
!       '  !!!!!!!!!!!LOI SPECIALE MILLAU MAQU ZREF= 0.8669  !!!!!', /& 
!       '  ---------------------------------------------------------' )
!    if (ivaru == 3) then
!       write(logger%iun_, fmt=10)
!    elseif(ivaru == 4) then
!       write(logger%iun_, fmt=11)
!    endif
! end subroutine LogWindSPeedProfileType





! subroutine LogWindZoneData(logger, iz, xlim, ub, alpha, al, ect, corr, expp)
!    class(logger_t), intent(in) :: logger
!    integer, intent(in) :: iz
!    real, intent(in)    :: xlim( : ), ub, alpha, ect( 3 )
!    real, intent(in), dimension(3,3) :: al, corr, expp 
!    ! local
!    integer i, j

!    10 format( /,&
!       ' ZONE ',I2,' DE X = ',F10.2,' A X= ',F10.2,/, & 
!       ' -----------------------------------------',/, & 
!       ' VITESSE DE BASE   : ', F10.3,/, & 
!       ' ALPHA / Z0        : ', F10.3,/, & 
!       ' I        LXI       LYI       LZI EC.TYPE I',/, &
!       ' U ',4F10.5 ,/ & 
!       ' V ',4F10.5 ,/ & 
!       ' W ',4F10.5 ,/ & 
!       ' I   CORR.X I  CORR.Y I  CORR.Z I  EXPP.X I  EXPP.Y I  EXPP.Z I', /, &
!       ' U ',3(1x,g9.4),3F10.5 ,/&
!       ' V ',3(1x,g9.4),3F10.5 ,/&
!       ' W ',3(1x,g9.4),3F10.5 &
!    )

!    write(logger%iun_, 10), iz, xlim( iz ), xlim( iz + 1 ), ub, alpha, &
!                      ( &
!                         ( al( i,j ), j = 1,3 ), &
!                           ect( i ),&
!                      i = 1, 3 ), &
!                      ( &
!                         ( corr( i,j ), j = 1,3 ), &
!                         ( expp( i,j ), j = 1,3 ), &
!                      i = 1, 3 )
! end subroutine LogWindZoneData







! subroutine LogTransFuncType(logger, itran)
!    class(logger_t), intent(in) :: logger
!    integer, intent(in) :: itran

!    10 format ( '  ITRAN = ', 1I1, & 
!                ' : TRANSFER MATRIX ASSUMED DIAGONAL (COUPLING COMES FROM THE ONE EXISTENT BETWEEN MODAL FORCES.)'/ )
!    11 format ( '  ITRAN = ', 1I1, & 
!                ' : TRANSFER MATRIX COMPUTED ACCOUNTING FOR COUPLING (COMING FROM NON-PROPORTIONAL DAMPING.)'/ )
!    12 format ( '  ITRAN = ', 1I1, & 
!                ' : TRANSFER MATRIX COMPUTED ACCOUNTING FOR COUPLING (SYMPLIFIED METHOD.)'/ )

!    ! stocap ( french translation )
!    13 format(  '  ITRAN = ', 1I1, & 
!                ' : MATRICE DE TRANSFERT SUPPOSEE DIAGONALE'/ & 
!                '  (LE COUPLAGE PROVIENT DE CELUI EXISTANT ENTRE F.GEN.)'/ )
!    14 format(  '  ITRAN = ', 1I1, & 
!                ' : MATRICE DE TRANSFERT CALCULEE AVEC COUPLAGE'/ & 
!                '  (PROVENANT DE L''AMORTISSEMENT NON PROPORTIONNEL)'/ )
!    15 format(  '  ITRAN = ', 1I1, & 
!                ' : MATRICE DE TRANSFERT CALCULEE AVEC COUPLAGE'/ & 
!                '  (METHODE SIMPLIFIEE)'/ )

   
!    if ( logger%iun_ == 6 ) then

!       if (itran==3) then
!          write(logger%iun_, 13), itran
!       elseif (itran==2) then
!          write(logger%iun_, 14), itran
!       elseif ( itran==4) then
!          write(logger%iun_, 15), itran
!       endif

!    else

!       if (itran==3) then
!          write(logger%iun_, 10), itran
!       elseif (itran==2) then
!          write(logger%iun_, 11), itran
!       elseif ( itran==4) then
!          write(logger%iun_, 12), itran
!       endif
!    endif
! end subroutine LogTransFuncType








! subroutine LogMaxValueType(logger, ityp, T)
!    class(logger_t), intent(in) :: logger
!    integer, intent(in) :: ityp
!    real, intent(in) :: T
!    ! local
!    character(len = 40) :: typ_max

!    10 format ( /'  TYPE DU MAXIMUM   : ', A40 )
!    11 format ( '  DUREE OBSERVATION : ', 1F8.3 )

!    if (ityp==0) then
!       typ_max='GAUSS : G = 3'
!    elseif (ityp==1) then
!       typ_max='POISSON - CALCUL DE LA LARGEUR DE BANDE'
!    elseif (ityp==2) then
!       typ_max='VANMARCKE - POISSON MODIFIE'
!    elseif (ityp==3) then
!       typ_max='POISSON - PROCESSUS EN BANDE ETROITE'
!    elseif (ityp==4) then
!       typ_max='G - IMPOSE'
!    endif

!    write( logger%iun_, 10 ) typ_max
!    IF (ityp==1.OR.ityp==2) write( logger%iun_, 11 ), T
! end subroutine LogMaxValueType





! subroutine LogElementWindLoad(logger, w, interv, y_c, z_c)
!    class(logger_t), intent(in) :: logger
!    real, intent(in)    :: w 			! element's width
!    real, intent(in)    :: interv		! delta_I ( incidence angle ) for derivative computation
!    real, intent(in)    :: y_c, z_c	! Y and Z coordinates of Center (cross-section) point

!    10 format( /, 25X,&
!       ' CHARGE DE VENT SUR L ELEMENT ',/, 25x, & 
!       ' ---------------------------- ',/, 25x, & 
!       ' MAITRE COUPLE :                  ',F10.4,/, 25x, & 
!       ' INTERVALLE POUR CALCUL DERIVEE : ',F10.4,/, 25x, & 
!       ' POSITON CENTRE M.C. :      YFV = ', F10.4, ' ZFV = ', F10.4 )

!    write(logger%iun_, 10), w, interv, y_c, z_c
! end subroutine LogElementWindLoad





! subroutine LOgDLMWindCoeffs(logger, n_coefs, ai, cx, cz, cm)
!    class(logger_t), intent(in) :: logger
!    integer, intent(in)  :: n_coefs              ! num of coeffs ( for given incidence angle )
!    real, intent(in), dimension(:) :: ai         ! angles of incidence ( given ) (  ERS )
!    real, intent(in), dimension(:) :: cx, cz, cm ! wind DRAG/LIFT/MOMENT coeffs (in GWRS )
!    ! local
!    integer i

!    10 format( &
!       25x, '    COEFFICIENTS DE TRAINEE ET PORTANCE',/, &
!       25X, '        I         Y         Z         M ', /, &
!       (25x, 4F10.4) )

!    write(logger%iun_, 10), (ai(i), cx(i), cz(i), cm(i), i = 1, n_coefs)
! end subroutine LOgDLMWindCoeffs







! subroutine LogElementWindIncidenceAngles(logger, ai0, ain1, ain2)
!    class(logger_t), intent(in) :: logger
!    real, intent(in)    :: ai0 	! mean incidence angle
!    real, intent(in)    :: ain1	! ??
!    real, intent(in)    :: ain2	! ??

!    10 format( 25X, ' ANGLES D INCIDENCE : I0= ', g12.5, ' I1= ', g12.5, ' I2= ', g12.5 )
!    write(logger%iun_, 10), ai0, ain1, ain2
! end subroutine LogElementWindIncidenceAngles







! subroutine LogElemWindNodalVel(logger, u1, u2)
!    implicit none
!    class(logger_t), intent(in) :: logger
!    real, intent(in)    :: u1, u2

!    10 format( 25X,' VITESSES :           U1= ',F10.2,'   U2= ',F10.2 )
!    write(logger%iun_, 10), u1, u2
! end subroutine LogElemWindNodalVel