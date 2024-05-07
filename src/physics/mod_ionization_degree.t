!> module ionization degree - get ionization degree for given temperature
module mod_ionization_degree
  implicit none
  double precision, private :: Te_H_CL(1:100)
  double precision, private :: iz_H_CL(1:100)
  double precision, dimension(:), allocatable :: Te_H_table
  double precision, dimension(:), allocatable :: iz_H_table
  double precision :: Te_table_min, Te_table_max, Te_table_step, inv_Te_table_step
  ! Carlsson, M., & Leenaarts, J. 2012, A&A, 539, A39
  data Te_H_CL /2.1000005d+03, 2.2349495d+03, 2.3785708d+03, 2.5314199d+03, 2.6940928d+03, &
                2.8672190d+03, 3.0514692d+03, 3.2475613d+03, 3.4562544d+03, 3.6783584d+03, &
                3.9147351d+03, 4.1662998d+03, 4.4340322d+03, 4.7189697d+03, 5.0222153d+03, &
                5.3449502d+03, 5.6884248d+03, 6.0539712d+03, 6.4430083d+03, 6.8570420d+03, &
                7.2976860d+03, 7.7666460d+03, 8.2657383d+03, 8.7969062d+03, 9.3622090d+03, &
                9.9638330d+03, 1.0604130d+04, 1.1285561d+04, 1.2010781d+04, 1.2782619d+04, &
                1.3604041d+04, 1.4478250d+04, 1.5408652d+04, 1.6398826d+04, 1.7452648d+04, &
                1.8574172d+04, 1.9767766d+04, 2.1038082d+04, 2.2390010d+04, 2.3828838d+04, &
                2.5360100d+04, 2.6989764d+04, 2.8724182d+04, 3.0570023d+04, 3.2534518d+04, &
                3.4625215d+04, 3.6850262d+04, 3.9218293d+04, 4.1738543d+04, 4.4420699d+04, &
                4.7275266d+04, 5.0313219d+04, 5.3546391d+04, 5.6987395d+04, 6.0649453d+04, &
                6.4546914d+04, 6.8694758d+04, 7.3109148d+04, 7.7807289d+04, 8.2807258d+04, &
                8.8128625d+04, 9.3791852d+04, 9.9819000d+04, 1.0623346d+05, 1.1306024d+05, &
                1.2032560d+05, 1.2805797d+05, 1.3628709d+05, 1.4504503d+05, 1.5436592d+05, &
                1.6428561d+05, 1.7484295d+05, 1.8607852d+05, 1.9803609d+05, 2.1076208d+05, &
                2.2430608d+05, 2.3872045d+05, 2.5406084d+05, 2.7038703d+05, 2.8776234d+05, &
                3.0625456d+05, 3.2593475d+05, 3.4688000d+05, 3.6917081d+05, 3.9289406d+05, &
                4.1814181d+05, 4.4501247d+05, 4.7360988d+05, 5.0404450d+05, 5.3643488d+05, &
                5.7090725d+05, 6.0759431d+05, 6.4663956d+05, 6.8819319d+05, 7.3241712d+05, &
                7.7948294d+05, 8.2957412d+05, 8.8288331d+05, 9.3961919d+05, 1.0000000d+06  /

  ! Carlsson, M., & Leenaarts, J. 2012, A&A, 539, A39
  data iz_H_CL /1.0000000d+00, 1.0000000d+00, 1.0000000d+00, 1.0000000d+00, 1.0000000d+00, &
                1.0000000d+00, 1.0000000d+00, 9.9999928d-01, 9.9999774d-01, 9.9999624d-01, &
                9.9999428d-01, 9.9999219d-01, 9.9998862d-01, 9.9998313d-01, 9.9997944d-01, &
                9.9995774d-01, 1.0000324d+00, 9.9822283d-01, 9.6391600d-01, 8.8995951d-01, &
                8.0025935d-01, 7.4219799d-01, 7.1378660d-01, 6.9565988d-01, 6.7326778d-01, &
                6.5032905d-01, 6.2487972d-01, 5.9925246d-01, 5.7028764d-01, 5.3673893d-01, &
                4.9144468d-01, 4.4332290d-01, 3.9158964d-01, 3.4621361d-01, 2.7826238d-01, &
                1.9633104d-01, 1.2593637d-01, 8.1383914d-02, 5.0944358d-02, 2.8770180d-02, &
                1.4449068d-02, 7.6071853d-03, 4.3617180d-03, 2.7694483d-03, 1.8341533d-03, &
                1.2467797d-03, 8.6159195d-04, 6.0232455d-04, 4.2862241d-04, 3.0666622d-04, &
                2.2256430d-04, 1.6328457d-04, 1.2179965d-04, 9.2065704d-05, 7.0357783d-05, &
                5.4168402d-05, 4.2026968d-05, 3.2807184d-05, 2.5801779d-05, 2.0414085d-05, &
                1.6271379d-05, 1.3041737d-05, 1.0538992d-05, 8.5329248d-06, 6.9645471d-06, &
                5.7069005d-06, 4.7155136d-06, 3.9203474d-06, 3.2734958d-06, 2.7415674d-06, &
                2.2911979d-06, 1.9045801d-06, 1.5692771d-06, 1.2782705d-06, 1.0281080d-06, &
                8.1567373d-07, 6.4200310d-07, 5.0668103d-07, 4.0577518d-07, 3.3245340d-07, &
                2.8006056d-07, 2.4119515d-07, 2.1033340d-07, 1.8429903d-07, 1.6189058d-07, &
                1.4264889d-07, 1.2616209d-07, 1.1193421d-07, 9.9408616d-08, 8.8335305d-08, &
                7.8593857d-08, 7.0042269d-08, 6.2532941d-08, 5.5898361d-08, 5.0049657d-08, &
                4.4888207d-08, 4.0317172d-08, 3.6192560d-08, 3.2330732d-08, 2.9387850d-08  /
    double precision :: Te_low_iz_He=5413.d0
  contains

    subroutine ionization_degree_init()
      use mod_global_parameters
  
      double precision :: fact1, fact2, fact3, dL1, dL2
      integer :: i, j, ntable, n_interpolated_table
      logical :: jump

      ntable=100
      n_interpolated_table=4000

      ! change neutral fraction to ionization degree
      iz_H_CL=1.d0-iz_H_CL
      
      allocate(Te_H_table(1:n_interpolated_table))
      allocate(iz_H_table(1:n_interpolated_table))
      Te_H_table=zero
      iz_H_table=zero
    
      ! create cooling table(s) for use in amrvac
      Te_table_step = (Te_H_CL(ntable)-Te_H_CL(1))/dble(n_interpolated_table-1)

      Te_H_table(1) = Te_H_CL(1)
      iz_H_table(1) = iz_H_CL(1)
      
      Te_H_table(n_interpolated_table) = Te_H_CL(ntable)
      iz_H_table(n_interpolated_table) = iz_H_CL(ntable)
      
      do i=2,n_interpolated_table        ! loop to create one table
        Te_H_table(i) = Te_H_table(i-1)+Te_table_step
        do j=1,ntable-1   ! loop to create one spot on a table
        ! Second order polynomial interpolation, except at the outer edge, 
        ! or in case of a large jump.
          if(Te_H_table(i) < Te_H_CL(j+1)) then
             if(j.eq. ntable-1 )then
               fact1 = (Te_H_table(i)-Te_H_CL(j+1))     &
                     /(Te_H_CL(j)-Te_H_CL(j+1)) 
      
               fact2 = (Te_H_table(i)-Te_H_CL(j))       &
                     /(Te_H_CL(j+1)-Te_H_CL(j)) 
      
               iz_H_table(i) = iz_H_CL(j)*fact1 + iz_H_CL(j+1)*fact2 
               exit
             else 
               dL1 = iz_H_CL(j+1)-iz_H_CL(j)
               dL2 = iz_H_CL(j+2)-iz_H_CL(j+1)
               jump =(max(dabs(dL1),dabs(dL2)) > 2.d0*min(dabs(dL1),dabs(dL2)))
             end if
               
             if( jump ) then
               fact1 = (Te_H_table(i)-Te_H_CL(j+1))     &
                     /(Te_H_CL(j)-Te_H_CL(j+1)) 
      
               fact2 = (Te_H_table(i)-Te_H_CL(j))       &
                     /(Te_H_CL(j+1)-Te_H_CL(j)) 
                      
               iz_H_table(i) = iz_H_CL(j)*fact1 + iz_H_CL(j+1)*fact2
               exit          
             else
               fact1 = ((Te_H_table(i)-Te_H_CL(j+1))     &
                     * (Te_H_table(i)-Te_H_CL(j+2)))   &
                     / ((Te_H_CL(j)-Te_H_CL(j+1)) &
                     * (Te_H_CL(j)-Te_H_CL(j+2)))
    
               fact2 = ((Te_H_table(i)-Te_H_CL(j))       &
                     * (Te_H_table(i)-Te_H_CL(j+2)))   &
                     / ((Te_H_CL(j+1)-Te_H_CL(j)) &
                     * (Te_H_CL(j+1)-Te_H_CL(j+2)))
    
               fact3 = ((Te_H_table(i)-Te_H_CL(j))       &
                     * (Te_H_table(i)-Te_H_CL(j+1)))   &
                     / ((Te_H_CL(j+2)-Te_H_CL(j)) &
                     * (Te_H_CL(j+2)-Te_H_CL(j+1)))
      
               iz_H_table(i) = iz_H_CL(j)*fact1 + iz_H_CL(j+1)*fact2 &
                        + iz_H_CL(j+2)*fact3
               exit
             endif
          endif
        enddo  ! end loop to find create one spot on a table
      enddo    ! end loop to create one table

      ! nondimensionalize temperature
      Te_H_table=Te_H_table/unit_temperature
      Te_table_min=Te_H_table(1)
      Te_table_max=Te_H_table(n_interpolated_table)
      inv_Te_table_step=unit_temperature/Te_table_step

      ! transition temperature
      Te_low_iz_He=Te_low_iz_He/unit_temperature
  
    end subroutine ionization_degree_init

    subroutine ionization_degree_from_temperature(ixI^L,ixO^L,Te,iz_H,iz_He)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: Te(ixI^S)
      double precision, intent(out) :: iz_H(ixO^S),iz_He(ixO^S)

      integer :: ix^D, i

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(Te(ix^D)<Te_table_min) then
          iz_H(ix^D)=0.d0
        else if (Te(ix^D)>Te_table_max) then
          iz_H(ix^D)=1.d0
        else
          i=int((Te(ix^D)-Te_table_min)*inv_Te_table_step)+1
          iz_H(ix^D)=iz_H_table(i)+(Te(ix^D)-Te_H_table(i))&
           *(iz_H_table(i+1)-iz_H_table(i))*inv_Te_table_step
        end if
        ! Ni, L. et al. A&A, 665, A116
        if(Te(ix^D)<Te_low_iz_He) then
          iz_He(ix^D)=1.0084814d-4
        else
          iz_He(ix^D)=1.d0-10**(0.322571d0-5.96d-5*Te(ix^D)*unit_temperature)
        end if
      {end do\}

    end subroutine ionization_degree_from_temperature

    ! gas constant R in ideal gas law for solar plasma
    function R_ideal_gas_law_partial_ionization(Te)
      use mod_global_parameters
      double precision, intent(in) :: Te
      double precision :: R_ideal_gas_law_partial_ionization

      double precision :: iz_H, iz_He
      integer :: i

      if(Te<Te_table_min) then
        iz_H=0.d0
      else if (Te>Te_table_max) then
        iz_H=1.d0
      else
        i=int((Te-Te_table_min)*inv_Te_table_step)+1
        iz_H=iz_H_table(i)+(Te-Te_H_table(i))&
         *(iz_H_table(i+1)-iz_H_table(i))*inv_Te_table_step
      end if
      if(Te<Te_low_iz_He) then
        iz_He=1.0084814d-4
      else
        iz_He=1.d0-10**(0.322571d0-5.96d-5*Te*unit_temperature)
      end if
      ! dimensionless: kB and mp are included in units 
      ! assume the first and second ionization of Helium have the same degree
      R_ideal_gas_law_partial_ionization=(1.d0+iz_H+0.1d0*(1.d0+iz_He*(1.d0+iz_He)))/2.3d0
      return

    end function R_ideal_gas_law_partial_ionization

end module mod_ionization_degree
