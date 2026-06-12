c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2026  by  Yanxing Wang        ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine readexfldatm  --  read atom-specific ext field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "readexfldatm" reads atom-specific external electric field
c     values from a file and stores them in the exfld_atm array for 3D-RISM
c
c     filename format: each line contains atom_id and three field components in V/nm
c     example: 1/Acetyl_Cap.C  2.96e-05 -6.56e-05  3.66e-05
c
c
      subroutine readexfldatm (filename)
      use atoms
      use extfld
      use iounit
      use units
      implicit none
      integer i,iatm,next
      integer iefatm,islash
      integer freeunit
      real*8 fx,fy,fz,fxx, fxy, fxz, fyy, fyz, fzz
      character*40 atmstr
      character*40 dummy
      character*240 filename
      character*240 record
      
c
c     open the external field file
c
      write (iout,10)  trim(filename)
   10 format (/,' Reading Atom-Specific External Field File: ',a)

      iefatm = freeunit ()
      open (unit=iefatm,file=filename,status='old')
      rewind (unit=iefatm)
c
c     read each line from the file
c
      do while (.true.)
         read (iefatm,20,err=40,end=40)  record
   20    format (a240)
c        
c        skip blank lines
c
         if (len_trim(record) .eq. 0)  goto 30
c
c        extract atom id number from first field (before the slash)
c
         islash = index(record,'/')
         if (islash .gt. 0) then
            atmstr = record(1:islash-1)
         else
            write (iout,25)  trim(record)
   25       format (/,' READEXFLDATM  --  Invalid Format,',
     &              ' Missing "/" in: ',/,a)
            call fatal
         end if
         read (atmstr,*,err=40,end=40)  iatm
c
c        skip the atom name after the slash and read the field values
c
         next = islash + 1
         call gettext (record,dummy,next)
         read (record(next:240),*,err=40,end=40) fx, fy, fz,
     &            fxx, fxy, fxz, fyy, fyz, fzz 
       
c        store the field components and gradients
c
         exfld_atm(1,iatm) = fx
         exfld_atm(2,iatm) = fy
         exfld_atm(3,iatm) = fz
         exfld_atm(4,iatm) = fxx
         exfld_atm(5,iatm) = fxy
         exfld_atm(6,iatm) = fxz
         exfld_atm(7,iatm) = fyy
         exfld_atm(8,iatm) = fyz
         exfld_atm(9,iatm) = fzz

c
c        convert external field from V/nm to atomic units
c
         exfld_atm(1,iatm) = exfld_atm(1,iatm) / (10 * elefield)  
         exfld_atm(2,iatm) = exfld_atm(2,iatm) / (10 * elefield)
         exfld_atm(3,iatm) = exfld_atm(3,iatm) / (10 * elefield)

c
c        convert external field gradients from V/nm^2 to atomic units
c

         exfld_atm(4,iatm) = exfld_atm(4,iatm) / (100 * elefield)  
         exfld_atm(5,iatm) = exfld_atm(5,iatm) / (100 * elefield)
         exfld_atm(6,iatm) = exfld_atm(6,iatm) / (100 * elefield)
         exfld_atm(7,iatm) = exfld_atm(7,iatm) / (100 * elefield)
         exfld_atm(8,iatm) = exfld_atm(8,iatm) / (100 * elefield)
         exfld_atm(9,iatm) = exfld_atm(9,iatm) / (100 * elefield)
   30    continue
      end do
   40 continue
c
c     close the external field file
c
      close (unit=iefatm)
c
c     print summary of external fields for all atoms
c
      write (iout,50)
   50 format (/,' External Field Summary for All Atoms:',/)
      do i = 1, n
         fx = exfld_atm(1,i) * (10 * elefield)   
         fy = exfld_atm(2,i) * (10 * elefield)
         fz = exfld_atm(3,i) * (10 * elefield)
         write (iout,60)  i, fx, fy, fz,
     &                    exfld_atm(1,i),
     &                    exfld_atm(2,i),
     &                    exfld_atm(3,i)
   60    format (' Atom',i6,' : Input=',3f12.8,
     &           ' (V/nm)  --> ',3f12.8,' (a.u.)')
      end do

c
c     print summary of external field gradients for all atoms
c

      write (iout,70)
   70 format (/,' External Field Gradients Summary for All Atoms:',/)

      do i = 1, n
         fxx = exfld_atm(4,i) * (100 * elefield) 
         fxy = exfld_atm(5,i) * (100 * elefield) 
         fxz = exfld_atm(6,i) * (100 * elefield) 
         fyy = exfld_atm(7,i) * (100 * elefield) 
         fyz = exfld_atm(8,i) * (100 * elefield) 
         fzz = exfld_atm(9,i) * (100 * elefield) 
         write (iout,80)  i, fxx, fxy, fxz,
     &                    fyy, fyz, fzz,
     &                    exfld_atm(4,i),
     &                    exfld_atm(5,i),
     &                    exfld_atm(6,i),
     &                    exfld_atm(7,i),
     &                    exfld_atm(8,i),
     &                    exfld_atm(8,i)
   80    format (' Atom',i6,' : Input=',6f12.8,
     &           ' (V/nm^2)  --> ',6f12.8,' (a.u.)')
      end do

      return
      end
