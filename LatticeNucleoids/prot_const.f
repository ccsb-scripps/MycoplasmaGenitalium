C prot_const.f
C 2021 David S. Goodsell, the Scripps Research Institute
C
C Licensed under the Apache License, Version 2.0 (the "License");
C you may not use this file except in compliance with the License.
C You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
C Unless required by applicable law or agreed to in writing, software
C distributed under the License is distributed on an "AS IS" BASIS,
C WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
C See the License for the specific language governing permissions and
C limitations under the License
C
C This work was supported in part by the US National Institutes of Health R01-GM120604
C
C Reads from unit 5:
C MG_149_145_1apr21.csv   # csv file with genomic and proteomic information (see file for format)
C MG149.segments          # input segment file
C proteins.dat            # input data file with radii and labels
C MG149.protconst         # output protein constraint file
C
C proteins.dat:
C DNA_GYRA   50.00 GYR    # label in csv file, radius (Angstroms), label
C
C Writes a file to unit 6 with positions of plectonemes:
C CONSTR   57793   57793   0.000   1.471 1 GYR  # constraint: two bead numbers,
C                                                 radius of sphere(s) (lattice units), 1 or 2 spheres, label
C                                                 --this will put a single sphere at position 57793
C CONSTR     206     250   1.324   0.662 2 SMC  # this will constrain beads 206 and 250 to be 1.324 units apart
C                                                and add two spheres to the model with radii 0.662 units
C
	character*3 name(100)
	real*4 rad(100)
	character*8 label,protlabel(100)
	character*80 line

C DNA length hardwired for mycoplasma model
	ncircle=580000/10

C protein CSV file
	read(5,*) line
	open(4,file=line)
	write(6,*) "input CSV file:  ",line

C segment file
	read(5,*) line
	open(2,file=line)
	write(6,*) "input segment file:  ",line

	ifork1=0
	ifork2=0

 300	read(2,301,end=309) line
 301	format(a80)
 302	format(4x,2i8)
	goto 300
 309	continue

C protein name file
	read(5,*) line
	open(1,file=line)
	write(6,*) "protein name file:  ",line

C output constraint file
	read(5,*) line
	open(3,file=line)
	write(6,*) "output file:  ",line

	nprot=1
 100	read(1,200, end=900) protlabel(nprot),rad(nprot),
     &   name(nprot)
 200	format(a8,f8.2,1x,a3)
	rad(nprot)=rad(nprot)/34.
	nprot=nprot+1
	goto 100
 900	continue
	nprot=nprot-1

C read in lines from the .csv file
 1	read(4,9,end=99) line
 9	format(a80)

	do iprot=1,nprot

	if (line(1:8).eq.protlabel(iprot)) then

 2	read(line,*,end=99) label,ipos,idir,imRNA

C this line ignores RNA polymerase that is transcribing
	if ((line(1:8).eq."RNA_POLY").and.(imRNA.ne.-1)) goto 230

	ipos=ipos/10
	   
	if (name(iprot).eq."SMC") then
         write(3,50) "CONSTR",ipos+1-22,ipos+1+22,
     &      rad(iprot), rad(iprot)/2., 2, name(iprot)
	else
	write(3,50) "CONSTR",ipos,ipos,
     &      0.0, rad(iprot), 1, name(iprot)
	endif

c constraints to hold the replisome together
	if (name(iprot).eq."DPL") then
	 iforkrep=ifork1
	 if ((idir.eq.1).or.(idir.eq.3)) iforkrep=ifork2
         write(3,50) "CONSTR",iforkrep,ipos+1,
     &      5.3, 0.1, 1, "REP"
	endif

C dnaA multimers
	if (name(iprot).eq."DAD") 
     &   write(3,50) "CONSTR",ipos+1,ipos+1,
     &      0.0, rad(iprot), 1, name(iprot)

	if (name(iprot).eq."DAT") then
         write(3,50) "CONSTR",ipos+1,ipos+1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-1,ipos-1,
     &      0.0, rad(iprot), 1, name(iprot)
	endif

	if (name(iprot).eq."DA4") then
         write(3,50) "CONSTR",ipos+1,ipos+1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-1,ipos-1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos+2,ipos+2,
     &      0.0, rad(iprot), 1, name(iprot)
	endif

	if (name(iprot).eq."DA5") then
         write(3,50) "CONSTR",ipos+1,ipos+1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-1,ipos-1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos+2,ipos+2,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-2,ipos-2,
     &      0.0, rad(iprot), 1, name(iprot)
	endif

	if (name(iprot).eq."DA6") then
         write(3,50) "CONSTR",ipos+1,ipos+1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-1,ipos-1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos+2,ipos+2,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-2,ipos-2,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos+3,ipos+3,
     &      0.0, rad(iprot), 1, name(iprot)
	endif

	if (name(iprot).eq."DA7") then
         write(3,50) "CONSTR",ipos+1,ipos+1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-1,ipos-1,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos+2,ipos+2,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-2,ipos-2,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos+3,ipos+3,
     &      0.0, rad(iprot), 1, name(iprot)
         write(3,50) "CONSTR",ipos-3,ipos-3,
     &      0.0, rad(iprot), 1, name(iprot)
	endif

 50	format(a6,2i8,2f8.3,i2,1x,a3)

	endif
C    from the RNA pol goto
 230 	continue
	enddo
	goto 1

 99	continue
	end
