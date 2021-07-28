C make_const.f
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
C Reads segment file from unit 5
C Writes constraint file to unit 6 with:
C DNA          1   58000                        # defines beads in DNA
C RNA      58001   58049                        # defines beads in each RNA strand
C CONSTR       1   58000   1.000   1.000 1 ORI  # constraints: two bead numbers, distance (lattice units),
C                                                 radius of sphere(s) (lattice units), 1 or 2 spheres, label
C                                                 --this constraint will close the circle
C CONSTR     334   58001   1.000   2.200 1 POL  # --this constraint connects RNA to DNA and adds a sphere for RNApol
C CONSTR   67851   67851   1.000   2.940 1 RIB  # --this constraint adds a large sphere for a ribosome at one position
C
	integer*4 node(10,6),nodeseg(10,0:10)
	integer*4 seg(10,2),segcircle(10),segnode(10,2)
	integer*4 plect(500),plectcircle(500),plectlength(500)
	integer*4 rna1(500,2),rna1root(500),rna1length(500)
	integer*4 rna2(500,2),rna2rib(500),nrib,nrib2
	character*200 line
	character*3 nodename(10)

	iseg=0
	iplect=0
	irna1=0
	irna2=0
	ifork1=0
	ifork2=0
	nrib2=0

 1	read(5,2,end=9) line
 2	format(a120)
	
	if (line(1:4).eq."SEGM") then
	  iseg=iseg+1
	  read(line,10) (seg(iseg,i),i=1,2),(segnode(iseg,j),j=1,2)
 10	  format(4x,2i8,16x,2i2,1x,3i1,a2)
	  goto 1
	endif

	if (line(1:4).eq."CIRC") then
	  iseg=iseg+1
	  read(line,16) ncircle
 16	  format(4x,8x,i8)
	  seg(iseg,1)=1
	  seg(iseg,2)=ncircle
	  goto 1
	endif

	if (line(1:4).eq."PLEC") then
	  iplect=iplect+1
	  read(line,11) plect(iplect),plectcircle(iplect),
     &          plectlength(iplect)
 11	  format(4x,i8,8x,2i8)
	  goto 1
	endif

	if (line(1:4).eq."RNA1") then
	  irna1=irna1+1
	  read(line,12) (rna1(irna1,i),i=1,2),
     &         rna1root(irna1),rna1length(irna1)
 12	  format(4x,4i8)
	  goto 1
	endif

	if (line(1:4).eq."RNA2") then
	  irna2=irna2+1
	  read(line,14) (rna2(irna2,i),i=1,2),nrib,
     &         (rna2rib(nrib2+j),j=1,nrib)
	  nrib2=nrib2+nrib
 14	format(4x,2i8,2x,22i8)
	  goto 1
	endif

	goto 1
  9	continue

	nseg=iseg
	nplect=iplect
	nrna1=irna1
	nrna2=irna2

C assign end of each segment

	do iseg=1,nseg
	write(6,50) "DNA   ",(seg(iseg,i),i=1,2)
	enddo 
	do irna1=1,nrna1
	write(6,50) "RNA   ",(rna1(irna1,i),i=1,2)
	enddo 
	do irna2=1,nrna2
	write(6,50) "RNA   ",(rna2(irna2,i),i=1,2)
	enddo 

C define constraints to connect segments
C hardwired circle constraint
	write(6,50) "CONSTR",1,ncircle,1.,1.,1,"ORI"

 50	format(a6,2i8,2f8.3,i2,1x,a3)

C constraint to attach transcribing RNA
	do irna1=1,nrna1

	iptot=0
	do iplect=1,nplect
	if (plectcircle(iplect).lt.rna1root(irna1)) then
	  iptot=iptot+plectlength(iplect)*2
	else
	  goto 60
	endif
	enddo
 60	continue

C hardwired for 10 bp/bead: d=150 nm -> r=2.2 lattice units
	write(6,50) "CONSTR",rna1root(irna1)+iptot, rna1(irna1,1),
     &      1.0, 2.2, 1, "POL"
	enddo

C add ribosomes
	do irna1=1,nrna1
	if (rna1length(irna1).gt.50) then
	iposition=rna1(irna1,1)+rna1length(irna1)/2
C hardwired for 10 bp/bead: d=200 nm -> r=2.94 lattice units
	endif
	enddo
	do irna2=1,nrib2
	iposition=rna2rib(irna2)
	write(6,50) "CONSTR",iposition,iposition, 1.0, 2.94, 1, "RIB"
	enddo

	end
