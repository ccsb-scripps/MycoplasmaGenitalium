C relax_circle.f
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
C Relaxes coordinates of DNA circle based on constraint file and cellular size
C Reads from unit 5:
C MG149_circle.1.pdb        # input idealized circle coordinates
C MG149_circle.1.const      # input idealized constraint file
C MG149_circle_relax.1.pdb  # output coordinate file
C 2.                        # weight for repulsion
C 1,42.6,42.6               # flag for constraint to cell volume, radius and half length of cell (lattice units)
C 20000,1234                # number of optimization steps, random seed
C
	character*80 line,filename
	character*30 pdb(200000)
	real*4 x(0:200000),y(0:200000),z(0:200000)
	integer*4 type(200000)
	integer*4 model(200000)
	integer*4 contact(0:200000,1000)
	integer*4 ncontact(0:200000)
	real*4 rmax(200000)
	integer RNAflag(200000)
	real*4 rdiag(2)
	integer*4 cspace(-100:100,-100:100,-100:100,50)
	integer*4 nspace(-100:100,-100:100,-100:100)
	integer*4 nbpstep
	integer*4 chain(500,2),circleflag(500),start(200000),end(200000)
	real*4 rconstraint(1500),rconstraintsq(1500)
	character*3 atomname(1500)
	real*4 radius(1500),radiussq(1500),radiuswRNAsq(1500)
	integer*4 nsphere(1500)
	integer*4 constraint(1500,2),imolconst(1500)
	integer*4 hist_constraint(-10:10),hist_sphsph(-10:10)
	integer*4 hist_clash(0:20),hist_bond(-10:10)
	integer*4 hist_sphrna(-10:10),dhist(0:20)
	integer*4 irandatom(5000)
	real*4 rrx(5000),rry(5000),rrz(5000)
	real*4 rxrepulsion(200000),ryrepulsion(200000),rzrepulsion(200000)

C	radius of nucleic acid beads, in lattice units
C	10 bp/bead, diameter of DNA = 2nm = 0.58 lattice units
	rNA=0.3

	read(5,88) filename
	open(2,file=filename)
	write(6,*) "input lattice coordinate file ",filename
	read(5,88) filename
	open(3,file=filename)
	write(6,*) "constraint file ",filename
	read(5,88) filename
	write(6,*) "output coordinate file ",filename
	open(1,file=filename)
 88	format(a80)

C ***  read start and stop for chains, constraints
	nchain=0
	ic=1
	do i=1,200000
	 rmax(i)=0.
	 type(i)=1
	enddo
 33	read(3,88,end=89) line

	if ((line(1:3).eq."DNA").or.(line(1:3).eq."RNA")) then
	  nchain=nchain+1
	  read(line,34) chain(nchain,1),chain(nchain,2)
	endif
C assign bead type=1 for DNA, type=2 for RNA
	if (line(1:3).eq."RNA") then
	do i=chain(nchain,1),chain(nchain,2)
	 type(i)=2
	enddo
	endif

	if (line(1:6).eq."CONSTR") then

	  read(line,76) constraint(ic,1),constraint(ic,2),
     &      rconstraint(ic),radius(ic),nsphere(ic),atomname(ic)
	  imolconst(ic)=ic
	  ic=ic+1
	endif
	goto 33
 89	continue

 	nconstraint=ic-1
	write(6,*) "number of chains ",nchain
	write(6,*) "number of constraints ",nconstraint
	do i=1,nchain
	  write(6,34) chain(i,1),chain(i,2)
	enddo
	do ic=1,nconstraint
	  write(6,76) constraint(ic,1),constraint(ic,2),
     &      rconstraint(ic),radius(ic),nsphere(ic),atomname(ic)
	enddo

 34	format(6x,3i8)
 76	format(6x,2i8,2f8.3,1x,i1,1x,a3)

	do ichain=1,nchain
	write(6,*) "chain ends ",chain(ichain,1),chain(ichain,2)
	enddo
	do ia=1,200000
	start(ia)=0
	end(ia)=0
	enddo
	do ichain=1,nchain
	start(chain(ichain,1))=1
	end(chain(ichain,2))=1
	enddo

C	currently doing everything in lattice units ~10 bp steps
	rscale=1.

	read(5,*) rwrepulsion
	write(6,*) "weight for repulsion: ",rwrepulsion

	read(5,*) icellflag,rend,rcyl
	if (icellflag.ne.0) then
	  write(6,*) "constrain to cell shape"
	  write(6,*) "rad and length (nm) ",rend,rcyl
	  rend=rend/rscale
	  rcyl=rcyl/rscale
	  rcyl=rcyl-rend
	  write(6,*) "rad and cyl half length (lattice unit) ",rend,rcyl
	endif
	rendsq=rend**2
	write(6,*) "rendsq ",rendsq

	read(5,*) nstep,nseed
	write(6,*) "relaxation steps ",nstep
	write(6,*) "random seed ",nseed
	do i=1,nseed
	r=rand()
	enddo

C  *** read atoms

	i=1
	imodel=0
 1	read(2,3,end=9) line
	if (line(1:4).eq."ATOM") then
         read(line,2) pdb(i),x(i),y(i),z(i)
	 model(i)=imodel
	 RNAflag(i)=0
	 if (pdb(i)(18:20).eq."MET") RNAflag(i)=1
	 i=i+1
	goto 1
	endif
 3	format(a80)
	if (line(1:4).eq."MODE") imodel=imodel+1
 2	format(a30,3f8.3,2f6.2)
	goto 1
 9	continue
	natom=i-1

	iclash=0
	write(6,*) "number of lattice points",natom

	write(6,*) "coordinates of start and end of chains"
	do ichain=1,nchain
	istart=chain(ichain,1)
	iend=chain(ichain,2)
	write(6,*) ichain,x(istart),y(istart),z(istart)
	write(6,*) ichain,x(iend),y(iend),z(iend)
	enddo

C **	fix up constraints
	do ic=1,nconstraint
	  rconstraint(ic)=rconstraint(ic)/rscale
	  rconstraintsq(ic)=rconstraint(ic)**2
	  radius(ic)=radius(ic)/rscale
	  radiussq(ic)=radius(ic)**2
C  radius of RNA 0.5 nm
	  rtmp=radius(ic)+rNA/rscale
	  radiuswRNAsq(ic)=rtmp**2
	  rmax(constraint(ic,1))=max(radius(ic),rmax(constraint(ic,1)))
	  rmax(constraint(ic,2))=max(radius(ic),rmax(constraint(ic,2)))
	  rmax(constraint(ic,1)+1)=rmax(constraint(ic,1))
	  rmax(constraint(ic,1)-1)=rmax(constraint(ic,1))
	  rmax(constraint(ic,2)+1)=rmax(constraint(ic,1))
	  rmax(constraint(ic,2)-1)=rmax(constraint(ic,1))
	write(6,*) ic,constraint(ic,1),rmax(constraint(ic,1))
	enddo

c --build contact list--

	do ix=-100,100
	do iy=-100,100
	do iz=-100,100
	nspace(ix,iy,iz)=0
	enddo
	enddo
	enddo
	write(6,*) "zeroed"

	do i=1,natom
	ix=int(x(i))
	iy=int(y(i))
	iz=int(z(i))
	nspace(ix,iy,iz)=nspace(ix,iy,iz)+1
	if (nspace(ix,iy,iz).gt.50) then
            write(6,*) " 1 nspace > 50"
	    nspace(ix,iy,iz)=50
	endif
	cspace(ix,iy,iz,nspace(ix,iy,iz))=i
	enddo

	do i=1,natom
	rxrepulsion(i)=0.
	ryrepulsion(i)=0.
	rzrepulsion(i)=0.
	enddo

	write(6,*) "finished first contact list"

	iramp=0
	itoggle=0

	rnstep=float(nstep)
	irandcount=1

	do istep=1,nstep


	iramp=iramp+1
	if (iramp.eq.200) then
	write(6,*) "cycle ",istep

	write(6,568) "        ",(float(idistance)/10.,idistance=-10,10)
	write(6,555) "constr  ",(hist_constraint(ih),ih=-10,10)
	write(6,555) "sphsph  ",(hist_sphsph(ih),ih=-10,10)
	write(6,555) "sphDNA  ",(hist_sphrna(ih),ih=-10,10)
	write(6,*)
	write(6,568) "        ",(float(idistance)/10.,idistance=0,20)
	write(6,555) "DNADNA  ",(hist_clash(ih),ih=0,20)
	write(6,555) "bonds   ",(hist_bond(ih),ih=0,20)
 555	format(a8,i8,9i8,i8,9i8,i8)
 568	format(a8,21f8.1)
 556	format(a8,21i8,i8)
 78	format(a17,3f7.4,2x,2f9.6)
	iramp=0

	maxcontact=0
	do ix=-100,100
	do iy=-100,100
	do iz=-100,100
	nspace(ix,iy,iz)=0
	enddo
	enddo
	enddo

	do i=1,natom
	ix=x(i)
	iy=y(i)
	iz=z(i)
	nspace(ix,iy,iz)=nspace(ix,iy,iz)+1
	if (nspace(ix,iy,iz).gt.48) then
	   write(6,*) "nspace ",nspace(ix,iy,iz),ix,iy,iz
	   write(6,*) i,natom,x(i),y(i),z(i)
	endif
	cspace(ix,iy,iz,nspace(ix,iy,iz))=i
	enddo

C	calculate simple mutual repulsion vector
	do i=1,natom
	rxsum=0.
	rysum=0.
	rzsum=0.
	rxrepulsion(i)=0.
	ryrepulsion(i)=0.
	rzrepulsion(i)=0.

	do j=1,natom
	if (abs(i-j).gt.10) then
	 rxdiff=x(j)-x(i)
	 rydiff=y(j)-y(i)
	 rzdiff=z(j)-z(i)
	 rweight=(rxdiff**2+rydiff**2+rzdiff**2)
	if (rweight.lt.42.) then
	     rxsum=rxsum+rxdiff/rweight
	     rysum=rysum+rydiff/rweight
	     rzsum=rzsum+rzdiff/rweight
	endif
	endif

C jatom loop
	enddo

	rxrepulsion(i)=rxsum/rwrepulsion
	ryrepulsion(i)=rysum/rwrepulsion
	rzrepulsion(i)=rzsum/rwrepulsion
	enddo

	do id=1,10
	dhist(id)=0
	enddo
	do ix=-100,100
	do iy=-100,100
	do iz=-100,100
	id=min(nspace(ix,iy,iz),10)
	dhist(id)=dhist(id)+1
	enddo
	enddo
	enddo
	write(6,*)
	write(6,130) "dhist  ",(dhist(id),id=1,10)
	write(6,*) "--------------------------------------------"
 130	format(a6,10i8)

	endif

C bond constraints
	 rbond=0.010
	 rllow=(0.95+float(istep)*0.04/rnstep)**2
	 rlhigh=(1.05-float(istep)*0.04/rnstep)**2

C radius of DNA=1.0 at end
c rclash is 2X radius of beads, squared
	 rclash=(rNA*2.*0.75+float(istep)*(rNA*2.*0.25)/rnstep)**2

C angle constraint
C   DNA
	 rdiag(1)=(1.8+float(istep)*0.1/rnstep)**2
C   RNA
	 rdiag(2)=(1.0+float(istep)*0.1/rnstep)**2

C constraint tolerances
	 rclow=-0.05+float(istep)*0.04/rnstep
	 rchigh=+0.05-float(istep)*0.04/rnstep
	 rcstep=0.040

C clash constraints

	rsphsph=0.004
	rsphrna=0.004
	rrnarna=0.004
	rmutual=0.004-float(istep)*0.0038/rnstep

C random step
	 rstep=0.15-float(istep)*0.140/rnstep
C reduction for local random move
	rlocal=0.1

C apply input constraints
	do ihist=-10,10
	hist_constraint(ihist)=0
	hist_sphsph(ihist)=0
	hist_sphrna(ihist)=0
	enddo
	do ihist=0,20
	hist_clash(ihist)=0
	hist_bond(ihist)=0
	enddo

	do ic=1,nconstraint

	ia1=constraint(ic,1)
	ia2=constraint(ic,2)

	if (ia1.ne.ia2) then
	rxv=x(ia1)-x(ia2)
	ryv=y(ia1)-y(ia2)
	rzv=z(ia1)-z(ia2)
	rdcsq=rxv**2+ryv**2+rzv**2
	if (rdcsq.eq.0) write(6,*) "rdcsq zero",ia1,ia2,
     &   x(ia1),x(ia2),y(ia1),y(ia2),z(ia1),z(ia2)
	rdc=sqrt(rdcsq)
	
	if (iramp.eq.199) then
	rdiff=rdc-rconstraint(ic)
	idiff=int(rdiff*10.)
	idiff=max(idiff,-10)
	idiff=min(idiff,10)
	hist_constraint(idiff)=hist_constraint(idiff)+1
	endif

	if (rdc.lt.rconstraint(ic)+rclow) then
	rxv=rxv/rdc*rcstep
	ryv=ryv/rdc*rcstep
	rzv=rzv/rdc*rcstep
	x(ia1)=x(ia1)+rxv
	y(ia1)=y(ia1)+ryv
	z(ia1)=z(ia1)+rzv
	x(ia2)=x(ia2)-rxv
	y(ia2)=y(ia2)-ryv
	z(ia2)=z(ia2)-rzv
	endif
	if (rdc.gt.rconstraint(ic)+rchigh) then
	rxv=rxv/rdc*rcstep
	ryv=ryv/rdc*rcstep
	rzv=rzv/rdc*rcstep

	x(ia1)=x(ia1)-rxv
	y(ia1)=y(ia1)-ryv
	z(ia1)=z(ia1)-rzv
	x(ia2)=x(ia2)+rxv
	y(ia2)=y(ia2)+ryv
	z(ia2)=z(ia2)+rzv

cenddo
	endif

	endif

	enddo

C sphere-sphere overlap for constraints
	do ic=1,nconstraint
	if (radius(ic).ne.0) then

	ia1=constraint(ic,1)
	ia2=constraint(ic,2)
	rxc=(x(ia1)+x(ia2))/2.
	ryc=(y(ia1)+y(ia2))/2.
	rzc=(z(ia1)+z(ia2))/2.

	do ic2=1,nconstraint
	if ((ic.ne.ic2).and.(radius(ic2).ne.0.).and.
     &      (imolconst(ic).ne.imolconst(ic2))) then

	ib1=constraint(ic2,1)
	ib2=constraint(ic2,2)
	rxc2=(x(ib1)+x(ib2))/2.
	ryc2=(y(ib1)+y(ib2))/2.
	rzc2=(z(ib1)+z(ib2))/2.

	rxv=rxc2-rxc
	ryv=ryc2-ryc
	rzv=rzc2-rzc
	rdcsq=(rxv**2+ryv**2+rzv**2)
	if (rdcsq.eq.0) write(6,*) "rdcsq s-s zero",
     &   ia1,ia2,ib1,ib2,
     &    rxc,ryc,rzc,rxc2,ryc2,rzc2
	rdc=sqrt(rdcsq)

	if (iramp.eq.199) then
	rdiff=rdc-radius(ic)-radius(ic2)
	idiff=int(rdiff*10.)
	idiff=max(idiff,-10)
	idiff=min(idiff,10)
	hist_sphsph(idiff)=hist_sphsph(idiff)+1
	endif

	if (rdc.lt.radius(ic)+radius(ic2)) then
	rxv=rxv/rdc*rsphsph
	ryv=ryv/rdc*rsphsph
	rzv=rzv/rdc*rsphsph
	x(ia1)=x(ia1)-rxv
	x(ia2)=x(ia2)-rxv
	y(ia1)=y(ia1)-ryv
	y(ia2)=y(ia2)-ryv
	z(ia1)=z(ia1)-rzv
	z(ia2)=z(ia2)-rzv
	x(ib1)=x(ib1)+rxv
	x(ib2)=x(ib2)+rxv
	y(ib1)=y(ib1)+ryv
	y(ib2)=y(ib2)+ryv
	z(ib1)=z(ib1)+rzv
	z(ib2)=z(ib2)+rzv
	endif

	endif	
	enddo
	endif	
	enddo

C small random moves
	ratom=float(natom)-6.

	do irand=1,2000

	ia=irandatom(irand)

	irandatom(irand)=int(rand()*ratom)+3
	rrandx=rand()-0.5
	rrandy=rand()-0.5
	rrandz=rand()-0.5
	rd=sqrt(rrandx**2+rrandy**2+rrandz**2)
	if (rd.eq.0) write(6,*) "rd zero 2"
	rrandx=rrandx/rd
	rrandy=rrandy/rd
	rrandz=rrandz/rd

	ia=irandatom(irand)

	rxv=x(ia-1)-x(ia+1)
	ryv=y(ia-1)-y(ia+1)
	rzv=z(ia-1)-z(ia+1)
	rd=sqrt(rxv**2+ryv**2+rzv**2)
	if (rd.eq.0) write(6,*) "rd zero 3"
	rxv=rxv/rd
	ryv=ryv/rd
	rzv=rzv/rd
	rxcross=ryv*rrandz-rzv*rrandy
	rycross=rzv*rrandx-rxv*rrandz
	rzcross=rxv*rrandy-ryv*rrandx
	rd=sqrt(rxcross**2+rycross**2+rzcross**2)
	rrx(irand)=rxcross/rd*rstep
	rry(irand)=rycross/rd*rstep
	rrz(irand)=rzcross/rd*rstep

cendif

C apply random move to local segment
	iseg=1
	rseg=float(iseg)+1.
	do ioff=-iseg,iseg
	roff=(rseg-(float(abs(ioff))))/rseg
C this is normal to chain direction
	x(ia+ioff)=x(ia+ioff)+rrx(irand)*roff*rlocal
	y(ia+ioff)=y(ia+ioff)+rry(irand)*roff*rlocal
	z(ia+ioff)=z(ia+ioff)+rrz(irand)*roff*rlocal
C this is isotropic random
cx(ia+ioff)=x(ia+ioff)+rrandx*roff*rstep
cy(ia+ioff)=y(ia+ioff)+rrandy*roff*rstep
cz(ia+ioff)=z(ia+ioff)+rrandz*roff*rstep
	enddo

	enddo
 876	format(3i6,9f9.5)

C *** step through RNA ****
	do ia=1,natom

C input constraints as spherical clash with RNA
	rmin=999.
	do ic=1,nconstraint
	if (radius(ic).ne.0) then
	ia1=constraint(ic,1)
	ia2=constraint(ic,2)
	if ((ia.ne.ia1).and.(ia.ne.ia2)) then
	rxc=(x(ia1)+x(ia2))/2.
	ryc=(y(ia1)+y(ia2))/2.
	rzc=(z(ia1)+z(ia2))/2.
	rxv=x(ia)-rxc
	ryv=y(ia)-ryc
	rzv=z(ia)-rzc
	rdcsq=rxv**2+ryv**2+rzv**2

	if (iramp.eq.199) then
	if ((abs(ia1-ia).gt.2.).and.(abs(ia2-ia).gt.2.)) then
	rmin=min(rmin,sqrt(rdcsq)-sqrt(radiuswRNAsq(ic)))
	endif
	endif

	if (rdcsq.lt.radiuswRNAsq(ic)) then
	x(ia)=x(ia)+rxv*rsphrna
	y(ia)=y(ia)+ryv*rsphrna
	z(ia)=z(ia)+rzv*rsphrna
	endif

	endif	
	endif	
	enddo

	if (iramp.eq.199) then
	idiff=int(rmin*10.)
	idiff=max(idiff,-10)
	idiff=min(idiff,10)
	hist_sphrna(idiff)=hist_sphrna(idiff)+1
	endif

C bond constraint

	if ((end(ia).eq.0).and.(itoggle.eq.1)) then
	itoggle=2
	rxv=x(ia)-x(ia+1)
	ryv=y(ia)-y(ia+1)
	rzv=z(ia)-z(ia+1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rllow) then
	x(ia)=x(ia)+rxv*rbond
	y(ia)=y(ia)+ryv*rbond
	z(ia)=z(ia)+rzv*rbond
	endif
	if (rd3sq.gt.rlhigh) then
	x(ia)=x(ia)-rxv*rbond
	y(ia)=y(ia)-ryv*rbond
	z(ia)=z(ia)-rzv*rbond
	endif
	endif

	if (start(ia).eq.0) then
	rxv=x(ia)-x(ia-1)
	ryv=y(ia)-y(ia-1)
	rzv=z(ia)-z(ia-1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rllow) then
	x(ia)=x(ia)+rxv*rbond
	y(ia)=y(ia)+ryv*rbond
	z(ia)=z(ia)+rzv*rbond
	endif
	if (rd3sq.gt.rlhigh) then
	x(ia)=x(ia)-rxv*rbond
	y(ia)=y(ia)-ryv*rbond
	z(ia)=z(ia)-rzv*rbond
	endif
	endif

	if ((end(ia).eq.0).and.(itoggle.eq.0)) then
	itoggle=1
	rxv=x(ia)-x(ia+1)
	ryv=y(ia)-y(ia+1)
	rzv=z(ia)-z(ia+1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rllow) then
	x(ia)=x(ia)+rxv*rbond
	y(ia)=y(ia)+ryv*rbond
	z(ia)=z(ia)+rzv*rbond
	endif
	if (rd3sq.gt.rlhigh) then
	x(ia)=x(ia)-rxv*rbond
	y(ia)=y(ia)-ryv*rbond
	z(ia)=z(ia)-rzv*rbond
	endif
	endif
	if (itoggle.eq.2) itoggle=0

	if (iramp.eq.199) then
	ibond=int(sqrt(rd3sq)*10.)
	ibond=min(ibond,20)
	hist_bond(ibond)=hist_bond(ibond)+1
	endif
C end bond constraint

C *** bond angles
	if ((start(ia).eq.0).and.(end(ia).eq.0)) then

	rxv=x(ia+1)-x(ia-1)
	ryv=y(ia+1)-y(ia-1)
	rzv=z(ia+1)-z(ia-1)
	rd3sq=rxv**2+ryv**2+rzv**2
	if (rd3sq.lt.rdiag(type(ia))) then
	x(ia-1)=x(ia-1)-rxv*0.004
	y(ia-1)=y(ia-1)-ryv*0.004
	z(ia-1)=z(ia-1)-rzv*0.004
	x(ia+1)=x(ia+1)+rxv*0.004
	y(ia+1)=y(ia+1)+ryv*0.004
	z(ia+1)=z(ia+1)+rzv*0.004
	endif

	rxv=x(ia+1)-x(ia-1)
        ryv=y(ia+1)-y(ia-1)
        rzv=z(ia+1)-z(ia-1)
        rd3sq=rxv**2+ryv**2+rzv**2
        if (rd3sq.lt.rdiag(type(ia))) then
        rxa=(x(ia+1)+x(ia-1))/2.
        rya=(y(ia+1)+y(ia-1))/2.
        rza=(z(ia+1)+z(ia-1))/2.
        x(ia)=x(ia)+(rxa-x(ia))*0.005
        y(ia)=y(ia)+(rya-y(ia))*0.005
        z(ia)=z(ia)+(rza-z(ia))*0.005
        x(ia+1)=x(ia+1)-(rxa-x(ia))*0.001
        y(ia+1)=y(ia+1)-(rya-y(ia))*0.001
        z(ia+1)=z(ia+1)-(rza-z(ia))*0.001
        x(ia-1)=x(ia-1)-(rxa-x(ia))*0.001
        y(ia-1)=y(ia-1)-(rya-y(ia))*0.001
        z(ia-1)=z(ia-1)-(rza-z(ia))*0.001
        endif

	endif
C *** end bond angles
c --check for clash
	iclash=0
	rmin=999.
C--this search box needs to be large enough to have at least two strands
	do ixo=-2,2
	do iyo=-2,2
	do izo=-2,2
	ix=x(ia)+ixo
	iy=y(ia)+iyo
	iz=z(ia)+izo
	do icontact=1,nspace(ix,iy,iz)
	ic=cspace(ix,iy,iz,icontact)
	if ((ic.lt.ia-2).or.(ic.gt.ia+2)) then
	rcsq=(x(ia)-x(ic))**2+(y(ia)-y(ic))**2+(z(ia)-z(ic))**2
	if (rcsq.lt.rclash) iclash=1
	if (rcsq.lt.rmin) then
	 rvx=x(ia)-x(ic)
	 rvy=y(ia)-y(ic)
	 rvz=z(ia)-z(ic)
	 rmin=rcsq
	endif
	endif
	enddo
	enddo
	enddo
	enddo

	if (iramp.eq.199) then
	imin=int(sqrt(rmin)*10.)
	imin=min(imin,20)
	hist_clash(imin)=hist_clash(imin)+1
	endif

	if (iclash.eq.1) then
	rnorm=sqrt(rvx**2+rvy**2+rvz**2)
	if (rnorm.eq.0.) rnorm=1.
	 x(ia)=x(ia)+rrnarna*rvx/rnorm
	 y(ia)=y(ia)+rrnarna*rvy/rnorm
	 z(ia)=z(ia)+rrnarna*rvz/rnorm
	
C this block adds a small motion to flanking regions with clash
	 if (start(ia).eq.0) then
	 x(ia-1)=x(ia-1)+rrnarna*rvx/rnorm/3.
	 y(ia-1)=y(ia-1)+rrnarna*rvy/rnorm/3.
	 z(ia-1)=z(ia-1)+rrnarna*rvz/rnorm/3.
	 endif
	 if (end(ia).eq.0) then
	 x(ia+1)=x(ia+1)+rrnarna*rvx/rnorm/3.
	 y(ia+1)=y(ia+1)+rrnarna*rvy/rnorm/3.
	 z(ia+1)=z(ia+1)+rrnarna*rvz/rnorm/3.
	 endif
	endif

C apply simple mutual repulsion
	if (RNAflag(ia).eq.0) then
	x(ia)=x(ia)-rxrepulsion(ia)*rmutual
	y(ia)=y(ia)-ryrepulsion(ia)*rmutual
	z(ia)=z(ia)-rzrepulsion(ia)*rmutual
	endif

C --constrain to cell sphere
	if (icellflag.ne.0) then

	radcell=sqrt(x(ia)**2+y(ia)**2+z(ia)**2)
	if (radcell.gt.rend-rmax(ia)) then
	  x(ia)=x(ia)-(x(ia)/radcell)*0.010
	  y(ia)=y(ia)-(y(ia)/radcell)*0.010
	  z(ia)=z(ia)-(z(ia)/radcell)*0.010
	else
C create a softer edge for constraint
	if (radcell.gt.rend-rmax(ia)-10.) then
	 rweight=0.002*(radcell-(rend-rmax(ia))+10.)/10.
	  x(ia)=x(ia)-(x(ia)/radcell)*rweight
	  y(ia)=y(ia)-(y(ia)/radcell)*rweight
	  z(ia)=z(ia)-(z(ia)/radcell)*rweight
	endif
	endif

	endif

	enddo
	enddo

	modellast=0
	ireslast=0
	imodel=1
	write(1,12) "MODEL ", imodel
	do i=1,natom

C this fails if residue numbers are messed up--need to fix
 567	format(22x,i4)

	if (pdb(i)(18:20).eq."ARG") modellast=1
	if ((pdb(i)(18:20).eq."GLY").and.(modellast.eq.1)) then
	  write(1,12) "ENDMDL"
	  imodel=imodel+1
	  write(1,12) "MODEL ", imodel
 	  modellast=0
	endif

	write(1,2) pdb(i),x(i)*rscale,y(i)*rscale,z(i)*rscale,
     &     rclash/2.*rscale
	enddo
	write(1,12) "ENDMDL"

	imodel=imodel+1
	write(1,12) "MODEL ", imodel
	do ic=1,nconstraint
	ia1=constraint(ic,1)
	ia2=constraint(ic,2)
	rxc=(x(ia1)+x(ia2))/2.
	ryc=(y(ia1)+y(ia2))/2.
	rzc=(z(ia1)+z(ia2))/2.
	pdb(ic)(18:20)=atomname(ic)
	pdb(ic)(1:6)="HETATM"
	pdb(ic)(14:16)="PRO"
	if (radius(ic).ne.0) then
	if (nsphere(ic).eq.1) then
	write(1,81) "HETATM",ic,atomname(ic),"PRO","P",ic,
     &     rxc*rscale,ryc*rscale,rzc*rscale,
     &     rconstraint(ic)*rscale
 81	format(a6,i5,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2)
	else
	pdb(ic)(14:16)=atomname(ic)
	rxc=x(ia1)+(x(ia2)-x(ia1))*0.25
	ryc=y(ia1)+(y(ia2)-y(ia1))*0.25
	rzc=z(ia1)+(z(ia2)-z(ia1))*0.25
	write(1,81) "HETATM",ic,atomname(ic),"PRO","P",ic,
     &     rxc*rscale,ryc*rscale,rzc*rscale,
     &     rconstraint(ic)*rscale
	rxc=x(ia1)+(x(ia2)-x(ia1))*0.75
	ryc=y(ia1)+(y(ia2)-y(ia1))*0.75
	rzc=z(ia1)+(z(ia2)-z(ia1))*0.75
	write(1,81) "HETATM",ic,atomname(ic),"PRO","P",ic,
     &     rxc*rscale,ryc*rscale,rzc*rscale,
     &     rconstraint(ic)*rscale
	endif
	endif

	enddo
	write(1,12) "ENDMDL"
	
 12	format(a6,i8)

	end
