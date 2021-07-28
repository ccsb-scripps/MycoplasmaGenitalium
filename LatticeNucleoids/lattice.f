C lattice.f
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
C Creates lattice model of DNA and RNA
C Reads from unit 5:
C MG149_lattice.1.pdb       # output coordinate file
C MG149_lattice.1.chain     # output chain specification file
C 2                         # number of passes for adding plectonemes
C MG149_circle_relax.1.pdb  # input circle coordinate file
C MG149.segments            # input segment file
C 1234                      # random seed
C 60,5                      # percent bias for straight plectoneme and RNA (0-100, high value=straighter)
C 6                         # lattice units per supercoil
C 42.6,42.6                 # radius and half length of cell (lattice units)
C
C Output coordinate file is in lattice units and uses the RESNAME to denote different parts of the model:
C DNA: GLY is non-plectonemic regions, ARG and LYS are two sides of plectonemic supercoil
C RNA: MET is nascent and CYS is free 
C
	character*80 filename,line
	integer*4 ringlast
	integer*4 plectlength(0:5000),plectroot(0:5000),plectflag(0:5000)
	integer*4 plectrootout(0:5000)
	integer*4 segout(0:100,2),plectout(0:5000)
	integer*4 coord(0:200000,3)
	real*4 rcoord(0:200000,3)
	integer*4 space(-120:120,-60:60,-60:60)
	integer*4 pspace(10,-120:120,-60:60,-60:60)
	integer*4 ptest(-120:120,-60:60,-60:60)
	integer*4 plect(5000,0:18000,3),plectspace(5000)
	integer*4 direction(0:18)
	integer*4 istraight(3)
	character chain(0:200000),segchain,pchain,chainout
	real*4 pv(0:200000,3),pnext(0:200000,3),pnormal(0:200000,3),pi

	pi=3.14159

	read(5,82) filename
 82	format(a80)
	open(2,file=filename)
	write(6,*) "output pdb file ",filename
	read(5,82) filename
	open(10,file=filename)
	write(6,*) "output dat file ",filename

	read(5,*) npass

c--input coordinates of circle
	read(5,82) filename
	open(20,file=filename)
	write(6,*) "input coordinates of non-plectonemic DNA: ",filename
c---plectoneme definition file
	read(5,82) filename
	open(3,file=filename)
	write(6,*) "plectoneme definition file: ",filename

	iseg=1
	nseg=0
	istartlast=1
	iendtot=999999
	ip=1
 345	read(3,82,end=346) line

	if ((line(1:4).eq."SEGM").or.(line(1:4).eq."CIRC")) then
	read(line,83) segout(iseg,1),segout(iseg,2),istart,iend,
     &      segchain
 83	format(4x,4i8,9x,a1)
	  do ia=segout(iseg,1),segout(iseg,2)
	  chain(ia)=segchain
	  enddo
	npt=iend-istart+1
	ncircle=segout(iseg,2)
	iendtot=npt
	write(6,*) "segment ",istart,iend,iendtot,
     &     (segout(iseg,i),i=1,2),segchain
	iseg=iseg+1
	nseg=nseg+1
	endif

	ifork1=0
	ifork2=0
	if (line(1:4).eq."FORK") then
	read(line,81) ifork1,ifork2
 81	format(4x,2i8)
	write(6,*) "fork numbers ",ifork1,ifork2,ncircle
	endif
C	if there is a fork, overwrite chain end with end of daughter strand
C	end of circle will be written because of first plectoneme in daughter
	if (line(1:4).eq."FOR2") then
	read(line,87) iendtot
 87	format(4x,24x,i8)
	endif

	if (line(1:4).eq."PLEC") then
	read(line,84) plectout(ip),plectroot(ip),plectlength(ip)
 84	format(4x,i8,8x,2i8)
	plectflag(ip)=1
	istartlast=max(plectroot(ip),istartlast)
	write(6,*) line(1:4),plectroot(ip),plectlength(ip),istartlast
	ip=ip+1
	endif

	if (line(1:4).eq."RNA1") then
	read(line,85) plectout(ip),ipend,plectroot(ip),plectlength(ip),
     &     pchain
 85	format(4x,4i8,1x,a1)
	  do ia=plectout(ip),ipend
	  chain(ia)=pchain
	  enddo
	plectflag(ip)=2
	iendtot=min(plectout(ip),iendtot)
	write(6,*) line(1:4),plectroot(ip),plectlength(ip),iendtot
	ip=ip+1
	endif

	if (line(1:4).eq."RNA2") then
	read(line,86) plectout(ip),ipend,pchain
	plectlength(ip)=ipend-plectout(ip)+1
 86	format(4x,2i8,1x,a1)
	  do ia=plectout(ip),ipend
	  chain(ia)=pchain
	  enddo
	plectflag(ip)=3
	plectroot(ip)=0
	iendtot=min(plectout(ip),iendtot)
	plectout(ip)=0
	write(6,*) line(1:4),plectout(ip),plectlength(ip),iendtot
	ip=ip+1
	endif


	goto 345
 346	nplect=ip-1


	write(6,*) "istartlast/iendtot ",istartlast,iendtot

	plectflag(0)=0
	plectroot(0)=1
	plectlength(0)=0
	do i=0,200
	plectrootout(i)=0
	enddo

	read(5,*) irandstart
	write(6,*) "random seed ",irandstart
	do ir=1,irandstart
	r=rand()
	enddo

	write(6,*) "number of plectonemes ", nplect
	write(6,*) "length of plectoneme (lattice units) ",
     &    (plectlength(ip),ip=1,nplect)
	write(6,*) "position of plectoneme roots ",
     &    (plectroot(ip),ip=1,nplect)

	read(5,*) istraight(1),istraight(2)
	istraigth=min(istraight(1),100)
	istraigth=max(istraight(1),0)
	istraigth=min(istraight(2),100)
	istraigth=max(istraight(2),0)
	write(6,*) "percent bias for straight plectoneme (DNA,RNA) ",
     &    istraight(1),istraight(2)
	istraight(3)=istraight(2)

	read(5,*) nsupercoil_bp
	rsupercoil=360./float(nsupercoil_bp)
	write(6,*) "number of lattice units per supercoil ",nsupercoil_bp
	write(6,*) "angle per lattice step (deg) ",rsupercoil
	
	read(5,*) rcell,rlcell

	write(6,*) "radius and half length of cell (lattice units) ",
     &             rcell,rlcell

	if (rcell.gt.59.) then
	write(6,*) "radius too large, setting to 59"
	rcell=59.
	endif
	if (rlcell.gt.119.) then
	write(6,*) "half length too large, setting to 119"
	rlcell=119.
	endif
c -- half length of center cylinder
	rend=rcell
	rcyl=rlcell-rcell
	ilength=int(rlcell)+1
	iwidth=int(rcell)+1
	write(6,*) "rend,rcyl",rend,rcyl

C----Jump to new read in of circle

	direction(0)=1
	direction(1)=2
	direction(2)=3
	direction(3)=4
	direction(4)=5
	direction(5)=6
	direction(6)=1
	direction(7)=2
	direction(8)=3
	direction(9)=4
	direction(10)=5
	direction(11)=6
	direction(12)=1
	direction(13)=2
	direction(14)=3
	direction(15)=4
	direction(16)=5
	direction(17)=6
	direction(18)=1

C----read in circle----
	iring=1
 457	read(20,455,end=458) line
 455	format(a80)
	if ((line(1:4).eq."ATOM").and.(line(14:14).eq."C")) then
	read(line,89) (rcoord(iring,j),j=1,3)
	coord(iring,1)=int(rcoord(iring,1))
	coord(iring,2)=int(rcoord(iring,2))
	coord(iring,3)=int(rcoord(iring,3))
	ix=int(coord(iring,1))
	iy=int(coord(iring,2))
	iz=int(coord(iring,3))
	iring=iring+1
	endif
	goto 457
 89	format(30x,3f8.3)
 458	ringlast=iring
	write(6,*) "input loop size:", ringlast
	write(6,*)
 10	format(41i3)

C---- calculate exclusion map and map in circle
	do i=-ilength,ilength
	do j=-iwidth,iwidth
	do k=-iwidth,iwidth
	r1=sqrt((float(j))**2+float(k)**2)
	r2=sqrt((float(i)-rcyl)**2+float(j)**2+float(k)**2)
	r3=sqrt((float(i)+rcyl)**2+float(j)**2+float(k)**2)
	space(i,j,k)=0
	if (r1.gt.rend) space(i,j,k)=999
	if ((r2.gt.rend).and.(i.gt.rcyl)) space(i,j,k)=999
	if ((r3.gt.rend).and.(i.lt.-rcyl)) space(i,j,k)=999
	enddo
	enddo
	enddo
	write(6,*) "RINGLAST(number of input coordinates) ",ringlast
	do i=1,ringlast

	space(coord(i,1),coord(i,2),coord(i,3))=99

	enddo
	do i=1,3
	coord(0,i)=coord(ringlast,i)
	coord(ringlast+1,i)=coord(1,i)
	enddo

	write(6,*) "Exclusion map for cell"
	do iy=-60,60,8
	write(6,10) (space(ix,iy,0),ix=-120,120,8)
	enddo
C---- plectonemes ---
	do j=1,3
	coord(0,j)=coord(ringlast,j)
	coord(ringlast+1,j)=coord(1,j)
	coord(ringlast+2,j)=coord(2,j)
	enddo
	write(6,*) "starting plectonemes"

	do i=1,nplect
	plectspace(i)=1
	enddo

	do i=-ilength,ilength
	do j=-iwidth,iwidth
	do k=-iwidth,iwidth
	do is=1,3
	if (space(i,j,k).ne.999) pspace(is,i,j,k)=0
	if (space(i,j,k).eq.999) pspace(is,i,j,k)=999
	enddo
	enddo
	enddo
	enddo

	ispacemax=-99

	do ipass=1,npass

	do ip=1,nplect

	plect(ip,0,1)=coord(plectroot(ip),1)
	plect(ip,0,2)=coord(plectroot(ip),2)
	plect(ip,0,3)=coord(plectroot(ip),3)
c--- free RNA
	if (plectflag(ip).eq.3) then
 887	ixrand=int((rand()-0.5)*ilength*2)
	iyrand=int((rand()-0.5)*iwidth*2)
	izrand=int((rand()-0.5)*iwidth*2)
	iclocal=0
	do i=-2,2
	do j=-2,2
	do k=-2,2
	if (space(ixrand+i,iyrand+j,izrand+k).ne.0) iclocal=iclocal+1
	enddo
	enddo
	enddo
	if (iclocal.gt.3) goto 887
	plect(ip,0,1)=ixrand
	plect(ip,0,2)=iyrand
	plect(ip,0,3)=izrand
	endif
	nlength=plectlength(ip)
	
	write(6,*) "starting plectoneme ",ip,(plect(ip,0,j),j=1,3),nlength
	do ispace=1,3

	if (ispace.ne.1) then 
	 write(6,*) "trying space ",ispace
	endif

	ispacemax=max(ispacemax,ispace)
	plectspace(ip)=ispace

	do ipzero=1,nlength
	if (ipass.gt.1) then
	pspace(ispace,plect(ip,ipzero,1),plect(ip,ipzero,2),
     &     plect(ip,ipzero,3))=0.
	endif
	enddo

	do iptry=0,1000

	do ipzero=1,nlength
	plect(ip,ipzero,1)=0.
	plect(ip,ipzero,2)=0.
	plect(ip,ipzero,3)=0.
	enddo

	do i=-ilength,ilength
	do j=-iwidth,iwidth
	do k=-iwidth,iwidth
	ptest(i,j,k)=0.
	if (space(i,j,k).ne.999) ptest(i,j,k)=0
	if (space(i,j,k).eq.999) ptest(i,j,k)=999
	enddo
	enddo
	enddo

	do iproot=1,nplect
	plect(iproot,0,1)=coord(plectroot(iproot),1)
	plect(iproot,0,2)=coord(plectroot(iproot),2)
	plect(iproot,0,3)=coord(plectroot(iproot),3)
c--- free RNA
	if (plectflag(iproot).eq.3) then
 886	ixrand=int((rand()-0.5)*ilength*2)
	iyrand=int((rand()-0.5)*iwidth*2)
	izrand=int((rand()-0.5)*iwidth*2)
	iclocal=0
	do i=-2,2
	do j=-2,2
	do k=-2,2
	if (space(ixrand+i,iyrand+j,izrand+k).ne.0) iclocal=iclocal+1
	enddo
	enddo
	enddo
	if (iclocal.gt.3) goto 886
	plect(iproot,0,1)=ixrand
	plect(iproot,0,2)=iyrand
	plect(iproot,0,3)=izrand
	endif
c -- make sure plectonemes in higher levels don't follow circle at root
	ptest(coord(plectroot(iproot),1),
     &        coord(plectroot(iproot),2),
     &        coord(plectroot(iproot),3))=999
	ptest(coord(plectroot(iproot)+1,1),
     &        coord(plectroot(iproot)+1,2),
     &        coord(plectroot(iproot)+1,3))=999
	ptest(coord(plectroot(iproot)-1,1),
     &        coord(plectroot(iproot)-1,2),
     &        coord(plectroot(iproot)-1,3))=999

	enddo
	
	do istep=0,nlength

	ix=plect(ip,istep,1)
	iy=plect(ip,istep,2)
	iz=plect(ip,istep,3)
	
c --	this will prefer random directions in plectonemes
	idir=int(rand()*6.)
	isuccess=0
c --this will prefer straight plectonemes
	itest=int(rand()*100.)
	if (itest.lt.istraight(plectflag(ip))) then
	if ((plect(ip,istep,1).eq.plect(ip,istep-1,1)+1)) idir=0
	if ((plect(ip,istep,1).eq.plect(ip,istep-1,1)-1)) idir=1
	if ((plect(ip,istep,2).eq.plect(ip,istep-1,2)+1)) idir=2
	if ((plect(ip,istep,2).eq.plect(ip,istep-1,2)-1)) idir=3
	if ((plect(ip,istep,3).eq.plect(ip,istep-1,3)+1)) idir=4
	if ((plect(ip,istep,3).eq.plect(ip,istep-1,3)-1)) idir=5
	endif

	do ircount=1,200
	if (ircount.ne.1) idir=int(rand()*6.)
	ir=direction(idir)

	if (ir.eq.1) then
	  if ((ispace.eq.1).and.(space(ix+1,iy,iz).ne.0)) goto 707
	  if (pspace(ispace,ix+1,iy,iz).ne.0) goto 707
	  if (ptest(ix+1,iy,iz).ne.0) goto 707
	   ix=ix+1
	   ptest(ix,iy,iz)=ip
	   isuccess=1
	   goto 700
	endif

	if (ir.eq.2) then
	  if ((ispace.eq.1).and.(space(ix-1,iy,iz).ne.0)) goto 707
	  if (pspace(ispace,ix-1,iy,iz).ne.0) goto 707
	  if (ptest(ix-1,iy,iz).ne.0) goto 707
	   ix=ix-1
	   ptest(ix,iy,iz)=ip
	   isuccess=1
	   goto 700
	endif

	if (ir.eq.3) then
	  if ((ispace.eq.1).and.(space(ix,iy+1,iz).ne.0)) goto 707
	  if (pspace(ispace,ix,iy+1,iz).ne.0) goto 707
	  if (ptest(ix,iy+1,iz).ne.0) goto 707
	   iy=iy+1
	   ptest(ix,iy,iz)=ip
	   isuccess=1
	   goto 700
	endif

	if (ir.eq.4) then
	  if ((ispace.eq.1).and.(space(ix,iy-1,iz).ne.0)) goto 707
	  if (pspace(ispace,ix,iy-1,iz).ne.0) goto 707
	  if (ptest(ix,iy-1,iz).ne.0) goto 707
	   iy=iy-1
	   ptest(ix,iy,iz)=ip
	   isuccess=1
	   goto 700
	endif

	if (ir.eq.5) then
	  if ((ispace.eq.1).and.(space(ix,iy,iz+1).ne.0)) goto 707
	  if (pspace(ispace,ix,iy,iz+1).ne.0) goto 707
	  if (ptest(ix,iy,iz+1).ne.0) goto 707
	   iz=iz+1
	   ptest(ix,iy,iz)=ip
	   isuccess=1
	   goto 700
	endif

	if (ir.eq.6) then
	  if ((ispace.eq.1).and.(space(ix,iy,iz-1).ne.0)) goto 707
	  if (pspace(ispace,ix,iy,iz-1).ne.0) goto 707
	  if (ptest(ix,iy,iz-1).ne.0) goto 707
	   iz=iz-1
	   ptest(ix,iy,iz)=ip
	   isuccess=1
	   goto 700
	endif

 707	continue
	enddo


 700	continue
	if (isuccess.eq.0) then
 720	format(a6,i3,i3,i5,4x,6i4)
	goto 701
	else

	plect(ip,istep+1,1)=ix
	plect(ip,istep+1,2)=iy
	plect(ip,istep+1,3)=iz
	if (istep.eq.nlength-1) goto 730
	endif

C istep segments in plectoneme
	enddo
 701	continue
C iptry redo failing plectoneme
	enddo
C ispace redo failing plectoneme
	enddo
 730	continue

	plectspace(ip)=ispace
	ispacemax=max(ispacemax,ispace)

	do istep=1,nlength
	pspace(plectspace(ip),plect(ip,istep,1),
     &     plect(ip,istep,2),plect(ip,istep,3))=999
	enddo

C ip plectonemes
	enddo
C ipass
	enddo

	write(6,*) "number of spaces ",ispacemax

	rs=0.27/float(ispacemax)
	
	chainout="A"
C---write out coordinates
	iatom=1
	iatomout=1
	icirclecount=0
	do ip=1,nplect
	nlength=plectlength(ip)

	if (ip.eq.1) then
	  istart=1
	  if (plectflag(ip).eq.1) iend=plectroot(ip)-1
	  if (plectflag(ip).gt.1) goto 777
	else
	  istart=plectroot(ip-1)+1
	  if (plectflag(ip).eq.1) iend=plectroot(ip)-1
	  if (plectflag(ip).gt.1) goto 777
	endif
	
	write(6,*) istart,iend

	rspace=float(plectspace(ip))
	rmax=float(ispacemax)+1.
	roffset=1./rmax*(rspace-1.)
c -- write circle segment
 40	format(a6,5i8)
	ires=0
C--these two if make the junctions match up on origin and on RNA
	do i=istart,iend

	do is=1,nseg
	if (segout(is,1).eq.iatom) write(2,40) "SEGMNT",is
	enddo

	rx=float(coord(i,1))
	ry=float(coord(i,2))
	rz=float(coord(i,3))
	ires=ires+1
C------ some stuff to label daughter strands correctly
	if (iatom.le.ncircle) then
	  iatomout=iatom
	else
	  chainout="B"
	  if (iatom.lt.ncircle+ifork1) then
	     iatomout=iatom-ncircle
	  else
	     iatomout=iatom-ncircle-ifork1+ifork2
	  endif
	endif
C---------------
	write(2,88) "ATOM",iatomout," CA ","GLY",chainout,
c    &    ires,rx,ry,rz,1.,10.
     &    ires,(rcoord(i,j),j=1,3),1.,10.
 88	format(a4,2x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2)
	iatom=iatom+1
	if (iatom.gt.99999) iatom=1
	icirclecount=icirclecount+1
	enddo
c -- write plectoneme
c -- check if it is RNA
	if (plectflag(ip).eq.2) plectrootout(ip)=iatom-1
	if (plectflag(ip).gt.1) goto 777
	do i=0,nlength
	do j=1,3
	pv(i,j)=0.
	pnormal(i,j)=0.
	pnext(i,j)=0.
	enddo
	enddo
c -- offset for first step --
	icircle=plectroot(ip)-1

	do j=1,3
	  pnext(0,j)=float(plect(ip,1,j)-plect(ip,0,1))
	  pnormal(0,j)=float(plect(ip,1,j)-plect(ip,0,1))
	  pv(0,j)=float(plect(ip,0,j))-rcoord(icircle,j)
	enddo

	if (pv(0,1)+pv(0,2)+pv(0,3).eq.0) then
	write(6,*) "zero"
	write(6,*) "pnext   ",(pnext(0,j),j=1,3)
	write(6,*) "pnormal ",(pnormal(0,j),j=1,3)
	write(6,*) "icircle ",icircle
	write(6,*) "plect0  ",(plect(ip,0,j),j=1,3)
	write(6,*) "coord   ",(coord(icircle,j),j=1,3)
	write(6,*) "plect0  ",(plect(ip,0,j),j=1,3),(plect(ip,1,j),j=1,3)
	write(6,*) "coord+2 ",(coord(icircle+2,j),j=1,3)
	write(6,*) "pv0     ",(pv(0,j),j=1,3)
	write(6,*)
	endif
c -- calculate offsets
	isharplast=0.
	do istep=1,nlength-1
	rl=0.
	rn=0.
	do j=1,3
	pnext(istep,j)=float(plect(ip,istep+1,j)-plect(ip,istep,j))
	pnormal(istep,j)=float(plect(ip,istep+1,j)-plect(ip,istep-1,j))
	rn=rn+(pnext(istep,j))**2
	rl=rl+(pnormal(istep,j))**2
	enddo
	if ((rn.eq.0.).or.(rl.eq.0.)) write(6,*) "rn rl zero",rn,rl
	do j=1,3
	pnext(istep,j)=pnext(istep,j)/sqrt(rn)
	pnormal(istep,j)=pnormal(istep,j)/sqrt(rl)
	enddo
	enddo

cwrite(6,*) ip,0,(pv(0,j),j=1,3)
	do is=1,nlength-1

	rparallel=(pnext(is-1,1)*pnext(is,1))+
     &            (pnext(is-1,2)*pnext(is,2))+
     &            (pnext(is-1,3)*pnext(is,3))
	if (rparallel.gt.0.01) goto 600

c --  sharp turn

	rdiff=(pv(is-1,1)*pnormal(is,1))+
     &        (pv(is-1,2)*pnormal(is,2))+
     &        (pv(is-1,3)*pnormal(is,3))
	if ((rdiff.gt.0.99).or.(rdiff.lt.-0.99)) then
	rxc=-pnext(is-1,2)*pnormal(is,3)+
     &       pnext(is-1,3)*pnormal(is,2)
	ryc=-pnext(is-1,3)*pnormal(is,1)+
     &       pnext(is-1,1)*pnormal(is,3)
	rzc=-pnext(is-1,1)*pnormal(is,2)+
     &       pnext(is-1,2)*pnormal(is,1)
	rn=sqrt(rxc*rxc+ryc*ryc+rzc*rzc)
	if (rn.eq.0) then
C --vectors are exactly parallel
	  write(6,*) "sharp turn 0 -- setting arbitrary normal "
	  if (pnext(is-1,1).ne.0.) then
	    ryc=1.
	  else
	    rxc=1.
	  endif
	  rn=1.
	endif
	rlast=1.
	if (isharplast.ne.0) rlast=-1.
	pv(is,1)=rlast*rxc/rn
	pv(is,2)=rlast*ryc/rn
	pv(is,3)=rlast*rzc/rn
	isharplast=1
	goto 601
	endif

 600	continue
	isharplast=0
	rxc=-pv(is-1,2)*pnormal(is,3)+pv(is-1,3)*pnormal(is,2)
	ryc=-pv(is-1,3)*pnormal(is,1)+pv(is-1,1)*pnormal(is,3)
	rzc=-pv(is-1,1)*pnormal(is,2)+pv(is-1,2)*pnormal(is,1)
	rn=sqrt(rxc*rxc+ryc*ryc+rzc*rzc)
	if (rn.eq.0) then
	   write(6,*) "zero"
	   write(6,*) (plect(ip,is-1,j),j=1,3)
	   write(6,*) (plect(ip,is,j),j=1,3)
	   write(6,*) (pv(is-1,j),j=1,3)
	   write(6,*) (pnormal(is,j),j=1,3)
	rn=1.
	endif
	pv(is,1)=rxc/rn
	pv(is,2)=ryc/rn
	pv(is,3)=rzc/rn

 601	continue

	enddo
	pv(nlength,1)=pv(nlength-1,1)
	pv(nlength,2)=pv(nlength-1,2)
	pv(nlength,3)=pv(nlength-1,3)
		 
c -- write first strand of plectoneme
c -- vectors are now twisted by 90 deg every step, change to superhelical density
	do i=1,nlength
	ra0=float(mod(i,4))*90.
	rai=mod(float(i)*rsupercoil,360.)
	ra=(rai-ra0)*pi/180.
	rcos=cos(ra)
	rsin=sin(ra)
	rx=pnormal(i,1)
	ry=pnormal(i,2)
	rz=pnormal(i,3)
	pvx=pv(i,1)*(rcos+rx*rx*(1.-rcos))+
     &      pv(i,2)*(rx*ry*(1.-rcos)-rz*rsin)+
     &      pv(i,3)*(rx*rz*(1.-rcos)+ry*rsin)
	pvy=pv(i,2)*(rcos+ry*ry*(1.-rcos))+
     &      pv(i,3)*(ry*rz*(1.-rcos)-rx*rsin)+
     &      pv(i,1)*(ry*rx*(1.-rcos)+rz*rsin)
	pvz=pv(i,3)*(rcos+rz*rz*(1.-rcos))+
     &      pv(i,1)*(rz*rx*(1.-rcos)-ry*rsin)+
     &      pv(i,2)*(rz*ry*(1.-rcos)+rx*rsin)
	pv(i,1)=pvx	
	pv(i,2)=pvy	
	pv(i,3)=pvz	

	enddo

	do i=0,nlength
C------ some stuff to label daughter strands correctly
	if (iatom.le.ncircle) then
	  iatomout=iatom
	else
	  chainout="B"
	  if (iatom.lt.ncircle+ifork1) then
	     iatomout=iatom-ncircle
	  else
	     iatomout=iatom-ncircle-ifork1+ifork2
	  endif
	endif
C---------------

	rox=pv(i,1)*rs
	roy=pv(i,2)*rs
	roz=pv(i,3)*rs
	rx=(float(plect(ip,i,1))+roffset-rox)
	ry=(float(plect(ip,i,2))+roffset-roy)
	rz=(float(plect(ip,i,3))+roffset-roz)
	ires=ires+1
	write(2,88) "ATOM",iatomout," CA ","ARG",chainout,
     &     ires,rx,ry,rz,1.,plectspace(ip)*10.
	iatom=iatom+1
	if (iatom.gt.99999) iatom=1
	enddo

	do i=nlength-1,0,-1
	rox=pv(i,1)*rs
	roy=pv(i,2)*rs
	roz=pv(i,3)*rs
	rx=(float(plect(ip,i,1))+roffset+rox)
	ry=(float(plect(ip,i,2))+roffset+roy)
	rz=(float(plect(ip,i,3))+roffset+roz)
	ires=ires+1
	write(2,88) "ATOM",iatomout+1," CA ","LYS",chainout,
     &    ires,rx,ry,rz,1.,plectspace(ip)*10.
	iatom=iatom+1
	iatomout=iatomout+1
	if (iatom.gt.99999) iatom=1
	enddo

 777	continue

	enddo
	

c---	write end segment of circle 
	write(6,*) "writing istartlast/iendtot ",istartlast,iendtot

	ires=0
C	fudge factor for cases with no plectonemes
	if (nplect.eq.0) iendtot=iendtot+1
	do i=istartlast,iendtot-1

	rx=float(coord(i,1))
	ry=float(coord(i,2))
	rz=float(coord(i,3))
	ires=ires+1
C------ some stuff to label daughter strands correctly
	if (iatom.le.ncircle) then
	  iatomout=iatom
	else
	  chainout="B"
	  if (iatom.lt.ncircle+ifork1) then
	     iatomout=iatom-ncircle
	  else
	     iatomout=iatom-ncircle-ifork1+ifork2
	  endif
	endif
C---------------
	do is=1,nseg
	if (segout(is,1).eq.iatom) write(2,40) "SEGMNT",is
	enddo
	write(2,88) "ATOM",iatomout," CA ","GLY",chainout,
c    &      ires,rx,ry,rz,1.,10.
     &      ires,(rcoord(i,j),j=1,3),1.,10.
	iatom=iatom+1
	if (iatom.gt.99999) iatom=1
	enddo


	write(10,40) "DNA   ",1,iatom-1,iatom-1
	
C-- write RNA
	iatom=1
	imod=nplect+2
	do ip=1,nplect
	if (plectflag(ip).ne.1) then
	write(2,40) "RNA   ",ip
	ires=0
	write(2,40) "ROOT  ",plectout(ip)
	do i=1,plectlength(ip)
	rx=(float(plect(ip,i,1))+roffset)
	ry=(float(plect(ip,i,2))+roffset)
	rz=(float(plect(ip,i,3))+roffset)
	ires=ires+1
	if (plectflag(ip).eq.2) then
	write(2,88) "ATOM",iatom," CA ","MET","T",
     &       ires,rx,ry,rz,1.,10.
	else
	write(2,88) "ATOM",iatom," CA ","CYS","F",
     &       ires,rx,ry,rz,1.,10.
	endif
	iatom=iatom+1
	if (iatom.gt.99999) iatom=1
	enddo
	write(10,40) "RNA   ",iatom-plectlength(ip),iatom-1,
c    &     plectrootout(ip)
     &     plectout(ip)
	imod=imod+1
	endif
	enddo
	end
