C make_segments.f
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
C reads transcript file from unit 5 with information on nucleic acids
C writes a segment file to unit 6 with segments of nucleic acid that need to be modeled:
C
C CIRC       1   58000       1    1886                    $ non-plectonemic circle: full genome start and end numbers,
C                                                            start and end number of beads in the circle
C PLEC      10     324      10     157                    $ plectoneme: start and end numbers,
C                                                            bead number in circle where it is inserted,
C                                                            length of plectoneme superhelix (==1/2 number of beads)
C RNA1   58001   58049      20      49 T                  $ nascent RNA: start and end numbers,
C                                                            bead number in circle where it is attached,
C                                                            length of RNA, label
C RNA2   69905   70208 R       3   69947   70123   70167  $ free RNA: start and end numbers,
C                                                            label, number of ribosomes, bead number of ribosomes
C
	integer*4 plect(5000,2),plectcircle(5000),plectlength(5000)
	integer*4 iptot(0:5000)
	integer*4 seg(10,2),segcircle(10,2)
	integer*4 rna1root(5000),rna1length(5000),rna1circle(5000)
	integer*4 rna2length(5000),rna2rib(5000),nrib(5000)
        integer*4 rib(5000,0:100)
	character*10 segspec(10)
	character*4 label
	integer*4 gene(5000,2),fork1,fork2,ndna
	character*200 line

 	ip=0
	ig=0
	irna1=0
	irna2=0
	itotal=0
C for future development of replication forks:
	fork1=0
	fork2=0
	nfork=0

 1	read(5,10,end=9) line
 10	format(a200)

	if (line(1:3).eq."DNA") then
	   read(line,11) ndna
	nseg=1
	endif
 11	format(4x,22i8)

	if (line(1:4).eq."PLEC") then
	   ip=ip+1
	   read(line,11) plect(ip,1),plect(ip,2)
	endif

	if (line(1:4).eq."GENE") then
	   ig=ig+1
	   read(line,11) gene(ig,1),gene(ig,2)
	endif

	if (line(1:4).eq."RNA1") then
	   irna1=irna1+1
	   read(line,11) rna1root(irna1),rna1length(irna1)
	endif

	if (line(1:4).eq."RNA2") then
	   irna2=irna2+1
	   read(line,11) rna2length(irna2),nrib(irna2)
	if (nrib(irna2).ne.0) then
	   read(line,11) rna2length(irna2),nrib(irna2),
     &     (rib(irna2,j),j=1,nrib(irna2))
	endif
	endif

	goto 1
 9	continue
	np=ip
	ng=ig
	nrna1=irna1
	nrna2=irna2

C segments 
	seg(1,1)=1
	seg(1,2)=ndna

C loop through plectonemes and find circle numbers
	iptot(0)=0
	do ip=1,np
	plectcircle(ip)=plect(ip,1)-iptot(ip-1)
	plectlength(ip)=plect(ip,2)-plect(ip,1)+1
c make sure length is divisible by 2
	plectlength(ip)=plectlength(ip)/2
	plectlength(ip)=plectlength(ip)*2
	iptot(ip)=iptot(ip-1)+plectlength(ip)

 12	format(a4,4i8)
	enddo

	do iseg=1,nseg

	segcircle(iseg,1)=seg(iseg,1)
	segcircle(iseg,2)=seg(iseg,2)
	do ip=1,np
	if (plect(ip,2).lt.seg(iseg,1))
     &      segcircle(iseg,1)=segcircle(iseg,1)-plectlength(ip)
	if (plect(ip,2).lt.seg(iseg,2))
     &      segcircle(iseg,2)=segcircle(iseg,2)-plectlength(ip)
	enddo

 31	continue

 32	continue

	enddo

	write(6,15) "CIRC",1,ndna,
     &       1,segcircle(1,2)

 15	format(a4,4i8,a10)

C ---write plectonemes
	do ip=1,np
	write(6,12) "PLEC",plect(ip,1),plect(ip,2),plectcircle(ip),
     &     plectlength(ip)/2
	enddo

C -- transcribing RNA
C --add "GENE" entries
	itotal=ndna
	if (fork1.ne.0) itotal=itotal+nfork
	do ig=1,ng
	nrna1=nrna1+1
	rna1length(nrna1)=(gene(ig,2)-gene(ig,1))/2+1
	rna1root(nrna1)=gene(ig,1)+rna1length(nrna1)
	enddo

	do irna=1,nrna1
	irnacircle=rna1root(irna)
	do ip=np,0,-1
	if (plect(ip,2).lt.rna1root(irna)) then
	  irnacircle=rna1root(irna)-iptot(ip)
	  goto 33
	endif
	enddo
 33	continue
	write(6,16) "RNA1",itotal+1,itotal+rna1length(irna),
     &    irnacircle,rna1length(irna),"T"
	itotal=itotal+rna1length(irna)
 16	format(a4,4i8,1x,a1)
	enddo
	
C -- free RNA

	do irna=1,nrna2
	write(6,17) "RNA2",itotal+1,itotal+rna2length(irna),
     &     "R",
     &     nrib(irna),(itotal+rib(irna,j),j=1,nrib(irna))
	itotal=itotal+rna2length(irna)
	enddo
 17	format(a4,2i8,1x,a1,22i8)

	end
