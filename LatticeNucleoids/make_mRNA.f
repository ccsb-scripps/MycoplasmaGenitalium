C make_mRNA.f
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
C reads csv file from unit 5 with positions of polymerases and nucleic acid binding proteins
C  (see file for format)
C writes a file to unit 6 with positions of superhelical plectonemes:
C DNA    58000                                    # DNA with 58000 beads
C FORK       0       0                            # not used (defines replication forks)
C PLEC      10     324                            # plectoneme start and end
C RNA1     334      49       1                    # nascent RNA: RNApol position, length, direction
C RNA2     304       3      43     219     263    # free RNA: length, number of ribosomes, positions of ribosomes
C
	integer*4 pol(0:1000),dir(0:1000),mrnalength(0:1000)
	integer*4 free(0:1000),rib(0:1000,0:100),nrib(0:1000)
	integer*4 nonplect(1000,2),iplectlength(1000)
	integer*4 plect(1000,2)
	integer*4 plectout(1000,2)
	character*20 label
	character*200 line

	ncircle=580000
	nplectmin=40

	iseed=1234
	do i=1,iseed
	  r=rand()
	enddo

	do i=0,1000
	nrib(i)=0
	do j=0,100
	rib(i,j)=0
	enddo
	enddo

C will be used when replication forks are implemented:
	ifork1=0
	ifork2=0

	rjitter=0.
	iplectmax=800

	npol=1
	nfree=0
 1	read(5,9,end=99) line
 9	format(a200)

	if (line(1:7).eq."RNA_POL") then
 2	read(line,*) label,pol(npol),dir(npol),mrnalength(npol)
	write(3,*) label, pol(npol),dir(npol),mrnalength(npol)
	if (mrnalength(npol).le.0) goto 1
	npol=npol+1
	endif

	if (line(1:4).eq."mRNA") then
  	read(line,*) label,index,label,ifree,icopy
	free(index)=ifree
	nfree=max(index,nfree)
	endif

	if (line(1:4).eq."RIBO") then
	read(line,*) label,irib,label,ilength,imrna,iposition
	if (iposition.eq.0) iposition=20
	nrib(imrna)=nrib(imrna)+1
	rib(imrna,nrib(imrna))=iposition
	endif

	goto 1
 99	continue
	npol=npol-1

C--everything read in now----

c add dummy RNA_POL so plectonomes are excluded at ends of circle segments
        npol=npol+1
        pol(npol)=1
        dir(npol)=-98
        npol=npol+1
        pol(npol)=ncircle
        dir(npol)=-98

C sort RNApol
	do ipol=1,npol
	do ipol2=ipol+1,npol
	if (pol(ipol).gt.pol(ipol2)) then
	  poltmp=pol(ipol)
	  dirtmp=dir(ipol)
	  mrnatmp=mrnalength(ipol)
	  pol(ipol)=pol(ipol2)
	  dir(ipol)=dir(ipol2)
	  mrnalength(ipol)=mrnalength(ipol2)
	  pol(ipol2)=poltmp
	  dir(ipol2)=dirtmp
	  mrnalength(ipol2)=mrnatmp
	endif
	enddo
	enddo

	do ipol=1,npol
	write(3,*)  pol(ipol),dir(ipol),mrnalength(ipol)
	enddo

C	write transcript file
	write(6,100) "DNA ",58000

C hardwired now to do segment break near RNA polymerase

C will be used when replication forks are implemented:
	write(6,100) "FORK",ifork1/10,ifork2/10

	plect(1,1)=10
	do ipol=1,npol-1

	if (pol(ipol+1)-pol(ipol).gt.100) then

	istart=pol(ipol)/10
	iend=pol(ipol+1)/10

	nplect=max(1,(iend-istart)/iplectmax)
	istep=(iend-istart)/nplect

	do iplect=1,nplect
	 plectout(iplect,1)=istart+(iplect-1)*istep+
     &       int(rand()*rjitter)+10
	 plectout(iplect,2)=istart+(iplect)*istep-
     &       int(rand()*rjitter)-10
	iplectlength(iplect)=plectout(iplect,2)-plectout(iplect,1)
	if ((plectout(iplect,1).lt.plectout(iplect,2)).and.
     &     (iplectlength(iplect).gt.nplectmin)) then
	  write(6,100) "PLEC",(plectout(iplect,j),j=1,2)
	endif
	enddo

	endif
	enddo
 100	format(a4,22i8)

	do ipol=1,npol
	if ((mrnalength(ipol).ne.0).and.(dir(ipol).gt.0)) then
	 mrnaout=max(6,mrnalength(ipol)/10)
	 write(6,100) "RNA1",pol(ipol)/10,mrnaout,dir(ipol)
	endif
	enddo
	
	do ifree=0,nfree
	if (free(ifree).ne.0) then
	 write(6,100) "RNA2",free(ifree)/10,nrib(ifree),
     &      (rib(ifree,j)/10,j=1,nrib(ifree))
	endif
	enddo
	
	end
