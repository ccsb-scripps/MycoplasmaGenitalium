C circle_segments.f
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
C Creates an idealized circular collection of closely-spaced points and a constraint file
C
C Reads from unit 5:
C MG149.segments           # input segment file
C MG149_circle.1.pdb       # output coordinate file
C MG149_circle.1.const     # output constraint file
C 15.,5.,1234              # max waviness, average number of waves, random seed
C
	character*80 filename,line

C	input segment file
	read(5,1) filename
	open(1,file=filename)

C	output coordinate file
	read(5,1) filename
	open(3,file=filename)

C	output constraint file
	read(5,1) filename
	open(2,file=filename)

 1	format(a80)
	read(5,*) rmax,rwave,iseed
	do i=1,iseed
	r=rand()
	enddo
	rnode1=4.+float(int(rand()*3.))
	rnode2=5.+float(int(rand()*4.))
	write(6,*) "nodes ",rnode1,rnode2

 20	read(1,1,end=9) line
	if (line(1:4).eq."CIRC") then
	  read(line,77) istart,npt
	endif
 77	format(4x,16x,2i8)
	goto 20
 9	continue

C	write constraint file
	write(2,10) "DNA   ",1,npt
	write(2,10) "CONSTR",1,npt,1.,1.,1,"ORI"
 10	format(a6,2i8,2f8.3,i2,1x,a3)
	
	icount=1
	rpi=3.14159
	do i=1,npt
	rang=float(i)/float(npt)*2.*rpi
	rang2=rang*rnode1
	rang3=rang*rnode2
	rxy=rmax+rwave*sin(rang2)
	rz=rwave*cos(rang3)
	xout=rxy*sin(rang)
	yout=rxy*cos(rang)
	zout=rz

	write(3,88) "ATOM",icount," CA ","PRO","A",
     &     icount,xout,yout,zout,1.,10.
88	format(a4,2x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2)
	icount=icount+1
	enddo

	end
