###############################################################################
### Function to count the number of lines in a file
def FileLength(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

############################################################################
### Write big.in or small.in file
def WriteObjInFile(here,whichdir,names,filename,Header,FirstLines,xv,s):
	import MercModule

	infile=open(here+'/'+whichdir+'/In/'+filename+'.in','w')
	# Header
	for i in range(len(Header)):
		infile.write(Header[i])
	# Data
	for i in range(len(names)):
		infile.write(FirstLines[i])
		infile.write("  "+xv[i][0]+"  "+xv[i][1]+"  "+xv[i][2]+"\n")
		infile.write("  "+xv[i][3]+"  "+xv[i][4]+"  "+xv[i][5]+"\n")
		infile.write(s[i])
	infile.close()

############################################################################
### Read info.out, put data into name/destination/time vectors and return
def ReadInfo(whichdir):
	assert type(whichdir) is str

### Load modules
	import MercModule
	import numpy as np

 	InfoFile=open(whichdir+'/Out/info.out','r')
	InfoLen=MercModule.FileLength(whichdir+'/Out/info.out')
	dest=	list(range((InfoLen-5)/4))
	time=	np.zeros(((InfoLen-5)/4))
	skip=[]
	flen=7
	header=True
	footer=False
	j=0

### Read through the file until reaching the start of integration
	while header==True:		#header, not needed
		j=j+1
		line=InfoFile.readline()
		if line=="   Beginning the main integration.\n":
			header=False
			j=j+1
			line=InfoFile.readline()
	hlen=j+1
	name1=	['']*(InfoLen-hlen-flen)
	dest1=	['']*(InfoLen-hlen-flen)
	time1=	[0]*(InfoLen-hlen-flen)

### Read through the integration data until reaching the file footer, 
###	parsing each line based on # of words to get collision data.
	while footer==False:	#body of file
		j=j+1
		line=InfoFile.readline()
		if line=="   Integration complete.\n":
			footer=True
			break
		splitline=line.split()
		if len(splitline)==8 and splitline[0]!='Continuing':
			name1[j-hlen],dest1[j-hlen],time1[j-hlen]=splitline[4],splitline[0],splitline[6]
		elif len(splitline)==9:
			name1[j-hlen],dest1[j-hlen],time1[j-hlen]=splitline[0], 'Sun',splitline[7]
		elif len(splitline)==5 and splitline[0]!='Fractional':
			name1[j-hlen],dest1[j-hlen],time1[j-hlen]=splitline[0],'Ejected',splitline[3]
		else:
			skip.append(splitline)
	
	InfoFile.close()
	return name1, dest1, time1

###############################################################################
# A program to read the collision info from info.out and add it to a running total
def CopyInfo(whichdir,whichtime,writegood):
	assert type(whichdir) is str
	assert type(whichtime) is str
	assert type(writegood) is bool

### Load modules
	import numpy as np
	import os
	import MercModule
#	from random import random
#	from math import pi, sin, cos
#	test change
	
	print('CopyInfo '+whichdir+'/Out/info.out, '+whichtime)

### Get collision info from info.out
	name1,dest1,time1 = MercModule.ReadInfo(whichdir)
	
### Summarize impacts for infosum.out
#	destname=['Sun','Mercury','Venus','Earth','Mars', 'Jupiter','Io','Europa','Ganymede','Callisto', 'Saturn','Enceladu','Rhea','Titan','Iapetus', 'Ejected']
	destname=['Sun','Mercury','Venus','Earth','Mars', 'Jupiter','Moon', 'Saturn','Ejected']
	InfoSum=[0]*(len(destname))
	if np.array(dest1) != 0:
		for k in range(len(destname)):
			InfoSum[k]=sum( np.array(dest1)==destname[k] )
#		InfoSum[13]=sum( np.array(dest1)=='ejected' )

### Get which timestep was used, and store in last column of InfoSum
	timestepfile=open(whichdir+'/timestep.txt','r')
	timestep=int(timestepfile.readline())
	timestepfile.close()

### Get total number of objects in simulation
	SmallInFileLength=MercModule.FileLength(whichdir+'/In/small.in')
	NTot=(SmallInFileLength-5)/4
	assert type(NTot) is int

### Write summed impacts to file
	InfoSumFile=open(whichdir+'/infosum.out','a')
	if os.path.getsize(whichdir+'/infosum.out')==0:
#		InfoSumFile.write('  Su  Me  Ve  Ea  Ma  Ju  Io  Eu  Ga  Ca  Sa  En  Rh  Ti  Ia  Ej Step\n')
		InfoSumFile.write('  Su  Me  Ve  Ea  Ma  Ju  Mn  Sa  Ej  Tot Step\n')
#	InfoSumFile.write(' {0:3d} {1:3d} {2:3d} {3:3d} {4:3d} {5:3d} {6:3d} {7:3d} {8:3d} {9:3d} {10:3d} {11:3d} {12:3d} {13:3d} {14:3d} {15:3d} {16:4d}\n'.format( \
#	InfoSum[0], InfoSum[1], InfoSum[2], InfoSum[3], InfoSum[4], InfoSum[5], InfoSum[6], InfoSum[7], InfoSum[8], InfoSum[9], InfoSum[10], InfoSum[11], InfoSum[12], InfoSum[13], InfoSum[14], InfoSum[15], InfoSum[16]))
	InfoSumFile.write(' {0:3d} {1:3d} {2:3d} {3:3d} {4:3d} {5:3d} {6:3d} {7:3d} {8:3d} {9:4d} {10:4d}\n'.format( \
	InfoSum[0], InfoSum[1], InfoSum[2], InfoSum[3], InfoSum[4], InfoSum[5], InfoSum[6], InfoSum[7], InfoSum[8], NTot, timestep))
	InfoSumFile.close()

### Get .in data for rocks that hit something and write to file
	if writegood==True:
		gooddest=['Jupiter','Io','Europa','Ganymede','Callisto','Moon', 'Saturn','Enceladu','Rhea','Titan','Iapetus']
		ind=np.array([
			(any(dest1[i]==gooddest[j] for j in range(len(gooddest)))) 
			for i in range(len(dest1)) ])
		if (len(name1) > 0):
			name=np.array(name1)[ind]
			dest=np.array(dest1)[ind]
			#print(name[0],dest[0],time1[0])
			goodin=open(whichdir+'/good.in','a')
			smallin=open(whichdir+'/In/small.in','r')
			SmallLen=MercModule.FileLength(whichdir+'/In/small.in')
			smalllines=['' for i in range(SmallLen)]
			for j in range(5):
				smalllines[j]= smallin.readline()
			for j in range(5,SmallLen):
				smalllines[j]= smallin.readline()
				if any(name[i]==smalllines[k].split()[0] and float(time1[i])>60. for i in range(len(name)) for k in range(j-3,j+1)):
					#for i in range(len(name)):
					#print(smalllines[j])
					#	for k in range(j-3,j+1):
					#		print(smalllines[k])
					#		print(name[i])
					#		print(time[i])
					goodin.write(smalllines[j])
			goodin.close()
			smallin.close()


############################################################################
### Select a timestep for the big objects and make a new big.in
def MakeMoon(whichdir,whichtime):
	assert type(whichdir) is str
	assert type(whichtime) is str
	assert int(whichtime) >= 5

	# Needed modules
	import MercModule
	import numpy as np
	import os
	from random import random
	from math import sqrt, pi, sin, cos
	
	here=os.getcwd()
	
	print('MakeMoon '+whichdir+'/In/big.in  '+whichtime)
	
	#constants/variables
	G = 6.674e-8 # cm^3/g/s^2
	mSun = 1.99e33					#g
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
	
#	big=['Mercury','Venus','Earth','Mars','Jupiter',
#	'Io','Europa','Ganymede','Callisto','Saturn','Uranus','Neptune']
	big=['Mercury','Venus','Earth','Mars','Jupiter','Moon',
	'Saturn','Uranus','Neptune']
	bigxv=['' for i in range(len(big))]

### Use chosen timestep
	timestep=int(whichtime)

# Find the correct timestep for the planets
	for i in range(5)+[6,7,8]:
		filename=here+'/'+whichdir+'/In/InitElemFiles/'+big[i]+'.aei'
		File=open(filename,'r')
		for j in range(timestep):
			thisline=File.readline().split()
		bigxv[i]=thisline[6:]
### Moon parameters based on which run
### HDir => smaller mass, #1 => smaller axis, currently all small density
	FakeFile=open('FakeMoons.txt','r')
	Fake=FakeFile.readlines()
	FakeFile.close()
	iFake=1
	if whichdir[0]=='H':
		jFake=1
	elif whichdir[0]=='L':
		jFake=2
	kFake=int(whichdir[-1])
	d=Fake[0].split()[iFake]
	m=str(float(Fake[1].split()[jFake])/mSun)
	a=Fake[2].split()[kFake]

### Generate Moon's position and velocity
	phi1=2*pi*random()				# planar angle for position
	theta1=0.0						# polar angle for position (i=0 => 0)
	v=sqrt(G*9.54266E-04*mSun/float(a))	# orbital velocity = sqrt(GM/a)
	phi2=phi1+pi/2					# planar angle for velocity
	theta2=0.0						# polar angle for velocity (i=0 => 0)

	x=float(a)*cos(phi1)*cos(theta1)/AU	# moon orbit relative to planet
	y=float(a)*sin(phi1)*cos(theta1)/AU
	z=float(a)*sin(theta1)          /AU
	u=v*cos(phi2)*cos(theta2)*day/AU
	v=v*sin(phi2)*cos(theta2)*day/AU
	w=v*sin(theta2)			 *day/AU
	xvmod=[x,y,z, u,v,w]

### Add Moon's position to Jupiter's
	bigxv[5]=[repr(float(bigxv[4][i])+xvmod[i]) for i in range(6)]

### Read density and mass of planets
	BigFirstLineFile=open('PlanetFirstLines.txt','r')
	BigFirstData=BigFirstLineFile.readlines()
	BigFirstData=[BigFirstData[i].split() for i in range(len(BigFirstData))]
	BigFirstData=np.array(BigFirstData)
	dpl=np.array([0. for i in range(len(big))])
	mpl=np.array([0. for i in range(len(big))])

	for j in range(len(big)):
		if np.any(BigFirstData[:,0]==big[j]):
			dpl[j]=BigFirstData[BigFirstData[:,0]==big[j],1][0]
			mpl[j]=BigFirstData[BigFirstData[:,0]==big[j],2][0]
	dpl[np.array(big)=='Moon']=d
	mpl[np.array(big)=='Moon']=m

	dpl=[str(dpl[i]) for i in range(len(dpl))]
	mpl=[str(mpl[i]) for i in range(len(dpl))]

### Format data as the first line of each big.in object entry
	BigFirstLines=[big[i].ljust(10)+'d= '+dpl[i]+'  m= '+mpl[i]+'\n'
                   for i in range(len(big))]

### Read generic big.in file header	
	BigHeadFile=open('bigheader.txt','r')
	BigHeader=BigHeadFile.readlines()
	BigHeadFile.close()

### No spin for all objects
	bigs=["  0.0  0.0  0.0\n" for i in range(len(BigFirstLines))]

	MercModule.WriteObjInFile(here,whichdir,big,'big',
                            BigHeader,BigFirstLines,bigxv,bigs)

### Record which timestep was used
	timestepfile=open(here+'/'+whichdir+'/timestep.txt','w')
	timestepfile.write(repr(timestep))
	timestepfile.close()

############################################################################
### Select a timestep for the big objects and make a new big.in
def MakeBigChoose(whichdir,whichtime):
	assert type(whichdir) is str
	assert type(whichtime) is str
	assert int(whichtime) >= 5

	# Needed modules
	import MercModule
	import numpy as np
	import os
#	from random import random
#	from math import pi, sin, cos
	
	here=os.getcwd()
	
	print('BakeBigChoose '+whichdir+'/In/big.in  '+whichtime)
	
	#constants/variables
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
	
#	big=['Mercury','Venus','Earth','Mars','Jupiter',
#	'Io','Europa','Ganymede','Callisto','Saturn','Uranus','Neptune']
	big=['Mercury','Venus','Earth','Mars', 		'Jupiter','Io','Europa','Ganymede','Callisto',
	'Saturn','Enceladu','Rhea','Titan','Iapetus','Uranus','Neptune']
	bigxv=['' for i in range(len(big))]

### Use chosen timestep
	timestep=int(whichtime)

# Find the correct timestep for each big thing
	for i in range(len(big)):
		filename=here+'/'+whichdir+'/In/InitElemFiles/'+big[i]+'.aei'
		File=open(filename,'r')
		for j in range(timestep):
			thisline=File.readline().split()
		bigxv[i]=thisline[6:]

### Read density and mass of planets
	BigFirstLineFile=open('PlanetFirstLines.txt','r')
	BigFirstData=BigFirstLineFile.readlines()
	BigFirstData=[BigFirstData[i].split() for i in range(len(BigFirstData))]
	BigFirstData=np.array(BigFirstData)
	dpl=np.array([0. for i in range(len(big))])
	mpl=np.array([0. for i in range(len(big))])

	for j in range(len(big)):
		if np.any(BigFirstData[:,0]==big[j]):
			dpl[j]=BigFirstData[BigFirstData[:,0]==big[j],1][0]
			mpl[j]=BigFirstData[BigFirstData[:,0]==big[j],2][0]

	dpl=[str(dpl[i]) for i in range(len(dpl))]
	mpl=[str(mpl[i]) for i in range(len(dpl))]

### Format data as the first line of each big.in object entry
	BigFirstLines=[big[i].ljust(10)+'d= '+dpl[i]+'  m= '+mpl[i]+'\n'
                   for i in range(len(big))]

### Read generic big.in file header	
	BigHeadFile=open('bigheader.txt','r')
	BigHeader=BigHeadFile.readlines()
	BigHeadFile.close()

### No spin for all objects
	bigs=["  0.0  0.0  0.0\n" for i in range(len(BigFirstLines))]

	MercModule.WriteObjInFile(here,whichdir,big,'big',
                            BigHeader,BigFirstLines,bigxv,bigs)

### Record timestep used
	timestepfile=open(here+'/'+whichdir+'/timestep.txt','w')
	timestepfile.write(repr(timestep))
	timestepfile.close()


############################################################################
### Pick a random timestep for the big objects and make a new big.in
def MakeBigRand(whichdir,whichtime):
	assert type(whichdir) is str
	assert type(whichtime) is str

	# Needed modules
	import MercModule
	import numpy as np
	import os
	from random import random
#	from math import pi, sin, cos
	
	here=os.getcwd()
	
	print('MakeBigRand '+whichdir+'/In/big.in  '+whichtime)
	
	#constants/variables
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
	
#	big=['Mercury','Venus','Earth','Mars','Jupiter',
#	'Io','Europa','Ganymede','Callisto','Saturn','Uranus','Neptune']
	big=['Mercury','Venus','Earth','Mars', 		'Jupiter','Io','Europa','Ganymede','Callisto',
	'Saturn','Enceladu','Rhea','Titan','Iapetus','Uranus','Neptune']
	bigxv=['' for i in range(len(big))]

### Pick a random timestep and get all big vectors at that point
	AEILen=MercModule.FileLength(here+'/'+whichdir+'/In/InitElemFiles/Jupiter.aei')-5
	timestep=5+int(AEILen*random())

# Find the correct timestep for each big thing
	for i in range(len(big)):
		filename=here+'/'+whichdir+'/In/InitElemFiles/'+big[i]+'.aei'
		File=open(filename,'r')
		for j in range(timestep):
			thisline=File.readline().split()
		bigxv[i]=thisline[6:]

### Read density and mass of planets
	BigFirstLineFile=open('PlanetFirstLines.txt','r')
	BigFirstData=BigFirstLineFile.readlines()
	BigFirstData=[BigFirstData[i].split() for i in range(len(BigFirstData))]
	BigFirstData=np.array(BigFirstData)
	dpl=np.array([0. for i in range(len(big))])
	mpl=np.array([0. for i in range(len(big))])

	for j in range(len(big)):
		if np.any(BigFirstData[:,0]==big[j]):
			dpl[j]=BigFirstData[BigFirstData[:,0]==big[j],1][0]
			mpl[j]=BigFirstData[BigFirstData[:,0]==big[j],2][0]

	dpl=[str(dpl[i]) for i in range(len(dpl))]
	mpl=[str(mpl[i]) for i in range(len(dpl))]

### Format data as the first line of each big.in object entry
	BigFirstLines=[big[i].ljust(10)+'d= '+dpl[i]+'  m= '+mpl[i]+'\n'
                   for i in range(len(big))]

### Read generic big.in file header	
	BigHeadFile=open('bigheader.txt','r')
	BigHeader=BigHeadFile.readlines()
	BigHeadFile.close()

### No spin for all objects
	bigs=["  0.0  0.0  0.0\n" for i in range(len(BigFirstLines))]

### Write data
	MercModule.WriteObjInFile(here,whichdir,big,'big',
                            BigHeader,BigFirstLines,bigxv,bigs)

### Record timestep used
	timestepfile=open(here+'/'+whichdir+'/timestep.txt','w')
	timestepfile.write(repr(timestep))
	timestepfile.close()

############################################################################
### Make a random cluster of small objects around preselected ones
### to write to small.in
def MakeSmall(whichdir,whichtime,n,whichpl,da,dv):
	assert type(whichdir) is str
	assert type(whichtime) is str
	assert type(n) is int
	assert n>=0

### Needed modules
	import MercModule
	import numpy as np
	import os
	from random import random
	from math   import pi, sin, cos

	here=os.getcwd()
	
	print('MakeSmall '+whichdir+'/In/small.in  '+whichtime+'  '+str(n))

### Constants/variables
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
	maxaspread=da/AU				#in AU
	maxvspread=dv*day/AU			#in AU/day
	
	small=['M'+str(i) for i in range(n)]
	smallxv=['' for i in range(n)]

### Generate slightly randomized rocks at different phases of Jupiter or Saturn's orbit
	for j in range(0,len(small)):
		### Pick a random timestep
		if whichpl=='J':
			filename=here+'/'+whichdir+'/In/InitElemFiles/Jupiter12Yr.aei'
		if whichpl=='S':
			filename=here+'/'+whichdir+'/In/InitElemFiles/Saturn29Yr.aei'
		AEILen=MercModule.FileLength(filename)-5
		timestep=5+int(AEILen*random())
		### Get Jupiter/Saturn's info at this point
		File=open(filename,'r')
		for k in range(timestep):
			thisline=File.readline().split()
		bigxv=thisline[6:]

		### Generate random variation
		phi1=2*pi*random()
		theta1=-pi/2+pi*random()
		r=maxaspread*random()
		v=maxvspread*random()
		phi2=2*pi*random()
		theta2=-pi/2+pi*random()

		x=r*cos(phi1)*cos(theta1)
		y=r*sin(phi1)*cos(theta1)
		z=r*sin(theta1)
		u=v*cos(phi2)*cos(theta2)
		v=v*sin(phi2)*cos(theta2)
		w=v*sin(theta2)

		### Coords = Jupiter/Saturn coords plus random variation
		smallxv[j]=[float(bigxv[0])+x, float(bigxv[1])+y, float(bigxv[2])+z,
		float(bigxv[3])+u, float(bigxv[4])+v, float(bigxv[5])+w]
		smallxv[j]=[repr(i) for i in smallxv[j]]

### Read density and mass of planetesimals
	SmallFirstLineFile=open('PlanetFirstLines.txt','r')
	SmallFirstData=SmallFirstLineFile.readlines()
	SmallFirstData=[SmallFirstData[i].split() 
                    for i in range(len(SmallFirstData))]
	SmallFirstData=np.array(SmallFirstData)
	d=str(SmallFirstData[SmallFirstData[:,0]=='Plantsml',1][0])
	m=str(SmallFirstData[SmallFirstData[:,0]=='Plantsml',2][0])

### Format data as the first line of each big.in object entry
	SmallFirstLines=[small[i].ljust(10)+'d= '+d+'  m= '+m+'  r= 0.001\n'
                   for i in range(len(small))]

### Read generic big.in file header	
	SmallHeadFile=open('SmallHeader.txt','r')
	SmallHeader=SmallHeadFile.readlines()
	SmallHeadFile.close()

### No spin for all objects
	smalls=["  0.0  0.0  0.0\n" for i in range(len(small))]

### Write data
	MercModule.WriteObjInFile(here,whichdir,small,'small',
                            SmallHeader,SmallFirstLines,smallxv,smalls)

############################################################################
### Copy the small rocks from good.in to small.in
def Good2Small(whichdir,whichtime,n):
	assert type(whichdir) is str
	assert type(whichtime) is str
	assert type(n) is int
	assert n>=0

### Needed modules
	import MercModule
	import numpy as np
	import os
	from random import random, sample
	from math import pi, sin, cos

	here=os.getcwd()	
	print('Good2Small '+whichdir+'/good.in  '+whichtime)

	#constants/variables
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
	
### Read in objects from good.in file
#	goodin=open(whichdir+'/good.in','r')
#	GoodLen=MercModule.FileLength(whichdir+'/good.in')
	goodin=open('good.in','r')
	GoodLen=MercModule.FileLength('good.in')
	if (GoodLen%4 != 0):
		print('??? good.in length = '+str(GoodLen))
	ngood=GoodLen/4
### If asked for more objects than available, give all and print warning
	if ngood < n:
		print('Warning: requested '+str(n)+' small objects. '+
		      'There are only '+str(ngood)+' good objects to use.')
		n=ngood
### New objects will be a random subset of the good list
	GoodInd=sample(xrange(ngood),n)
### Create and fill vectors with the data from good.in
	header, pos, vel, s = [],[],[],[]
	for j in range(ngood):
		header.append(goodin.readline() )		# not used
		pos.append(goodin.readline().split() )
		vel.append(goodin.readline().split() )
		s.append(goodin.readline() )

### Generate new, sequential names for the objects
	name=['M'+str(i) for i in range(n)]
	smallxv=['' for i in range(n)]
	smalls =['' for i in range(n)]

### Fill data of object j in new list with ind[j] from old list
	for j in range(n):
		smallxv[j]=pos[GoodInd[j]]+vel[GoodInd[j]]
		smalls[j] =s[GoodInd[j]]

### Read density and mass of planetesimals
	SmallFirstLineFile=open('PlanetFirstLines.txt','r')
	SmallFirstData=SmallFirstLineFile.readlines()
	SmallFirstData=[SmallFirstData[i].split() 
                    for i in range(len(SmallFirstData))]
	SmallFirstData=np.array(SmallFirstData)
	d=str(SmallFirstData[SmallFirstData[:,0]=='Plantsml',1][0])
	m=str(SmallFirstData[SmallFirstData[:,0]=='Plantsml',2][0])

### Format data as the first line of each big.in object entry
	SmallFirstLines=[name[i].ljust(10)+'d= '+d+'  m= '+m+'  r= 0.001\n'
                   for i in range(len(name))]

### Read generic big.in file header	
	SmallHeadFile=open('SmallHeader.txt','r')
	SmallHeader=SmallHeadFile.readlines()
	SmallHeadFile.close()

### Write data
	MercModule.WriteObjInFile(here,whichdir,name,'small',
                            SmallHeader,SmallFirstLines,smallxv,smalls)









