###############################################################################
### Function to count the number of lines in a file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

###############################################################################
# A program to read the collision info from info.out and add it to a running total
def copyinfo(whichdir,whichtime,writegood):
	import numpy, os
	import MercModule
#	from random import random
#	from math import pi, sin, cos
#	test change
	
	print('copyinfo '+whichdir+'/Out/info.out, '+whichtime)

 	InfoFile=open(whichdir+'/Out/info.out','r')
	InfoLen=MercModule.file_len(whichdir+'/Out/info.out')
	dest=	list(range((InfoLen-5)/4))
	time=	numpy.zeros(((InfoLen-5)/4))
	#name=['M'+str(i) for i in range(50)]
	skip=[]
	flen=7
	header=True
	footer=False
	j=0
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
	
	### Summarize impacts
	destname=['Sun','Mercury','Venus','Earth','Mars', 'Jupiter','Io','Europa','Ganymede','Callisto', 'Saturn','Enceladu','Rhea','Titan','Iapetus', 'Ejected']
	InfoSum=[0]*(len(destname)+1)
	if numpy.array(dest1) != 0:
		for k in range(len(destname)):
			InfoSum[k]=sum( numpy.array(dest1)==destname[k] )
#		InfoSum[13]=sum( numpy.array(dest1)=='ejected' )
	timestepfile=open(whichdir+'/timestep.txt','r')
	InfoSum[len(destname)]=int(timestepfile.readline())
	timestepfile.close()

	### Write summed impacts to file
	InfoSumFile=open(whichdir+'/infosum.out','a')
	if os.path.getsize(whichdir+'/infosum.out')==0:
		InfoSumFile.write('  Su  Me  Ve  Ea  Ma  Ju  Io  Eu  Ga  Ca  Sa  En  Rh  Ti  Ia  Ej  Stp\n')
	InfoSumFile.write(' {0:3d} {1:3d} {2:3d} {3:3d} {4:3d} {5:3d} {6:3d} {7:3d} {8:3d} {9:3d} {10:3d} {11:3d} {12:3d} {13:3d} {14:3d} {15:3d} {16:4d}\n'.format( \
	InfoSum[0], InfoSum[1], InfoSum[2], InfoSum[3], InfoSum[4], InfoSum[5], InfoSum[6], InfoSum[7], InfoSum[8], InfoSum[9], InfoSum[10], InfoSum[11], InfoSum[12], InfoSum[13], InfoSum[14], InfoSum[15], InfoSum[16]))
	InfoSumFile.close()

### Get .in data for rocks that hit something and write to file
	if writegood==True:
		gooddest=['Jupiter','Io','Europa','Ganymede','Callisto', 'Saturn','Enceladu','Rhea','Titan','Iapetus']
		ind=numpy.array([any(dest1[i]==gooddest[j] for j in range(len(gooddest)) ) for i in range(len(dest1))])
		if (len(name1) > 0):
			name=numpy.array(name1)[ind]
			dest=numpy.array(dest1)[ind]
			#print(name[0],dest[0],time1[0])
			goodin=open(whichdir+'/good.in','a')
			smallin=open(whichdir+'/In/small.in','r')
			SmallLen=MercModule.file_len(whichdir+'/In/small.in')
			smalllines=['' for i in range(SmallLen)]
			for j in range(5):
				smalllines[j]= smallin.readline()
			for j in range(5,SmallLen):
				smalllines[j]= smallin.readline()
				if any(name[i]==smalllines[k].split()[0] and float(time1[i])>30 for i in range(len(name)) for k in range(j-3,j+1)):
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
def makebigchoose(whichdir,whichtime):
	# Needed modules
	import numpy, os
#	from random import random
#	from math import pi, sin, cos
	
	here=os.getcwd()
	
	print('makebig '+whichdir+'/In/big.in  '+whichtime)
	
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

### Write big.in
	bigheader=["Mercury   d=5.427  m= 1.660E-07\n",
	"Venus     d=5.204  m= 2.4476E-06\n",
	"Earth     d=5.515  m= 3.0032E-06\n",
	"Mars      d=3.9335  m= 3.2268E-07\n",
	"Jupiter   d=1.326  m= 9.54266E-04\n",
	"Io        d=3.530  m= 4.491E-08\n",
	"Europa    d=2.99   m= 2.412E-08\n",
	"Ganymede  d=1.94   m= 7.451E-08\n",
	"Callisto  d=1.851  m= 5.409E-08\n",
	"Saturn    d=0.687  m= 2.85717E-04\n",
	"Enceladus d=1.606  m= 5.4321E-11\n",
	"Rhea      d=1.233  m= 1.161E-09\n",
	"Titan     d=1.880  m= 6.76452E-08\n",
	"Iapetus   d=1.088  m= 9.0790E-10\n",
	"Uranus    d=1.318  m= 4.36430E-05\n",
	"Neptune   d=1.638  m= 5.1486E-05\n"]
	bigs=["   0.0  0.0 0.0\n" for i in range(len(bigheader))]	
	bigfile=open(here+'/'+whichdir+'/In/big.in','w')
	bigfile.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
	bigfile.write(") Lines beginning with `)' are ignored.\n")
	bigfile.write(")---------------------------------------------------------------------\n")
	bigfile.write("style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
	bigfile.write(" epoch (in days) = 0.0\n")
	bigfile.write(")---------------------------------------------------------------------\n")
	for i in range(len(big)):
		bigfile.write(bigheader[i])
		bigfile.write("  "+bigxv[i][0]+"  "+bigxv[i][1]+"  "+bigxv[i][2]+"\n")
		bigfile.write("  "+bigxv[i][3]+"  "+bigxv[i][4]+"  "+bigxv[i][5]+"\n")
		bigfile.write(bigs[i])
	bigfile.close()

	timestepfile=open(here+'/'+whichdir+'/timestep.txt','w')
	timestepfile.write(repr(timestep))
	timestepfile.close()


############################################################################
### Pick a random timestep for the big objects and make a new big.in
def makebigrand(whichdir,whichtime):
	# Needed modules
	import MercModule
	import numpy, os
	from random import random
#	from math import pi, sin, cos
	
	here=os.getcwd()
	
	print('makebig '+whichdir+'/In/big.in  '+whichtime)
	
	#constants/variables
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
	
#	big=['Mercury','Venus','Earth','Mars','Jupiter',
#	'Io','Europa','Ganymede','Callisto','Saturn','Uranus','Neptune']
	big=['Mercury','Venus','Earth','Mars', 		'Jupiter','Io','Europa','Ganymede','Callisto',
	'Saturn','Enceladu','Rhea','Titan','Iapetus','Uranus','Neptune']
	bigxv=['' for i in range(len(big))]

### Pick a random timestep and get all big vectors at that point
	AEILen=MercModule.file_len(here+'/'+whichdir+'/In/InitElemFiles/Jupiter.aei')-5
	timestep=5+int(AEILen*random())

# Find the correct timestep for each big thing
	for i in range(len(big)):
		filename=here+'/'+whichdir+'/In/InitElemFiles/'+big[i]+'.aei'
		File=open(filename,'r')
		for j in range(timestep):
			thisline=File.readline().split()
		bigxv[i]=thisline[6:]

### Write big.in
	bigheader=["Mercury   d=5.427  m= 1.660E-07\n",
	"Venus     d=5.204  m= 2.4476E-06\n",
	"Earth     d=5.515  m= 3.0032E-06\n",
	"Mars      d=3.9335  m= 3.2268E-07\n",
	"Jupiter   d=1.326  m= 9.54266E-04\n",
	"Io        d=3.530  m= 4.491E-08\n",
	"Europa    d=2.99   m= 2.412E-08\n",
	"Ganymede  d=1.94   m= 7.451E-08\n",
	"Callisto  d=1.851  m= 5.409E-08\n",
	"Saturn    d=0.687  m= 2.85717E-04\n",
	"Enceladus d=1.606  m= 5.4321E-11\n",
	"Rhea      d=1.233  m= 1.161E-09\n",
	"Titan     d=1.880  m= 6.76452E-08\n",
	"Iapetus   d=1.088  m= 9.0790E-10\n",
	"Uranus    d=1.318  m= 4.36430E-05\n",
	"Neptune   d=1.638  m= 5.1486E-05\n"]
	bigs=["   1.792010121756643E-03  3.102769342987413E-01 -1.834367089631395E-03\n",
	"   4.157775963413559E-03  7.198965919703363E-01 -7.277392122177882E-05\n",
	"   5.681242937724309E-03  9.836767216445764E-01 -7.177966387212160E-06\n",
	"   8.941883860042635E-03  1.548239196455944E+00 -1.272582251831381E-03\n",
	"   2.975007110148294E-02  5.151065133208689E+00  3.596210195046981E-04\n",
	"   2.975344763881742E-02  5.151649762525215E+00 -9.454545719755713E-03\n",
	"   2.975661205240542E-02  5.152197663753571E+00  8.077321734487707E-03\n",
	"   2.976997769835456E-02  5.154511853612028E+00  5.847323539274621E-03\n",
	"   2.980999358082247E-02  5.161440388883488E+00 -2.305346759016693E-03\n",
 	"   5.790966270834483E-02  1.002674727852874E+01 -7.096561661114646E-05\n",
 	"   5.790214154262439E-02  1.002544502905271E+01 -2.924573993367133E-03\n",
 	"   5.791751526309882E-02  1.002810690623725E+01 -4.097004346993902E-03\n",
 	"   5.786595119628821E-02  1.001917886482960E+01  2.515571900187165E-04\n",
 	"   5.792032785500802E-02  1.002859389143027E+01  1.721040392753832E-03\n",
	"   1.119280984160482E-01  1.937974948733879E+01  1.791737569649848E-04\n",
	"   1.744789613275039E-01  3.021009567025112E+01 -2.281823905170589E-05\n"]	
	bigfile=open(here+'/'+whichdir+'/In/big.in','w')
	bigfile.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
	bigfile.write(") Lines beginning with `)' are ignored.\n")
	bigfile.write(")---------------------------------------------------------------------\n")
	bigfile.write("style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
	bigfile.write(" epoch (in days) = 0.0\n")
	bigfile.write(")---------------------------------------------------------------------\n")
	for i in range(len(big)):
		bigfile.write(bigheader[i])
		bigfile.write("  "+bigxv[i][0]+"  "+bigxv[i][1]+"  "+bigxv[i][2]+"\n")
		bigfile.write("  "+bigxv[i][3]+"  "+bigxv[i][4]+"  "+bigxv[i][5]+"\n")
		bigfile.write(bigs[i])
	bigfile.close()

	timestepfile=open(here+'/'+whichdir+'/timestep.txt','w')
	timestepfile.write(repr(timestep))
	timestepfile.close()




############################################################################
### Make a random cluster of small objects around preselected ones
### to write to small.in
def makesmall(whichdir,whichtime,n,whichpl,da,dv):
	n=int(n)
### Needed modules
	import MercModule
	import numpy, os
	from random import random
	from math   import pi, sin, cos

	here=os.getcwd()
	
	print('makesmall '+whichdir+'/In/small.in  '+whichtime+'  '+str(n))

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
		AEILen=MercModule.file_len(filename)-5
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

### Write small.in
	smallheader="    m=1.0012066e-17  r=0.001 d=2.0\n"
	smalls="  0.0 0.0 0.0\n"
	smallfile=open(whichdir+'/In/small.in','w')
	smallfile.write(")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
	smallfile.write(") Lines beginning with `)' are ignored.\n")
	smallfile.write(")---------------------------------------------------------------------\n")
	smallfile.write("style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
	smallfile.write(")---------------------------------------------------------------------\n")
	for j in range(len(small)):
		smallfile.write(small[j]+smallheader)
		smallfile.write("  "+repr(smallxv[j][0])+"  "+repr(smallxv[j][1])+"  "+repr(smallxv[j][2])+"\n")
		smallfile.write("  "+repr(smallxv[j][3])+"  "+repr(smallxv[j][4])+"  "+repr(smallxv[j][5])+"\n")
		smallfile.write(smalls)

	smallfile.close()


############################################################################
### Copy the small rocks from good.in to small.in
def good2small(whichdir,whichtime,n):
	n=int(n)
	# Needed modules
	import MercModule
	import numpy, os
	from random import random
	from math import pi, sin, cos

	here=os.getcwd()
	
	print('good2small '+whichdir+'/good.in  '+whichtime)

	#constants/variables
	AU = 1.496e13					#cm/AU
	day = 24.*3600.					#s/day
#	maxaspread=.57446e9/AU			#in AU
#	maxvspread=5.7446*day/AU		#in AU/day
	
	goodin=open(whichdir+'/good.in','r')
	GoodLen=MercModule.file_len(whichdir+'/good.in')
	if (GoodLen%4 != 0):
		print('??? good.in length = '+str(GoodLen))
	ngood=GoodLen/4
	header, pos, vel, s = [],[],[],[]
	for j in range(ngood):
		header.append(goodin.readline() )
		pos.append(goodin.readline().split() )
		vel.append(goodin.readline().split() )
		s.append(goodin.readline() )
	posf=[ [float(pos[i][j]) for j in range(3)] for i in range(len(pos)) ]
	velf=[ [float(vel[i][j]) for j in range(3)] for i in range(len(vel)) ]

	name=['M'+str(i) for i in range(ngood)]
	small=['' for i in range(ngood)]

	for j in range(len(small)):
		if (j < ngood):
			small[j]=[name[j], pos[j], vel[j]]
### For two good files:
#	print('good2small '+whichdir+'/good2.in  '+whichtime+'  '+str(n2))
#	good2in=open(whichdir+'/good2.in','r')
#	Good2Len=MercModule.file_len(whichdir+'/good2.in')
#	if (Good2Len%4 != 0):
#		print('??? good2.in length = '+str(Good2Len))
#	ngood2=Good2Len/4
#	header2, pos2, vel2, s2 = [],[],[],[]
#	for j in range(ngood2):
#		header2.append(good2in.readline() )
#		pos2.append(good2in.readline().split() )
#		vel2.append(good2in.readline().split() )
#		s2.append(good2in.readline() )
# 	posf2=[ [float(pos2[i][j]) for j in range(3)] for i in range(len(pos2)) ]
#	velf2=[ [float(vel2[i][j]) for j in range(3)] for i in range(len(vel2)) ]
#
#	name2=['M'+str(i+n1) for i in range(ngood2)]
#	small2=['' for i in range(ngood2)]
#
#	for j in range(len(small2)):
#		if (j < ngood2):
#			small2[j]=[name2[j], pos2[j], vel2[j]]
### Can eventually add some randomness here?
#		else:
#			phi1=2*pi*random()
#			theta1=-pi/2+pi*random()
#			r=maxaspread*random()
#			v=maxvspread*random()
#			phi2=2*pi*random()
#			theta2=-pi/2+pi*random()
#			x=r*cos(phi1)*cos(theta1)
#			y=r*sin(phi1)*cos(theta1)
#			z=r*sin(theta1)
#			u=v*cos(phi2)*cos(theta2)
#			v=v*sin(phi2)*cos(theta2)
#			w=v*sin(theta2)
#			if ( j-ngood< (n-ngood)//ngood):
#				small[j]=[name[j],pos[0][0]+x,pos[0][1]+y,pos[0][2]+z,
#				vel[0][0]+u,vel[0][1]+v,vel[0][2]+w]
# ?????????????????????????????


### Write small.in
	smallheader="    m=1.0012066e-17  r=0.001 d=2.0\n"
	smalls="  0.0 0.0 0.0\n"
	smallfile=open(whichdir+'/small.in','w')
	smallfile.write(")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
	smallfile.write(") Lines beginning with `)' are ignored.\n")
	smallfile.write(")---------------------------------------------------------------------\n")
	smallfile.write("style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
	smallfile.write(")---------------------------------------------------------------------\n")
	for j in range(len(small)):
		smallfile.write(name[j]+smallheader)
		smallfile.write("  "+str(small[j][1][0])+"  "+str(small[j][1][1])+"  "+str(small[j][1][2])+"\n")
		smallfile.write("  "+str(small[j][2][0])+"  "+str(small[j][2][1])+"  "+str(small[j][2][2])+"\n")
		smallfile.write(smalls)
#	for j in range(len(small2)):
#		smallfile.write(name2[j]+smallheader)
#		smallfile.write("  "+str(small2[j][1][0])+"  "+str(small2[j][1][1])+"  "+str(small2[j][1][2])+"\n")
#		smallfile.write("  "+str(small2[j][2][0])+"  "+str(small2[j][2][1])+"  "+str(small2[j][2][2])+"\n")
#		smallfile.write(smalls)

	smallfile.close()










