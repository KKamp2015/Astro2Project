#Packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import copy
import sys
#Variables
H0= 1/13.8e9#Hubble constant needs to be in s^-1
iMISM=10**12 #Iniital mass of ISM
iMStar=10**12
StarburstTime=5
Starburst=False                            #Here
StarburstFrac=0.5
#Star Class
class Star():
	'''
	inputs
	------

	'''
	def __init__(self, Mass,weight): # initaionliing the class with mass
		self.mass=Mass #giving star class atribute mass
		self.StartMass=Mass
		self.calc_rest() # profoming calc_rest fucntion definded below
		self.TLife=1.1*self.TMS
		self.Dead=False # Giving star atribute dead and initalizing it to false i.e. star is alive
		self.Blown=False
		self.Weight=weight
		self.SNIa=False
	
	def calc_rest(self): #defining calc_rest 
		self.TMS=10**10 *(self.mass)**(-3.5) # calculating the time on the main sequnce 
	
	def Wind(self): # setting up wind procedure for stars past MS
		if self.Blown==True: #checking for wind 
			print("Star's wind already accounted for") # rpinting error message
			return 0,0,0,0,0,0 # retunring 0 mass to ism
		else:
			if self.mass>=3: # for stars greater than or equal to 3 solar masses
				self.Blown==True
				winda=0.015*self.mass
				windb=0.01*self.mass
				self.mass=self.mass-sum([0,0,winda,winda,windb,0.00])
				return 0,0,winda*self.Weight,winda*self.Weight,windb*self.Weight,0.00
			else:
				self.Blown==True
				return 0,0,0,0,0,0

	def Kill(self): # setting up Kill funciton
		iMass=copy.deepcopy(self.mass)
		rem=0
		H=0
		He=0
		He2=0
		if self.Dead==True: # checking if star has not already been killed
			print("Star already Dead ", self.mass) # Printing warning
			return 0,0,0,0,0,0 # retuening 0 mass to ism
		else:
			if self.StartMass>=8: # for massive stars
				self.Dead=True #setting the dead atribute to true
				if self.StartMass<50: # for less massive set mass for neutron star as reminent
					self.mass=1.4 # setting mass for reminent
				else:
					self.mass=3.0 # setting mass of remeninat for very massive stars.
				rem=iMass-(0.1+0.1+0.1+0.03)-self.mass
				return 0.76*rem*self.Weight,0.24*rem*self.Weight,0.1*self.Weight,0.1*self.Weight,0.1*self.Weight,0.03*self.Weight # returning masses of elemens given to ISM  [H, He, C,N,O,Fe]
			elif self.StartMass<=3:
				self.Dead=True# setting dead to True
				self.mass=0.6
				H=(self.StartMass-self.mass) *0.76
				He=(self.StartMass-self.mass) *0.24
				excess=0.000
				excess=(H+He+self.mass-self.StartMass)
				return H*self.Weight,He*self.Weight+excess*self.Weight,0,0,0,0
			elif self.StartMass<=8: # for non massive stars
				self.Dead=True# setting dead to True
				self.mass=0.6
				H=(s.StartMass-self.mass) *0.76
				He=(s.StartMass-self.mass) -H
				if (H+He+self.mass)-self.StartMass !=0:
					dale=0
					#print('\nMed')
					#print(self.StartMass)
					#print((H+He)-self.StartMass)
					#sys.exit()
				return H*self.Weight,He*self.Weight,0,0,0,0# returning masses of elemens given to ISM  [H,He, C,N,O,Fe]

#Functions
LookBack= lambda z: (2/(3*H0))*(1-(1/(1+z)**(3/2))) # creating fucntion to convert redshift z to look back time

def StarBurst(Gal,f,ISM):
	print('Pre H mass is: ',ISM)
	test=ISM
	Masses1,Weights1=IMF(f*ISM)
	ISM=ISM-(f*ISM)
	BurstMass=0
	for m in range(len(Masses1)):
		tStar=Star(Masses1[m],Weights1[m])
		BurstMass+=(Masses1[m]*Weights1[m])
		tStar.TMS+=t
		tStar.TLife+=t
		Gal.append(tStar)
	print('Burst Mass is:',BurstMass)
	print('Post H mass is : ',ISM)
	return ISM, Gal

#time array
Z=np.linspace(0,12,1000) #creating array of reshifts staring at z=12 and going to 0 with 100000 elements
Time=LookBack(Z) #convering redshifts to lookback times using funciton amde above
#Time=np.abs(Time1-max(Time1))[::-1]
#print(Time)
#IMF
def IMF(Mass):
	Mass=np.abs(Mass)
	Masses=[]
	Weights=[]
	Mran=np.linspace(0.1,100,10000)
	#Const=36354.560545987355
	C= 3.35*Mass/(100**1.35)/10
	N= lambda M1: C*M1**(-(1.35))
	for i in range(1,len(Mran)):
	    #Masses.append((Mran[i-1]+Mran[i])/2)
	    Masses.append(Mran[i])
	    Weights.append(N(Mran[i]))
	#print(Masses)
	return Masses,Weights # returning IMF mass array
print('Creating IMF')
MassArray,WeightArray=IMF(iMStar) #Getting mass array from above function 
plt.bar(MassArray,WeightArray)
plt.yscale('log')
plt.savefig('InitalMassFunction.png')
plt.close()
totalmass=sum([a*b for a,b in zip(MassArray,WeightArray)])
'''
M=Table()
M["Masses"]=MassArray
M.write('Masses1000.csv')
'''
print('Populating Galaxy')
Galaxy=[] #creating an empty galaxy
for i in range(len(MassArray)): #populating galazy with stars from IMF
	Galaxy.append(Star(MassArray[i],WeightArray[i])) #appending Galaxyarray with star objects
#ISM loop
MISM=[iMISM] #mass arrawy over time of the ISM
cMISM=iMISM # current ISM mass set to inital mass of ISM for t=0
MSTARS=[totalmass] #cumluative mass array of Stars
H=[0.76*iMISM]
He=[0.24*iMISM]
C=[0]
N=[0]
O=[0]
Fe=[0]
Tarr=[0]
Zarr=[0]
ctempH=0
ctempHe=0
ctempc=0
ctempn=0
ctempo=0
ctempfe=0
print('Running Simulation')
for t in Time[1:]: # running the time
	tMSTARS=0 #setting the timestep mass of stars to 0
	if t>LookBack(12-StarburstTime) and Starburst==False:
		print('\nSTARBURST!!!')
		print('Start ISM Mass: ',sum([ctempH,ctempHe,ctempc,ctempn,ctempo,ctempfe]))
		ctempH,Galaxy=StarBurst(Galaxy,StarburstFrac,ctempH)
		print('End ISM Mass: ',sum([ctempH,ctempHe,ctempc,ctempn,ctempo,ctempfe]))
		Tarr.append(t)
		Zarr.append(Z[Time==t][0])
		H.append(ctempH)
		He.append(ctempHe)
		C.append(ctempc)
		N.append(ctempn)
		O.append(ctempo)
		Fe.append(ctempfe)
		MISM.append(cMISM)# appedning ISM contribitions time steps
		MSTARS.append(tMSTARS)
		Starburst=True
	for star in range(len(Galaxy)): # running loop over all stars in galaxy
		s=Galaxy[star]
		iMass=copy.deepcopy(s.mass)
		temph=0
		temphe=0
		tempc=0 
		tempn=0 
		tempo=0 
		tempfe=0
		temph1=0
		temphe1=0
		tempc1=0 
		tempn1=0 
		tempo1=0 
		tempfe1=0
		temph2=0
		temphe2=0
		tempc2=0 
		tempn2=0 
		tempo2=0 
		tempfe2=0 # initaling temp for ISM to 0 if star doesnt die or wind
		tempISM=0
		if s.TLife<=t and s.Dead==False: #case for stars that are dead
			if s.Dead==False:
				temph1,temphe1,tempc1,tempn1,tempo1,tempfe1=s.Kill() # killing stars and getting contrubition to ISM
			if s.Blown==False:
				temph2,temphe2,tempc2,tempn2,tempo2,tempfe2=s.Wind()
			temph=temph1+temph2
			temphe=temphe1+temphe2
			tempc=tempc1+tempc2
			tempn=tempn1+tempn2
			tempo=tempo1+tempo2
			tempfe=tempfe1+tempfe2
			tempISM=sum([temph,temphe,tempc,tempn,tempo,tempfe])
		elif s.TMS<=t and s.Blown==False:
				(temph,temphe,tempc,tempn,tempo,tempfe)==s.Wind() # accounting for stars winds' contribtion to ISM
				tempISM=sum([temph,temphe,tempc,tempn,tempo,tempfe])
		if s.StartMass<=8 and s.StartMass>=3 and s.TLife+2e8<t and s.SNIa==False: #SNe Ia
			s.SNIa=True
			s.mass=0
			tempc+=0.1*s.Weight*10**-5
			tempn+=0.1*s.Weight*10**-5
			tempo+=0.1*s.Weight*10**-5
			tempfe+=0.3*s.Weight*10**-5
			tempISM=sum([temph,temphe,tempc,tempn,tempo,tempfe])
		cMISM+=tempISM #adding contrubuted masses to ism:
		ctempH+=temph
		ctempHe+=temphe
		ctempc+=tempc
		ctempn+=tempn
		ctempo+=tempo
		ctempfe+=tempfe
		tMSTARS+=(s.mass*s.Weight)
		Galaxy[star]=s
	H.append(ctempH)
	He.append(ctempHe)
	C.append(ctempc)
	N.append(ctempn)
	O.append(ctempo)
	Fe.append(ctempfe)
	Tarr.append(t)
	Zarr.append(Z[Time==t][0])
	MISM.append(cMISM)# appedning ISM contribitions time steps
	MSTARS.append(tMSTARS) #appedning to mass of stars for time steps
	spot=np.where(Time==t)
	print(' '+str(spot[0][0])+'/'+str(len(Time)),end="\r")
print('\nCreating Output Table')
T=Table()
T['z']=Zarr
T['Time']=Tarr
T['ISM']=MISM
T['Stars']=MSTARS
T['H']=H
T['He']=He
T['C']=C
T['N']=N
T['O']=O
T['Fe']=Fe
T.write("OutputBigSB.csv",overwrite=True)                               #Here
print('Creating Plots',end="\n")
print('Plotting ISM mass vs. z',end="\r")
plt.clf()
plt.plot(T["z"][::-1],T['ISM'],label='Total')
plt.plot(T["z"][::-1],T['H'],label='Hydrogen')
plt.plot(T["z"][::-1],T['He'],label='Helium')
plt.plot(T["z"][::-1],T['C'],label='Carbon')
plt.plot(T["z"][::-1],T['N'],label='Nitrogen')
plt.plot(T["z"][::-1],T['O'],label='Oxygen')
plt.plot(T["z"][::-1],T['Fe'],label='Iron')
plt.xlabel('Redshift')
plt.xlim([12.1,0])
plt.yscale('log')
plt.legend()
plt.ylabel(r'ISM Mass [M$_{\odot}]$')
plt.title('ISM Mass vs. Redshift')
plt.savefig('ISMvzSB.png',dpi=300)                               #Here
print('Plotting Stellar mass vs. z')
plt.clf()
plt.plot(T["z"][::-1],T['Stars'],label='Stars')
plt.xlabel('Redshift')
plt.xlim([12.1,0])
plt.yscale('log')
plt.ylabel(r'Stellar Mass [M$_{\odot}]$')
plt.title('Stellar Mass vs. Redshift')
plt.savefig('StarsvzSB.png',dpi=300)                               #Here
plt.close()

'''
endmass=[]
for s in Galaxy:
	endmass.append(s.mass)

M=Table()
M['StartMass']=MassArray
M['EndMass']=endmass
M.write('Masses.csv',overwrite=True)
plt.clf()
Diff=[np.abs(MassArray[i]-endmass[i]) for i in range(len(MassArray)) ]
plt.bar(range(len(Diff)),Diff)
plt.savefig('MassDiff.png')
'''