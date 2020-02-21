#Packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
#Variables
H0= 1/13.8e9#Hubble constant needs to be in s^-1
iMISM=10**12 #Iniital mass of ISM
#Star Class
class Star():
	'''
	inputs
	------

	'''
	def __init__(self, Mass,weight): # initaionliing the class with mass
		self.mass=Mass #giving star class atribute mass
		self.calc_rest() # profoming calc_rest fucntion definded below
		self.Dead=False # Giving star atribute dead and initalizing it to false i.e. star is alive
		self.Blown=False
		self.Weight=weight
	
	def calc_rest(self): #defining calc_rest 
		self.TMS=10**10 *(self.mass)**(-3.5) # calculating the time on the main sequnce 
	
	def Wind(self): # setting up wind procedure for stars past MS
		if self.Blown==True: #checking for wind 
			print("Star's wind already accounted for") # rpinting error message
			return 0,0,0,0 # retunring 0 mass to ism
		else:
			if self.mass>=3: # for stars greater than or equal to 3 solar masses
				self.Blown==True
				winda=0.015*self.mass*self.Weight
				windb=0.01*self.mass*self.Weight
				return winda,winda,windb,0.00
			else:
				return 0,0,0,0

	def Kill(self): # setting up Kill funciton
		if self.Dead==True: # checking if star has not already been killed
			print("Star already Dead ", self.mass) # Printing warning
			return 0,0,0,0 # retuening 0 mass to ism
		else:
			if self.mass>8: # for massive stars
				self.Dead=True #setting the dead atribute to true
				if self.mass<50: # for less massive set mass for neutron star as reminent
					self.mass=1.4 # setting mass for reminent
				else:
					self.mass=3.0 # setting mass of remeninat for very massive stars.
				return 0.1*self.Weight,0.1*self.Weight,0.1*self.Weight,0.3*self.Weight # returning masses of elemens given to ISM  [C,N,O,Fe]
			elif self.mass<=8: # for non massive stars
				self.Dead=True# setting dead to True
				self.mass=0
				return 0.3*self.Weight,0.3*self.Weight,0.3*self.Weight,0.3*self.Weight# returning masses of elemens given to ISM  [C,N,O,Fe]

#Functions
LookBack= lambda z: (2/(3*H0))*(1-(1/(1+z)**(3/2))) # creating fucntion to convert redshift z to look back time


#time array
Z=np.linspace(0,12,1000) #creating array of reshifts staring at z=12 and going to 0 with 1000 elements
Time=LookBack(Z) #convering redshifts to lookback times using funciton amde above
#print(Time)
#IMF
def IMF():
	Masses=[]
	Weights=[]
	Mran=np.linspace(0.1,100,10000)
	#Const=36354.560545987355
	C= lambda Min, Max,total: total/(-0.35)*(Max**(-(0.35))-Min**(-(0.35)))
	N= lambda M1, M2: C(0.1,100,10**12)/(-1.35) *(M2**(-(1.35))-M1**(-(1.35)))
	for i in range(1,len(Mran)):
	    Masses.append((Mran[i-1]+Mran[i])/2)
	    Weights.append(N(Mran[i-1],Mran[i]))
	        
	return Masses,Weights # returning IMF mass array
print('Creating IMF')
MassArray,WeightArray=IMF() #Getting mass array from above function 
plt.bar(MassArray,WeightArray)
plt.yscale('log')
plt.savefig('InitalMassFunction.png')
'''
M=Table()
M["Masses"]=MassArray
M.write('Masses1000.csv')
'''
print(' Populating Galaxy')
Galaxy=[] #creating an empty galaxy
for i in range(len(MassArray)): #populating galazy with stars from IMF
	Galaxy.append(Star(MassArray[i],WeightArray[i])) #appending Galaxyarray with star objects
#ISM loop
MISM=[] #mass arrawy over time of the ISM
cMISM=iMISM # current ISM mass set to inital mass of ISM for t=0
MSTARS=[] #cumluative mass array of Stars
C=[]
N=[]
O=[]
Fe=[]
ctempc=0
ctempn=0
ctempo=0
ctempfe=0
print('Running Simulation')
for t in Time: # running the time
	tMSTARS=0 #setting the timestep mass of stars to 0
	for s in Galaxy: # running loop over all stars in galaxy
		tempc=0 
		tempn=0 
		tempo=0 
		tempfe=0 # initaling temp for ISM to 0 if star doesnt die or wind
		if s.TMS*1.1<t: #case for stars that are dead
			if s.Dead==False:
				tempc1,tempn1,tempo1,tempfe1=s.Kill() # killing stars and getting contrubition to ISM
			if s.Blown==False:
				tempc2,tempn2,tempo2,tempfe2=s.Wind()
			tempc=tempc1+tempc2
			tempn=tempn1+tempn2
			tempo=tempo1+tempo2
			tempfe=tempfe1+tempfe2
		elif s.TMS<t and s.Blown==False:
			(tempc,tempn,tempo,tempfe)==s.Wind() # accounting for stars winds' contribtion to ISM
		cMISM+=sum([tempc,tempn,tempo,tempfe]) #adding contrubuted masses to ism:
		ctempc+=tempc
		ctempn+=tempn
		ctempo+=tempo
		ctempfe+=tempfe
		tMSTARS+=(s.mass*s.Weight)
	C.append(ctempc)
	N.append(ctempn)
	O.append(ctempo)
	Fe.append(ctempfe)
	MISM.append(cMISM)# appedning ISM contribitions time steps
	MSTARS.append(tMSTARS) #appedning to mass of stars for time steps
	spot=np.where(Time==t)
	print(' '+str(spot[0][0])+'/'+str(len(Time)),end="\r")
print('\nCreating Output Table')
T=Table()
T['z']=Z
T['Time']=Time
T['ISM']=MISM
T['Stars']=MSTARS
T['C']=C
T['N']=N
T['O']=O
T['Fe']=Fe
T.write("OutputBig.csv",overwrite=True)
print('Creating Plots',end="\n")
print(' Plotting ISM mass vs. z',end="\r")
plt.clf()
plt.plot(T["z"][::-1],T['ISM'],label='ISM')
plt.xlabel('Redshift')
plt.xlim([12.1,0])
plt.yscale('log')
plt.ylabel(r'ISM Mass [M$_{\odot}]$')
plt.title('ISM Mass vs. Redshift')
plt.savefig('ISMvz.png',dpi=300)
print(' Plotting Stellar mass vs. z')
plt.clf()
plt.plot(T["z"][::-1],T['Stars'],label='Stars')
plt.xlabel('Redshift')
plt.xlim([12.1,0])
plt.yscale('log')
plt.ylabel(r'Stellar Mass [M$_{\odot}]$')
plt.title('Stellar Mass vs. Redshift')
plt.savefig('Starsvz.png',dpi=300)
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