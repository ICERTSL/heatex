import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import solveblockLU

def lengthsolver(T,P,T_hot_in,T_cold_in,P_hot_in,P_cold_in,m,N_ch,Length,Me):
	print("hello length")	
	
	d = 2e-3        # Channel square side in metres
	res = 10
	res1old=0
	#Length = 1
	Recold_in = 50000;
	s = 1.5e-3;             # Spacing between channels in metres
	width=s;

	nc = 150/m;
	r = Me-1;               #Number of elements in a row
	N_layers= 4*N_ch +1     # Total number of layers 
	enthold=np.zeros([N_layers,r+1])
	enth=np.zeros([N_layers,r+1])
	deltaz=Length/r
	delx=deltaz
	K_solid= 21 # Thermal conductivity of Aluminium in W/(m K)
	R_axial=deltaz/(K_solid*s*d)
	R_longi= 0.5*s/(K_solid*deltaz*d)

	R_hot=np.zeros([N_ch,r]);
	R_cold=np.zeros([N_ch,r]);
	ehin = 0
	ecin = 0
	for i in range(1,N_layers-1,4):
		ehin = ehin+PropsSI("H","T",T_hot_in,"P",P_hot_in,"CO2")
	for i in range(3,N_layers-1,4):
		ecin = ecin+PropsSI("H","T",T_cold_in,"P",P_cold_in,"CO2")
		
	A = np.zeros([N_layers,N_layers,r,r],dtype = "float16")

	for i in range (0,N_layers,1):
		for j in range (0,r-1,1):
			A[i,i,j+1,j] = -(1/R_axial)
			A[i,i,j,j+1] = -(1/R_axial)
	
	b = np.zeros([N_layers,r])


	#Turbulent fully developed smooth walled flow
	def friction_factor(x):
		f=np.power((0.79*(np.log(x))-1.64),-2) # use blassius eqn?
		return(f)

#Finds the Nusselt number given the Reynolds number
	def Nu(x,y,z):
		f=friction_factor(z);
		mu=PropsSI('V','T',x,'P',y,'co2');
		cp=PropsSI('C','T',x,'P',y,'co2');
		k=PropsSI('L','T',x,'P',y,'co2');
		# Pr=refpropm('^','T',T,'P',P,'co2');
		Pr = (mu*cp)/k;
		Pr2 = np.power(Pr, 0.66666)
		f2 = f/8
		f2 = np.power(f2,0.5)
		return(((f/8)*(z-1000)*Pr)/(1+12.7*f2*(Pr2-1)))

#finds Cp
	def Cp(x,y):
		C=PropsSI('C','T',x,'P',y,'co2')
		return(C)

	def K_liquid(x,y):
		k=PropsSI('L','T',x,'P',y,'co2')
		return(k)
	
	Told = np.zeros([N_layers,r])	
	count=0
	while(res>0.0001): #abs(np.max(T-Told))>0.1
		print("Iteration no.");print(count);print("residual");print(np.max(abs(T-Told)));
		for i in range(1,N_layers,4):
			for j in range (0,r):
			       	if(j==0):
		        		T_mean=0.5*(T[i,j]+T_hot_in)
				

	        
		       		if(j>0):
		       			T_mean=0.5*(T[i,j]+T[i,j-1])
        	#print(i);print(j)
		
		       		density_hot=PropsSI('D','T',T_mean,'P',0.5*(P[i,j]+P[i,j+1]),'CO2')
		       		Re_hot=m/(d*PropsSI('V','T',T_mean,'P',0.5*(P[i,j]+P[i,j+1]),'co2'))
		       		deltaphot=friction_factor(Re_hot)*(0.5*m*m*deltaz/(density_hot*(np.power(d,5))))   #this value is in Pa
		       		P[i,j+1]= P[i,j] - deltaphot   # Working with KPa
		       		P_mean=0.5*(P[i,j]+P[i,j+1]); I1 = int((i-1)/4); 	
			
		       		R_hot[I1,j]=R_longi+(1/(K_liquid(T_mean,P_mean)*Nu(T_mean,P_mean,Re_hot)*(deltaz))) # define functions Nu and ffactor
		       		#Re_duph[i,j]=Re_hot #See about this later
		       		A[i,i-1,j,j]= -1/R_hot[I1,j] #L(i-1).maindiagonal(j)
		       		A[i,i,j,j]=(1/R_hot[I1,j])+m*Cp(T_mean,P_mean) #DEFINE Cp FUNCTION# M(i).maindiagonal(j)
		       		if(j>0):
		       			A[i,i,j,j-1]=(1/R_hot[I1,j])-m*Cp(T_mean,P_mean) #M(i).lowerdiagonal(j-1)
	        	#print("hello")
		       		A[i+1,i,j,j]=-0.5/(R_hot[I1,j]) #L(i).maindiagonal(j)
   
	
			b[i,0]=T_hot_in*((-1/R_hot[I1,1])+m*Cp(0.5*(T_hot_in+T[i,1]),0.5*(P_hot_in+P[i,1])))
			A[i,i+1,:,:]=A[i,i-1,:,:]; #print(i);print(A[i,i+1,:,:])#U(i)=L(i-1);
		      #Replace i with i+2
			for j in range(r-1,-1,-1):
				if(j==r-1):
			            T_mean=0.5*(T[i+2,j]+T_cold_in);
        	
				if(j<r-1):
	        			T_mean=0.5*(T[i+2,j]+T[i+2,j+1]); #print(i);print(j)
        	
				density_cold=PropsSI('D','T',T_mean,'P',0.5*(P[i+2,j]+P[i+2,j+1]),'co2');I1 = int((i-1)/4)
				Re_cold=m/(d*PropsSI('V','T',T_mean,'P',0.5*(P[i+2,j]+P[i+2,j+1]),'co2'))
				deltapcold=friction_factor(Re_cold)*(0.5*m*m*deltaz/(density_cold*(d*d*d*d*d)))  #this value is in Pa
				P[i+2,j]=P[i+2,j+1]-deltapcold # take care of Pa or Kpa
				P_mean=0.5*(P[i+2,j]+P[i+2,j+1])
				R_cold[I1,j]=R_longi+(1/(K_liquid(T_mean,P_mean)*Nu(T_mean,P_mean,Re_cold)*(deltaz)))
				#Re_dupc[i,j]=Re_cold
				A[i+2,i+1,j,j]=-1/R_cold[I1,j] #L(i+1).maindiagonal(j)
				A[i+3,i+2,j,j]=-0.5/R_cold[I1,j] #L(i+2).maindiagonal(j)
				A[i+2,i+2,j,j]=(1/R_cold[I1,j])+m*Cp(T_mean,P_mean) #M(i+2).maindiagonal(j)
				if(j<r-1):
					A[i+2,i+2,j,j+1]=(1/R_cold[I1,j])-m*Cp(T_mean,P_mean) #M(i+2).upperdiagonal(j)
				#end if
				A[i+1,i+2,j,j]=-0.5/(R_cold[I1,j]) #U(i+1).maindiagonal(j)
				if(j>=1 and j<=r-1):
	        			A[i+1,i+1,j,j]=(1/R_cold[I1,j])+(1/R_hot[I1,j])+2/R_axial #M(i+1).maindiagonal(j)
	        		#end if
				if(j==0 or j==r-1):
					A[i+1,i+1,j,j]=(1/R_cold[I1,j])+(1/R_hot[I1,j])+1/R_axial #M(i+1).maindiagonal(j)
	        		#end if
	        		#if(i<N_layers-3)
	        		#L(i+3).maindiagonal(j)=-0.5*R_cold((i+2)/4,j);
	        		#end if
			#print("hello1")
				if(i>1): 
					A[i-1,i,j,j]=-0.5/R_hot[int((i-1)/4),j] #U(i-1).maindiagonal(j)
					if(j>=1 and j<=r-2):
						A[i-1,i-1,j,j]=(1/R_cold[int((i-5)/4),j])+(1/R_hot[int((i-1)/4),j])+2/R_axial #M(i-1).maindiagonal(j)
        	#end if
					if(j==0 or j==r-1):
						A[i-1,i-1,j,j]=(1/R_cold[int((i-5)/4),j])+(1/R_hot[int((i-1)/4),j])+1/R_axial #M(i-1).maindiagonal(j)
				#print("hello2")
        			#endif

			if(i>1):
				for j in range (0,r-1):		
					A[i-1,i,j+1,j]=A[i-1,i,j+1,j+1] #U(i-1).lowerdiagonal=U(i-1).maindiagonal(2:r);

				b[i-1,0]=0.5*T_hot_in/(R_hot[int((i-1)/4),1])
				b[i-1,r-1]=0.5*T_cold_in/(R_cold[int((i-5)/4),1])
	#end if
	

	
			A[i+2,i+3,:,:] = A[i+2,i+1,:,:] #U(i+2)=L(i+1)

  
			b[i+2,r-1]=-T_cold_in*(1/R_cold[int((i+2)/4),r-1]-m*Cp(0.5*(T_cold_in+T[i+2,r-1]),0.5*(P_cold_in+P[i+2,r-1])));	
	
			for j in range (0,r-1):
				A[i+1,i,j+1,j]= A[i+1,i,j+1,j+1] #L(i).lowerdiagonal=L(i).maindiagonal(2:r);	
				A[i+3,i+2,j,j+1]= A[i+3,i+2,j,j] #L(i+2).upperdiagonal=L(i+2).maindiagonal(1:r-1);
				A[i+1,i+2,j,j+1]= A[i+1,i+2,j,j] #U(i+1).upperdiagonal=U(i+1).maindiagonal(1:r-1);
    
			b[i+1,0]=0.5*T_hot_in/(R_hot[int((i-1)/4),0]);
			b[i+1,r-1]=0.5*T_cold_in/(R_cold[int((i-1)/4),r-1]);
			b[i+2,r-1]=T_cold_in*(m*Cp(0.5*(T_cold_in+T[i+2,r-1]),0.5*(P_cold_in+P[i+2,r-1]))-1/(R_cold[int((i-1)/4),r-1]))
	#end for loop

	#first layer input and last layer input
		for j in range (0,r):
			A[0,1,j,j] = -0.5/(R_hot[0,j]) #U(1).maindiagonal(j)=-0.5/(R_hot(1,j));
			A[N_layers-1,N_layers-2,j,j]=-0.5/(R_cold[N_ch-1,j]) #L(N_layers-1).maindiagonal(j)=-0.5/(R_cold(N_ch,j));
			if(j>=1 and j<=r-2):
				A[0,0,j,j] = (1/R_hot[0,j])+2/R_axial #M(1).maindiagonal(j)=(1/R_hot(1,j))+2/R_axial;
				A[N_layers-1,N_layers-1,j,j]= (1/R_cold[N_ch-1,j])+2/R_axial #M(N_layers).maindiagonal(j)=(1/R_cold(N_ch,j))+2/R_axial;
		#end if
			if(j==0 or j==r-1):
				A[0,0,j,j] = (1/R_hot[0,j])+1/R_axial #M(1).maindiagonal(j)=1/R_hot(1,j)+1/R_axial;
				A[N_layers-1,N_layers-1,j,j]= (1/R_cold[N_ch-1,j])+1/R_axial #M(N_layers).maindiagonal(j)=(1/R_cold(N_ch,j))+1/R_axial;
		#end if
		for j in range (0,r-1):
			A[N_layers-1,N_layers-2,j,j+1]= A[N_layers-1,N_layers-2,j,j]#L(N_layers-1).upperdiagonal=L(N_layers-1).maindiagonal(1:r-1);
			A[0,1,j+1,j]=A[0,1,j+1,j+1] #U(1).lowerdiagonal=U(1).maindiagonal(2:r);

		b[0,0]=0.5*T_hot_in/R_hot[0,0]
		b[N_layers-1,r-1]=T_cold_in*0.5/R_cold[N_ch-1,r-1]
		Told=T.copy()
	#T=solveblockLU(M,L,U,b,r,N_layers);
	#print("Told")
	#print(Told)
	#print("limit 1")
		for i in range(1,N_layers-1,4): #CHOR FIX, correct this later
			for j in range(0,r-1):	
				A[i,i,j,j+1] = 0
				A[i+2,i+2,j+1,j]=0
	#print("limit 2")

		solveblockLU.solve(A,r,N_layers,b,T,P,m,ehin,ecin,res)
	#print("limit 3")
		ecout=0
		ehout=0
		for i in range(1,N_layers-1,4):
			ehout= ehout+PropsSI('H','T',T[i,r-1],'P',0.5*(P[i,r-1]+P[i,r]),'co2')
		for i in range(3,N_layers-1,4):
			ecout= ecout+PropsSI('H','T',T[i,0],'P',0.5*(P[i,r-1]+P[i,r]),'co2')
		
		res1= abs(abs((ehin-ehout)/(ecout-ecin))-1) 
	
		res = abs(res1-res1old)
		res1old=res1
	#print("final ehin ehout ecin ecout")
	#print(ehin);print(ehout);print(ecin);print(ecout)
	#print("deltah deltac")
	#print(ehin-ehout);print(ecout-ecin)
	
		count = count+1
	print("T final")
	print(T)
	print("final residual")
	print(res)
	print("final ehin ehout ecin ecout")
	print(ehin);print(ehout);print(ecin);print(ecout)
	print("Total heat transfer (KW)")
	print(0.001*m*(ehin - ehout))
	return(0.001*m*(ehin - ehout))

	
