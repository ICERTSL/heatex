import numpy as np  #### DIAG4 I.E I=3
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
	
def solve(A,r,N_layers,b,T):
	X=np.zeros([1,N_layers-1	,r,r])
	diagonal=np.zeros([1,N_layers,r,r])
	lower=np.zeros([1,N_layers-1,r,r])
	m=np.zeros([1,N_layers,r,r])	
	#diagonal[0,k-1,:,:] = A[k-1,k-1,:,:]
	for k in range (0,N_layers-1):	
		diagonal[0,k,:,:] = A[k,k,:,:]
		lower[0,k,:,:] = A[k+1,k,:,:]
		m[0,k,:,:] = A[k,k,:,:]
	m[0,N_layers-1,:,:]= A[N_layers-1,N_layers-1,:,:]
	diagonal[0,N_layers-1,:,:]= A[N_layers-1,N_layers-1,:,:]
		
	for k in range (1,N_layers):
				
		X[0,k-1,:,:]= np.matmul(np.linalg.inv(diagonal[0,k-1,:,:]),A[k-1,k,:,:]) #multiply(inv(diagonal(:,:,k-1)),U(k-1),r);	
		diagonal[0,k,:,:]=m[0,k,:,:]-np.matmul(lower[0,k-1,:,:],X[0,k-1,:,:])#diagonal(:,:,k)=m(:,:,k)-lower(:,:,k-1)*P(:,:,k-1);		
				
		#X[0,k-1,:,:]= np.matmul(A[k,k-1,:,:],A[k-1,k,:,:])	
		
	Tpseudo=np.zeros([N_layers,r]);
	#for k in range (0,r):
	Tpseudo[0,:] = np.matmul((np.linalg.inv(diagonal[0,0,:,:])),np.transpose(b[0,:]))

#for k=2:N_layers
#    a=diagonal(:,:,k);
#    Tpseudo(k,1:r)=inv(a)*(transpose(b(k,1:r))-lower(:,:,k-1)*transpose(Tpseudo(k-1,1:r)));
#end

	for k in range(1,N_layers):
		#a=diagonal[0,k,:,:]
		#a = np.linalg.inv(a)
		#print("a no.");print(k);print(np.matmul(a,np.transpose(b[k,:])));print(np.matmul(lower[0,k-1,:,:],(np.transpose(Tpseudo[k-1,:]))))
 
		#Tpseudo[k,:] =a*(np.transpose(b[k,:]) - np.matmul(lower[0,k-1,:,:],(np.transpose(Tpseudo[k-1,:]))))
		Tpseudo[k,:]=np.matmul((np.linalg.inv(diagonal[0,k,:,:])),(np.transpose(b[k,:]))-np.matmul(lower[0,k-1,:,:],np.transpose(Tpseudo[k-1,:])))	

	T[N_layers-1,:]=Tpseudo[N_layers-1,:]
	for k in range (N_layers-2,-1,-1):
		T[k,:]=np.transpose(Tpseudo[k,:])- np.matmul(X[0,k,:,:],np.transpose(T[k+1,:]))	
	#res = residual(T,P,N_layers,r,ehin,ecin)
			
	
	#print("T")
	return T
	#print("Tpseudo")
	#print(Tpseudo)
	#print("diagonal")
	#print(diagonal)
	#print("lower")
	#print(lower)
	#print("m")
	#print(m)
	




