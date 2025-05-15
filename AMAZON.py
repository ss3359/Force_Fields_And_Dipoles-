import math 
import numpy as np 
import matplotlib.pyplot as plt 

#Constants
PI=np.pi
mu0=(4*PI)*10e-7 # Permeability of Free space N/A^2

#Unit Vectors
i=np.array([1,0,0])
j=np.array([0,1,0])
k=np.array([0,0,1])

def cross(U,V):
   

   """
      +    -    +
  |   i    j    k  |
  |  u0   u1   u2  | = i*(u1v2-u2v1)-j*(u0*v2-u2*v0)+k*(u0v1-u1v0)
  |  v0   v1   v2  |
   """
   UcrossV=((U[1]*V[2]-U[2]*V[1]))*i-((U[0]*V[2]-U[2]*V[0]))*j+(U[0]*V[1]-V[0]*U[1])*k

   return (UcrossV)

class MagneticField: 
   def __init__(self,q=1.0):

      self.q=q

   def B(self,R,M):
     
      denom=np.linalg.norm(R)**3
      if(abs(denom)<1e-7): 
         return np.array([1e-7,1e-7,1e-7])
      numerator=3*np.dot(M,R)*R-M*(np.linalg.norm(R)**2)
      return (mu0/(4*PI))*(numerator/denom)
   
   def F(self,B,V):
      q=self.q
      return q*(cross(B,V)) 
      
      
   def PlotMagneticField(self): 
      

     x=np.linspace(-5,5,15)
     y=np.linspace(-5,5,15)
     z=np.linspace(-5,5,15)
     X,Y,Z=np.meshgrid(x,y,z)
     
     Bx=np.zeros_like(X)
     By=np.zeros_like(Y)
     Bz=np.zeros_like(Z)
     Fx=np.zeros_like(X)
     Fy=np.zeros_like(Y)
     Fz=np.zeros_like(Z)
     V=np.array([1,0,0])
     for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
            
                R=np.array([X[i,j,k],Y[i,j,k], Z[i,j,k] ])
                if np.linalg.norm(R)<1e-3: 
                    continue
                M=np.array([0,1,0])
                B=self.B(R,M)
                F=self.F(B,V)
                Bx[i,j,k]=B[0]
                By[i,j,k]=B[1]
                Bz[i,j,k]=B[2]
                Fx[i,j,k]=F[0]
                Fy[i,j,k]=F[1]
                Fz[i,j,k]=F[2]
    
     B_norm=np.sqrt((Bx**2)+(By**2)+(Bz**2))
     Bx/=B_norm+100
     By/=B_norm+100
     Bz/=B_norm+100

     fig=plt.figure(figsize=(10,5))
     fig.suptitle("Magnetic Field and Lorentz Force Field")

     ax1=fig.add_subplot(121,projection='3d')
     ax1.quiver(X,Y,Z,Bx,By,Bz,color="red",length=0.5,normalize=True)
     ax1.set_title("Magnetic Field Of A Diplole")
     ax1.set_xlabel("x (meters)")
     ax1.set_ylabel("y (meters)")
     ax1.set_zlabel("z (meters)")
     ax1.set_box_aspect([1,1,1])
     


     ax2=fig.add_subplot(122,projection='3d')
     ax2.quiver(X,Y,Z,Fx,Fy,Fz,length=0.5,color="green",normalize=True)
     ax2.set_title("Plot Of A Lorentz Force Field")
     ax2.set_xlabel("x (meters)")
     ax2.set_ylabel("y (meters)")
     ax2.set_zlabel("z (meters)")
     ax2.set_box_aspect([1,1,1])

     plt.tight_layout()
     plt.show()



          


      

def main():

  
    MF=MagneticField()
    MF.PlotMagneticField()

main()