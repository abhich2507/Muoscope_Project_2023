import csv
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors
# from scipy.stats import gaussian_kde
# import seaborn as sns
# from sklearn.metrics import confusion_matrix
import torch
from torch import Tensor
#from sklearn.linear_model import LinearRegression

# from skspatial.objects import Line
# from skspatial.objects import Points
# from skspatial.plotting import plot_3d

#from torch import linalg

def compute_xy_min_max(X_discrete:Tensor,Y_discrete:Tensor)->Tensor:

    mask_x = X_discrete[0,:]>X_discrete[1,:]
    mask_y = Y_discrete[0,:]>Y_discrete[1,:]

    X_max = torch.where(mask_x,X_discrete[:-1,:],X_discrete[1:,:])
    X_min = torch.where(mask_x,X_discrete[1:,:],X_discrete[:-1,:])

    Y_max = torch.where(mask_y,Y_discrete[:-1,:],Y_discrete[1:,:])
    Y_min = torch.where(mask_y,Y_discrete[1:,:],Y_discrete[:-1,:])
    
    return X_max,X_min,Y_max,Y_min

def csv_to_tensor(filename: str, Nplane: int = 6,nb_event: int=None) -> Tensor:
    '''
        Create a tensor containing hits with shape [Nevent][Nplane][3]
    '''
    data = pd.read_csv(filename)
    Nevent = len(data)
    Hits = torch.zeros(3,Nplane,Nevent,dtype=torch.double) #[x,y,z][6 planes][N events]
    for i in range(Nplane):
#         Hits[0][i] = torch.tensor(data['X'+str(i)].values,dtype = torch.double)+1000
#         Hits[1][i] = torch.tensor(data['Y'+str(i)].values,dtype = torch.double)+1000
        Hits[0][i] = torch.tensor(data['X'+str(i)].values,dtype = torch.double)
        Hits[1][i] = torch.tensor(data['Y'+str(i)].values,dtype = torch.double)
        Hits[2][i] = torch.tensor(data['Z'+str(i)].values,dtype = torch.double)
       
    # Fixing the z hit coordinate to match detector placement
#     Hits[2,0,:]=1.
#     Hits[2,1,:]=101.
#     Hits[2,2,:]=201.
#     Hits[2,3,:]=801.
#     Hits[2,4,:]=901.
#     Hits[2,5,:]=1001.
    
    if (nb_event is not None):
            Hits = Hits[:,:,:nb_event]

    return Hits[:,:,:]

def norm(vect: Tensor)->Tensor:
    '''
    Compute vector norm. Vector must have shape [Nevent][3] (x,y,z) or [3].
    '''
    x,y,z = None,None,None
    if((vect.size(dim=-1)==3) & (len(vect)==3)):
        x,y,z = vect[0],vect[1],vect[2]
    else:
        x,y,z = vect[:,0],vect[:,1],vect[:,2]
    return torch.sqrt(x**2+y**2+z**2)



class VOI:
    
    def __init__(self, position:list, dimension:list, voxel_width:float=10., extra_gap = 50.):
        '''
        position = [x,y,z] in mm
        dimension = [dx,dy,dz] in mm 
        Voxel width = 10 mm default
        '''
        #extend VOI dimension 
        self.Extra = extra_gap #mm
        
        #VOI dimensions
        self.Center = position
        
        self.dX = dimension[0]
        self.dY = dimension[1]
        self.dZ = dimension[2]
        
        #X
        self.Xmin = position[0]-dimension[0]/2-self.Extra
        self.Xmax = position[0]+dimension[0]/2+self.Extra
        #Y
        self.Ymin = position[1]-dimension[1]/2-self.Extra
        self.Ymax = position[1]+dimension[1]/2+self.Extra
        #Z
        self.Zmin = position[2]-dimension[2]/2-self.Extra
        self.Zmax = position[2]+dimension[2]/2+self.Extra
             
        
        #VOXELIZATION
        self.Voxel_width = voxel_width
        
        self.Nvoxel_x = None
        self.Nvoxel_y = None
        self.Nvoxel_z = None
        
        self.Voxels_centers = None
        self.Voxels_edges = None
        
                
        
    def Compute_N_voxels(self)->None:
        self.Nvoxel_x=(self.dX+2*self.Extra)/self.Voxel_width
        self.Nvoxel_y=(self.dY+2*self.Extra)/self.Voxel_width
        self.Nvoxel_z=(self.dZ+2*self.Extra)/self.Voxel_width
        
        if((self.Nvoxel_x%1!=0)|(self.Nvoxel_y%1!=0)|(self.Nvoxel_z%1!=0)):
            print('Caution voxel size doesnt match VOI size')
             
        self.Nvoxel_x=int(self.Nvoxel_x)
        self.Nvoxel_y=int(self.Nvoxel_y)
        self.Nvoxel_z=int(self.Nvoxel_z)
        
        
    def Compute_voxel_centers(self,x_min_: float, x_max_:float,Nvoxel_:int)->torch.tensor:
        
        '''
        x_min,max border of the volume of interset for a given coordinate
                
        return voxels centers along given coordinate
        '''
        xs_ = torch.linspace(x_min_,x_max_,Nvoxel_+1)
        xs_ += self.Voxel_width/2
        return xs_[:-1]
        
    def Generate_voxels(self)->None:
                
        Voxels_centers = torch.zeros((self.Nvoxel_x,self.Nvoxel_y,self.Nvoxel_z,3),dtype=torch.double)


        xs_ = self.Compute_voxel_centers(x_min_=self.Xmin, x_max_=self.Xmax,Nvoxel_= self.Nvoxel_x)
        ys_ = self.Compute_voxel_centers(x_min_=self.Ymin, x_max_=self.Ymax,Nvoxel_= self.Nvoxel_y)
        zs_ = self.Compute_voxel_centers(x_min_=self.Zmin, x_max_=self.Zmax,Nvoxel_= self.Nvoxel_z)

        xs,ys,zs = xs_[:,None,None],ys_[:,None,None],zs_[:,None,None]
        
        xs = xs.expand(len(xs_),len(ys_),len(zs_))
        ys = torch.transpose(ys.expand(len(ys_),len(xs_),len(zs_)),0,1)
        zs = torch.transpose(zs.expand(len(zs_),len(ys_),len(xs_)),0,2)

        Voxels_centers[:,:,:,0],Voxels_centers[:,:,:,1],Voxels_centers[:,:,:,2] = xs,ys,zs
        
        self.Voxels_centers = Voxels_centers
        
        Voxels_edges = torch.zeros((self.Nvoxel_x,self.Nvoxel_y,self.Nvoxel_z,2,3))
        Voxels_edges[:,:,:,0,:] = Voxels_centers-self.Voxel_width/2
        Voxels_edges[:,:,:,1,:] = Voxels_centers+self.Voxel_width/2
        self.Voxels_edges = Voxels_edges
        

class Detector:
    
    def __init__(self, z_position_:list, dimension_:list, xy_position_:list):
        
        '''
        z_position = [z0,...,Zn]
        dimension = [dX,dY]
        xy_position = [x,y]
        '''
        self.Z = z_position_
        self.X = xy_position_[0]
        self.Y = xy_position_[1]
        
        self.dX = dimension_[0]
        self.dY = dimension_[1]
        self.dZ = None
        
        self.Xmin= self.X - self.dX/2
        self.Ymin= self.Y - self.dY/2
        self.Zmin= min(self.Z)
        
        self.Xmax= self.X + self.dX/2
        self.Ymax= self.Y + self.dY/2
        self.Zmax= max(self.Z)
            
    
class Tracks(VOI,Detector):
    
    def __init__(self, Hits:torch.tensor, aVOI, aDetector)->None:
        
        #HITS
        self.hits = Hits
        self.X = Hits[0,:,:]
        self.Y = Hits[1,:,:]
        self.Z = Hits[2,:,:]
        self.Nevent = Hits.shape[2]
        self.events = torch.linspace(0,self.Nevent-1,self.Nevent)
        
        #TRACKS
        self.tracks_up = torch.zeros(self.Nevent,3,dtype=torch.double)
        self.tracks_down = torch.zeros(self.Nevent,3,dtype=torch.double)
        
        self.t_point_up = None
        self.t_point_down = None

        self.phis_up = None
        self.phis_down = None
        
        self.thetas_up = None
        self.thetas_down = None
        
        self.dthetas = None
        self.dphis = None
        
        #SCATTERING
        self.theta_msc = None
        self.theta_msc_y = None
        self.theta_msc_x = None
        
        
        #muon hit VOI mask
        self.mask = None
        
        #tracking computation
        self.track_in_discrete = None
        self.track_out_discrete = None
        
        #VOXELIZATION
        self.Voxels_centers = aVOI.Voxels_centers
        self.Voxels_edges = aVOI.Voxels_edges
        
        self.VOI_Xmin = aVOI.Xmin
        self.VOI_Xmax = aVOI.Xmax
        
        self.VOI_Ymin = aVOI.Ymin
        self.VOI_Ymax = aVOI.Ymax
        
        self.VOI_Zmin = aVOI.Zmin
        self.VOI_Zmax = aVOI.Zmax
        
        self.Voxel_width = aVOI.Voxel_width
        
        self.Nvoxel_x = aVOI.Nvoxel_x
        self.Nvoxel_y = aVOI.Nvoxel_y
        self.Nvoxel_z = aVOI.Nvoxel_z
        
        self.Voxels_hit_indices_up = None
        self.Voxels_hit_indices_down = None
        
        self.Voxels_hit_indices_up = None
        self.Voxels_hit_indices_down = None
        
        #Detector boundaries (use for ploting function)
        self.Xmin = aDetector.Xmin
        self.Xmax = aDetector.Xmax
        
        self.Ymin = aDetector.Ymin
        self.Ymax = aDetector.Ymax
        
        self.Zmin = aDetector.Zmin
        self.Zmax = aDetector.Zmax
        
        self.Nplane = Hits.shape[1]
        
        #POCA
        self.POCA_points = None
        self.POCA_mask = None
        
    def Compute_tracks_3planes(self)->None:
        
        '''
            Compute upper and lower track from hits[3][Nplane][Nevent]
        '''
        
        from skspatial.objects import Line
        from skspatial.objects import Points
        from skspatial.plotting import plot_3d
        
        track_up,track_down = torch.zeros(self.Nevent,3,dtype=torch.double),torch.zeros(self.Nevent,3,dtype=torch.double)
        point_up,point_down = torch.zeros(self.Nevent,3,dtype=torch.double),torch.zeros(self.Nevent,3,dtype=torch.double)
        
        
        for i in range(self.Nevent):
            up_hits = torch.transpose(self.hits[:,:3,i],0,1)
            down_hits = torch.transpose(self.hits[:,3:,i],0,1)
            
            pointsUP, pointsDOWN = Points(up_hits),Points(down_hits)
            fit_UP, fit_DOWN = Line.best_fit(pointsUP),Line.best_fit(pointsDOWN)
            track_up[i], track_down[i] = torch.tensor(fit_UP.direction),torch.tensor(fit_DOWN.direction)
            point_up[i], point_down[i] = torch.tensor(fit_UP.point),torch.tensor(fit_DOWN.point)
        
        self.tracks_up,self.tracks_down = track_up,track_down
        self.t_point_up,self.t_point_down = point_up,point_down
        
    def Compute_tracks_2planes(self)->None:
        
        track_out,track_in = torch.zeros((self.Nevent,3)),torch.zeros((self.Nevent,3))

        track_in[:,0],track_in[:,1],track_in[:,2] = self.X[1]-self.X[0],self.Y[1]-self.Y[0],self.Z[1]-self.Z[0]
        track_in[:,0]/=(torch.sqrt(track_in[:,0]**2+track_in[:,1]**2+track_in[:,2]**2))
        track_in[:,1]/=(torch.sqrt(track_in[:,0]**2+track_in[:,1]**2+track_in[:,2]**2))
        track_in[:,2]/=(torch.sqrt(track_in[:,0]**2+track_in[:,1]**2+track_in[:,2]**2))

        track_out[:,0],track_out[:,1],track_out[:,2] = self.X[3]-self.X[2],self.Y[3]-self.Y[2],self.Z[3]-self.Z[2]
        track_out[:,0]/=(torch.sqrt(track_out[:,0]**2+track_out[:,1]**2+track_out[:,2]**2))
        track_out[:,1]/=(torch.sqrt(track_out[:,0]**2+track_out[:,1]**2+track_out[:,2]**2))
        track_out[:,2]/=(torch.sqrt(track_out[:,0]**2+track_out[:,1]**2+track_out[:,2]**2))


        t_point_up,t_point_down = torch.zeros((self.Nevent,3)),torch.zeros((self.Nevent,3))
        t_point_up[:,0],t_point_up[:,1],t_point_up[:,2] = self.X[1],self.Y[1],self.Z[1]
        t_point_down[:,0],t_point_down[:,1],t_point_down[:,2] = self.X[2],self.Y[2],self.Z[2]
        
        self.tracks_up,self.tracks_down = track_in,track_out
        self.t_point_up,self.t_point_down = t_point_up,t_point_down
        
        
    def Compute_scattering_angle(self)->None:
        '''
        Compute spatial scattering angle theta_msc, and x,y projected scattering angle theta_msc_x, theta_msc_y, for each event.
        INPUT: tracks_up, tracks_down
        OUTPUT: theta_msc_x,theta_msc_y,theta_msc
        '''
        
        xup, yup, zup = self.tracks_up[:,0],self.tracks_up[:,1],self.tracks_up[:,2]
        xdown, ydown, zdown = self.tracks_down[:,0],self.tracks_down[:,1],self.tracks_down[:,2]
        
        self.theta_msc = torch.acos((xup*xdown + yup*ydown + zup*zdown)/(norm(self.tracks_up)*norm(self.tracks_down)))
        self.theta_msc_x = torch.acos((xup*xdown + zup*zdown)/((torch.sqrt(xup**2+zup**2))*(torch.sqrt(xdown**2+zdown**2))))
        self.theta_msc_y = torch.acos((yup*ydown + zup*zdown)/((torch.sqrt(yup**2+zup**2))*(torch.sqrt(ydown**2+zdown**2))))
        # scattering amgles close to zero computed as NaN
        #Nans are changed to 0 scattering angles
        self.theta_msc[torch.isnan(self.theta_msc)]=0.0001
        #Random error during track fitting occurs and change z track coordinate to -z
        #which gives an angle close to pi instead of close to zero
        #We change these angles to pi - theta
        self.theta_msc[self.theta_msc==0] = 0.0001
        self.theta_msc[self.theta_msc>2] = torch.abs(torch.pi-self.theta_msc[self.theta_msc>2])
        self.theta_msc_x[self.theta_msc_x>2] = torch.abs(torch.pi-self.theta_msc_x[self.theta_msc_x>2])
        self.theta_msc_y[self.theta_msc_y>2] = torch.abs(torch.pi-self.theta_msc_y[self.theta_msc_y>2])
        self.theta_msc_x[torch.isnan(self.theta_msc_x)]=0.0001
        self.theta_msc_y[torch.isnan(self.theta_msc_y)]=0.0001
        self.theta_msc_x = torch.where(((self.X[-1,:]-self.X[0,:])>0),self.theta_msc_x,-self.theta_msc_x)
        
        
        
    def Compute_theta_phi(self)->None:
        
        '''
        Compute zenith (theta) and azimuthal (phi) angle of the incoming track
        INPUT: tracks_up, tracks_down
        OUTPUT: thetas, phis
        '''
        
        for i in [0,1]:
            if(i==0):
                track = self.tracks_up
            if(i==1):
                track = self.tracks_down
                
            x, y, z = track[:,0],track[:,1],track[:,2]
            r = torch.sqrt(x**2 +y**2+z**2)

            theta = torch.acos(z/r)
            phi = torch.tensor(np.shape(track[:,0]),dtype=torch.double)

            mask1 = (x>0)
            mask2 = (x<0)&(y>=0)
            mask3 = (x<0)&(y<0)
            mask4 = (x==0)&(y>0)
            mask5 = (x==0)&(y<0)

            phi = np.where(mask1,torch.atan(y/x),phi)
            phi = np.where(mask2,torch.atan(y/x) + math.pi,phi)
            phi = np.where(mask3,torch.atan(y/x) - math.pi,phi)
            phi = np.where(mask4,math.pi/2,phi)
            phi = np.where(mask5,-math.pi/2,phi)
            
            if(i==0):
                self.thetas_up = theta
                self.phis_up = torch.tensor(phi)
            if(i==1):
                self.thetas_down = theta
                self.phis_down = torch.tensor(phi)
                
        
    def Compute_dtheta_dphi(self)->None:
        self.dphis = torch.abs(self.phis_up - self.phis_down)
        self.dthetas = torch.abs(self.thetas_up - self.thetas_down)

        
        
    def Compute_POCA_points(self)->None:
        from numpy import array, cross
        '''
            @MISC {3334866,
        TITLE = {Closest points between two lines},
        AUTHOR = {Brian (https://math.stackexchange.com/users/72614/brian)},
        HOWPUBLISHED = {Mathematics Stack Exchange},
        NOTE = {URL:https://math.stackexchange.com/q/3334866 (version: 2019-08-26)},
        EPRINT = {https://math.stackexchange.com/q/3334866},
        URL = {https://math.stackexchange.com/q/3334866}
}
        
        Compute POCA points. 
        INPUT: incoming and outgoing reconstructed tracks (Tensor[Nevent][3], Tensor[Nevent][3])
        OUTPUT: POCA points (Tensor[Nevent][3])
        
        Given 2 lines L1, L2 aka incoming and outgoing tracks with parametric equation:
        L1 = P1 + t*V1
        
        1- A segment of shortest length between two 3D lines L1 L2 is perpendicular to both lines (if L1 L2 are neither parallele or in the same plane). One must compute V3, vector perpendicular to L1 and L2
        2- Search for points where L3 = P1 + t1*V1 +t3*V3 crosses L2. One must find t1 and t2 for which:
        L3 = P1 + t1*V1 +t3*V3 = P2 + t2*V2
        
        3- Then POCA location M is the middle of the segment Q1-Q2 where Q1,2 = P1,2 +t1,2*V1,2
        '''
        P1, P2 = self.t_point_up[:], self.t_point_down[:]
        V1, V2 = self.tracks_up[:], self.tracks_down[:]
        
        V3 = torch.tensor(cross(V2,V1))
        
        RES = P2 - P1
        LES = torch.transpose(torch.stack([V1,-V2,V3]),0,1)
        LES = torch.transpose(LES,-1,1)
        if RES!=0:
            (ts = torch.linalg.solve(LES,RES))

        t1 = torch.stack([ts[:,0],ts[:,0],ts[:,0]],-1)
        t2 = torch.stack([ts[:,1],ts[:,1],ts[:,1]],-1)

        Q1s,Q2s = P1+t1*V1, P2+t2*V2
        M = (Q2s-Q1s)/2+Q1s
        
        self.POCA_points = M
        
    def POCA_inside_VOI(self)->None:
        
        '''
        Only keep events for which POCA point is located INSIDE the VOI.
        
        INPUT: self.POCA_points, the POCA locations
        
        OUTPUT: self.POCA_mask, a 1D mask to be applied on events 
        self.POCA_mask rejects events which do not reach the VOI & POCA points outside the VOI
        '''
        #POCA point has to be inside VOI
        mask_x = (self.POCA_points[:,0]>=self.VOI_Xmin) & (self.POCA_points[:,0]<=self.VOI_Xmax)
        mask_y = (self.POCA_points[:,1]>=self.VOI_Ymin) & (self.POCA_points[:,1]<=self.VOI_Ymax)
        mask_z = (self.POCA_points[:,2]>=self.VOI_Zmin) & (self.POCA_points[:,2]<=self.VOI_Zmax)

        #interaction point inside VOI condition
        mask_VOI = mask_x & mask_y & mask_z

        #POCA points inside VOI + track hit VOI
        mask = mask_VOI & self.mask
        self.POCA_mask = mask
        
        
    def assign_voxel_POCA(self)->list:

        '''
        - Assign a voxel to a POCA point for each event
        
        - INPUT: self.POCA_points, the POCA locations
        - STATIC: Voxel edges positions
        - OUTPUT: INDICES, list of 1D tensor with size 3, containing the X,Y and Z index of the hit voxel
        len(INDICES) = len(self.events[self.POCA_mask])
        
        '''
        
        POCA_in = self.POCA_points[self.POCA_mask]
        INDICES=[]
        for i in range(len(POCA_in)):

            mask_vox_x1 = (POCA_in[i,0]>=self.Voxels_edges[:,:,:,0,0]) 
            mask_vox_x2 = (POCA_in[i,0]<=self.Voxels_edges[:,:,:,1,0])

            mask_vox_y1 = (POCA_in[i,1]>=self.Voxels_edges[:,:,:,0,1]) 
            mask_vox_y2 = (POCA_in[i,1]<=self.Voxels_edges[:,:,:,1,1])

            mask_vox_z1 = (POCA_in[i,2]>=self.Voxels_edges[:,:,:,0,2]) 
            mask_vox_z2 = (POCA_in[i,2]<=self.Voxels_edges[:,:,:,1,2])

            mask_x_ = mask_vox_x1 & mask_vox_x2
            mask_y_ = mask_vox_y1 & mask_vox_y2    
            mask_z_ = mask_vox_z1 & mask_vox_z2

            mask_vox = mask_x_ & mask_y_ & mask_z_

            indices = (mask_vox==True).nonzero()
            INDICES.append(indices[0])
        return INDICES
        
        
    def Compute_POCA_scores(self):
        #identify interation point voxels

        #POCA point coordinate has to be inside VOI
        mask_x = (self.POCA_points[:,0]>=self.VOI_Xmin) & (self.POCA_points[:,0]<=self.VOI_Xmax)
        mask_y = (self.POCA_points[:,1]>=self.VOI_Ymin) & (self.POCA_points[:,1]<=self.VOI_Ymax)
        mask_z = (self.POCA_points[:,2]>=self.VOI_Zmin) & (self.POCA_points[:,2]<=self.VOI_Zmax)

        #interaction point inside VOI condition
        mask_VOI = mask_x & mask_y & mask_z

        #POCA points inside VOI
        POCA_in = self.POCA_points[mask_VOI]

        #ratio of POCA points inside volume


        INDICES = []
        Voxels_scores = torch.zeros((self.Nvoxel_x,self.Nvoxel_y,self.Nvoxel_z,2))
        for i in range(len(self.events[mask_VOI])):

            mask_vox_x1 = (POCA_in[i,0]>=self.Voxels_edges[:,:,:,0,0]) 
            mask_vox_x2 = (POCA_in[i,0]<=self.Voxels_edges[:,:,:,1,0])

            mask_vox_y1 = (POCA_in[i,1]>=self.Voxels_edges[:,:,:,0,1]) 
            mask_vox_y2 = (POCA_in[i,1]<=self.Voxels_edges[:,:,:,1,1])

            mask_vox_z1 = (POCA_in[i,2]>=self.Voxels_edges[:,:,:,0,2]) 
            mask_vox_z2 = (POCA_in[i,2]<=self.Voxels_edges[:,:,:,1,2])

            mask_x = mask_vox_x1 & mask_vox_x2
            mask_y = mask_vox_y1 & mask_vox_y2    
            mask_z = mask_vox_z1 & mask_vox_z2

            mask = mask_x & mask_y & mask_z

            indices = (mask==True).nonzero()
            ix,iy,iz = indices[0,0],indices[0,1],indices[0,2]
            INDICES.append(indices)
            if(self.theta_msc[mask_VOI][i]*(180/math.pi)>1):
                Voxels_scores[ix,iy,iz,0]+=(self.theta_msc[mask_VOI][i]**2/(self.Voxel_width/self.thetas_up[mask_VOI][i]))
                Voxels_scores[ix,iy,iz,1]+=1

        Voxels_scores[:,:,:,0][Voxels_scores[:,:,:,0]!=0]=Voxels_scores[:,:,:,0][Voxels_scores[:,:,:,0]!=0]/Voxels_scores[:,:,:,1][Voxels_scores[:,:,:,0]!=0]
        self.Voxels_scores = Voxels_scores
        
        
    def Compute_discrete_tracks(self)->None:
        
        '''
        Compute x,y,z position at Zmax and Zmin of each voxel layer (for incoming and outgoing tracks)
        INPUT: - Zmin Zmax of each voxel layer
               - theta_xy_in_x,y, theta_xy_out_x,y: zentih angle in x,y for incoming and outgoing track
               - xyz_in_x,y,z, xyz_out_x,y,z: muon position entering and exiting the Volume of Interest
               
        OUTPUT: - track_in_discrete, track_out_discrete: x,y,z position at Zmax and Zmin of each voxel layer (for incoming and outgoing tracks), size = [coordinate,Nlayer_along_Z + 1, Nevents] ([3,7,9999])
        '''

        Z_discrete = torch.linspace(torch.min(self.Voxels_edges[:,:,:,:,2]).item(),
                                    torch.max(self.Voxels_edges[:,:,:,:,2]).item(),
                                    self.Nvoxel_z+1)

        Z_discrete.unsqueeze_(1)
        Z_discrete = Z_discrete.expand(len(Z_discrete),len(self.events))
        
        x_track_in = self.xyz_in_x + (self.xyz_in_z[0]-Z_discrete)*torch.tan(self.theta_xy_in_x)
        x_track_out = self.xyz_out_x + (self.xyz_out_z[0]-Z_discrete)*torch.tan(self.theta_xy_out_x)

        y_track_in = self.xyz_in_y + (self.xyz_in_z[0]-Z_discrete)*torch.tan(self.theta_xy_in_y)
        y_track_out = self.xyz_out_y + (self.xyz_out_z[0]-Z_discrete)*torch.tan(self.theta_xy_out_y)

        size = list(x_track_in.size())
        size.insert(0,3)

        self.track_in_discrete = torch.ones(size)
        self.track_out_discrete = torch.ones(size)

        self.track_in_discrete[0],self.track_in_discrete[1],self.track_in_discrete[2] = x_track_in,y_track_in,Z_discrete
        self.track_out_discrete[0],self.track_out_discrete[1],self.track_out_discrete[2] = x_track_out,y_track_out,Z_discrete
        
    def find_hit_voxel(self):
        
        '''
        Identify voxels compatibles with both incoming and outgoing track.
        The indicies of hit voxels are then stored in self.Voxel_hit_indices
        INPUT: Incoming, outgoing track, Voxels_edges, Voxel_width
        OUTPUT: self.Voxels_hit_indices_up (voxel hit by incoming track)
                self.Voxels_hit_indices_down (voxel hit by outgoing track)
                self.Voxels_hit_indices (voxels hit by both)
        '''
        
        #Compute x,y track position at z entry and exit of each voxel 
        t=torch.linspace(0,2600,5000)#to be modified for other VOI dimensions
        zmin, zmax, zrange = 200.,800.,600. #mm
        z = torch.linspace(self.VOI_Zmin,self.VOI_Zmax,int(self.Nvoxel_z)+1)
        z.unsqueeze_(1)
        z = z.expand(len(z),len(self.thetas_up[:]))

        x_up = self.t_point_up[:,0] + (z-self.t_point_up[0,2])/torch.cos(self.thetas_up[:])*self.tracks_up[:,0]
        y_up = self.t_point_up[:,1] + (z-self.t_point_up[0,2])/torch.cos(self.thetas_up[:])*self.tracks_up[:,1]
        
        x_down = self.t_point_down[:,0] + (z-self.t_point_down[0,2])/torch.cos(self.thetas_down[:])*self.tracks_down[:,0]
        y_down = self.t_point_down[:,1] + (z-self.t_point_down[0,2])/torch.cos(self.thetas_down[:])*self.tracks_down[:,1]

        #Find if track goes from low x,y to large x,y
        mask_min_xup = x_up[0,:]>x_up[1,:]
        mask_min_yup = y_up[1,:]>y_up[0,:]
        
        mask_min_xdown = x_down[0,:]>x_down[1,:]
        mask_min_ydown = y_down[1,:]>y_down[0,:]

        empty = torch.zeros((np.shape(x_up[1:,:])),dtype=torch.double)
        #x_min,max at the entry/exit of each voxel
        x_max_up = torch.where(mask_min_xup,x_up[:-1,:],x_up[1:,:])
        x_min_up = torch.where(mask_min_xup,x_up[1:,:],x_up[:-1,:])

        y_min_up = torch.where(mask_min_yup,y_up[:-1,:],y_up[1:,:])
        y_max_up = torch.where(mask_min_yup,y_up[1:,:],y_up[:-1,:])
        
        x_max_down = torch.where(mask_min_xdown,x_down[:-1,:],x_down[1:,:])
        x_min_down = torch.where(mask_min_xdown,x_down[1:,:],x_down[:-1,:])

        y_min_down = torch.where(mask_min_ydown,y_down[:-1,:],y_down[1:,:])
        y_max_down = torch.where(mask_min_ydown,y_down[1:,:],y_down[:-1,:])

        list_indices_up=[]
        list_indices_down=[]
        list_indices = []
        
        #find voxel hit by the track
        for i in (self.events[self.mask]):# loop over events to be changed
            i=int(i.item())
            X_mask_up = (self.Voxels_edges[:,:,:,0,0]>=x_min_up[:,i]-self.Voxel_width) & (self.Voxels_edges[:,:,:,1,0]<=x_max_up[:,i]+self.Voxel_width)
            Y_mask_up = (self.Voxels_edges[:,:,:,0,1]>=y_min_up[:,i]-self.Voxel_width) & (self.Voxels_edges[:,:,:,1,1]<=y_max_up[:,i]+self.Voxel_width)
            mask_up = X_mask_up & Y_mask_up
            list_indices_up.append(((mask_up==True).nonzero()))
            
            X_mask_down = (self.Voxels_edges[:,:,:,0,0]>=x_min_down[:,i]-self.Voxel_width) & (self.Voxels_edges[:,:,:,1,0]<=x_max_down[:,i]+self.Voxel_width)
            Y_mask_down = (self.Voxels_edges[:,:,:,0,1]>=y_min_down[:,i]-self.Voxel_width) & (self.Voxels_edges[:,:,:,1,1]<=y_max_down[:,i]+self.Voxel_width)
            mask_down = X_mask_down & Y_mask_down
            list_indices_down.append(((mask_down==True).nonzero()))
            
            list_indices.append(((mask_up & mask_down)==True).nonzero())
            
            
        #return a list containing indices of hit voxels for each event
        self.Voxels_hit_indices_up = list_indices_up
        self.Voxels_hit_indices_down = list_indices_down
        self.Voxels_hit_indices = list_indices
        
    def mask_hit_voxel(self)->None:
        
        '''
        Create a mask which reject muons who do not reach the VOI.
        
        - Compute track x,y positions for each voxel layer alomg z axis
        - Check if among x,y positions, at least one x AND one y is located inside the VOI
        - Assign the mask as a class attribute
        
        '''
        #Compute track position
        z = torch.linspace(self.VOI_Zmin,self.VOI_Zmax,int(self.Nvoxel_z)+1)
        z.unsqueeze_(1)
        z = z.expand(len(z),len(self.thetas_up[:]))

        x_up = self.t_point_up[:,0] + (z-100.)/torch.cos(self.thetas_up[:])*self.tracks_up[:,0]
        y_up = self.t_point_up[:,1] + (z-100.)/torch.cos(self.thetas_up[:])*self.tracks_up[:,1]
        
        #check if inside the VOI
        mask_x = (x_up[:,:]>=self.VOI_Xmin) & (x_up[:,:]<=self.VOI_Xmax) 
        mask_y = (y_up[:,:]>=self.VOI_Ymin) & (y_up[:,:]<=self.VOI_Ymax) 
        mask = (mask_x & mask_y)
        mask = torch.where(mask,1,0)
        mask = torch.sum(mask,0)
        mask = torch.where(mask>0,True,False)

        self.mask = mask
        
    def plot_tracks(self, Event: int, filename:str=None)->None:
        
        from matplotlib.patches import Rectangle
        print('Event = ', Event)
        print('scattering angle = ', self.theta_msc[Event].item())
        
#         fig = plt.figure(figsize = (10,8))
#         ax = plt.axes(projection='3d')
#         ax.set_xlim3d(self.Xmin, self.Xmax)
#         ax.set_ylim3d(self.Ymin, self.Ymax)
#         ax.set_zlim3d(self.Zmin, self.Zmax)
#         ## Plot hits
#         ax.scatter(self.X[:int(self.Nplane/2),Event].tolist(),
#                    self.Y[:int(self.Nplane/2),Event].tolist(),
#                    self.Z[:int(self.Nplane/2),Event].tolist())
        
#         ax.scatter(self.X[int(self.Nplane/2):,Event].tolist(),
#                    self.Y[int(self.Nplane/2):,Event].tolist(),
#                    self.Z[int(self.Nplane/2):,Event].tolist())
        
#         ## Plot reconstructed track
        cd_up, cd_down = self.tracks_up[Event], self.tracks_down[Event]
        p_up, p_down = self.t_point_up[Event], self.t_point_down[Event]
        
        t = np.linspace(-8000,8000,100)
        x0,y0,z0 = cd_up[0]*t+p_up[0], cd_up[1]*t+p_up[1], cd_up[2]*t+p_up[2]
        x1,y1,z1 = cd_down[0]*t+p_down[0], cd_down[1]*t+p_down[1], cd_down[2]*t+p_down[2]
        
#         ax.plot(x0,y0,z0,label = 'incoming track')
#         ax.plot(x1,y1,z1, label = 'outgoing track')
#         if(self.POCA_points is not None):
#             ax.scatter(self.POCA_points[Event,0].item(),
#                        self.POCA_points[Event,1].item(),
#                        self.POCA_points[Event,2].item(), label = 'POCA point')
#         ax.legend()
#         plt.show()
        
        fig, ax = plt.subplots(1, 2, figsize=(20, 7))
        ax[0].set_xlabel('X [mm]')
        ax[0].set_ylabel('Z [mm]')
        ax[0].set_title('reconstructed track in XZ plane')
        #plot xz frame
        ax[0].scatter(x0,z0,label = 'track up',alpha = 0.3,s=10)
        ax[0].scatter(x1,z1,label = 'track down',alpha = 0.3,s=10)
#         ax[0].scatter(x0,z0,label = 'outgoing track',alpha = 0.3,s=10)
#         ax[0].scatter(x1,z1,label = 'incoming track',alpha = 0.3,s=10)
        ax[0].set_xlim(self.Xmin-50,self.Xmax+50)
        ax[0].set_ylim(self.Zmin-50,self.Zmax+50)
        
        for i in range(self.Nplane):
            ax[0].axhline(y=self.Z[i,Event], color='r', linestyle='-')
            ax[0].scatter(self.X[i,Event],self.Z[i,Event],marker = 'x',s=100,color='blue')

        
        ax[0].scatter(self.X[0,Event],self.Z[0,Event], label='recorded hit',marker = 'x',s=100,color='blue')

        if(self.POCA_points is not None):
            ax[0].scatter(self.POCA_points[Event,0], self.POCA_points[Event,2],marker='o',color='green',label='POCA point')
        
        ax[0].add_patch( Rectangle((self.VOI_Xmin,self.VOI_Zmin),
        self.VOI_Xmax-self.VOI_Xmin, self.VOI_Zmax-self.VOI_Zmin,
        fc ='none', 
        ec ='g',
        lw = 4,
        label = 'Volume of interest'))
        
        ax[0].legend()
        
        
        ax[1].set_xlabel('Y [mm]')
        ax[1].set_ylabel('Z [mm]')
        ax[1].set_title('reconstructed track in YZ plane')
        #plot xz frame
        ax[1].scatter(y0,z0,label = 'track up',alpha = 0.3,s=10)
        ax[1].scatter(y1,z1,label = 'track down',alpha = 0.3,s=10)
        ax[1].set_xlim(self.Ymin-50,self.Ymax+50)
        ax[1].set_ylim(self.Zmin-50,self.Zmax+50)
        
        for i in range(self.Nplane):
            ax[1].axhline(y=self.Z[i,Event], color='r', linestyle='-')
            ax[1].scatter(self.Y[i,Event],self.Z[i,Event],marker = 'x',s=100,color='blue')
            
        ax[1].scatter(self.Y[2,Event],self.Z[2,Event], label='recorded hit',marker = 'x',s=100,color='blue')

        ax[1].add_patch( Rectangle((self.VOI_Ymin,self.VOI_Zmin),
        self.VOI_Ymax-self.VOI_Ymin, self.VOI_Zmax-self.VOI_Zmin,
        fc ='none', 
        ec ='g',
        lw = 4,
        label = 'VOI'))
        
        if(self.POCA_points is not None):
            ax[1].scatter(self.POCA_points[Event,1], self.POCA_points[Event,2],marker='o',color='green',label='POCA point')
        
        ax[1].legend()
        plt.show()
        
        if(filename is not None):
            fig.savefig('Figures/'+filename+'.png')



    def plot_voxel_track(self, Event: int, up:bool=True,down:bool=True,filename:str=None):
        
        print('Event = ', Event, ', Scattering angle = ', self.theta_msc[Event].item())
        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        ax=axs.ravel()
        ax[0].set_title('Reconstructed voxels hits XZ projection')
        ax[0].set_xlabel('X [mm]')
        ax[0].set_ylabel('Z [mm]')
        
        #adjut axis limit
        ax[0].set_xlim(0,self.VOI_Xmax)
        ax[0].set_ylim(self.Z[0,Event]-4*self.Voxel_width,self.Z[5,Event]+4*self.Voxel_width)
        #plot detector layers
        ax[0].axhline(y=self.Z[2,Event], color='r', linestyle='-')
        ax[0].axhline(y=self.Z[3,Event], color='r', linestyle='-')
        ax[0].axhline(y=self.Z[1,Event], color='r', linestyle='-')
        ax[0].axhline(y=self.Z[4,Event], color='r', linestyle='-')
        ax[0].axhline(y=self.Z[0,Event], color='r', linestyle='-')
        ax[0].axhline(y=self.Z[5,Event], color='r', linestyle='-')
        #plot recorded hits
        ax[0].scatter(self.X[0,Event],self.Z[0,Event],marker = 'x',s=100,color='blue')
        ax[0].scatter(self.X[1,Event],self.Z[1,Event],marker = 'x',s=100,color='blue')
        ax[0].scatter(self.X[2,Event],self.Z[2,Event], label='recorded hit',marker = 'x',s=100,color='blue')
        ax[0].scatter(self.X[3,Event],self.Z[3,Event],marker = 'x',s=100,color='blue')
        ax[0].scatter(self.X[4,Event],self.Z[4,Event],marker = 'x',s=100,color='blue')
        ax[0].scatter(self.X[5,Event],self.Z[5,Event],marker = 'x',s=100,color='blue')
        #plot centers of hit voxels in XZ plane
        
        for i in range(len(self.Voxels_hit_indices_up[Event])):
            ix_up,iz_up = self.Voxels_hit_indices_up[Event][i,0].item(),self.Voxels_hit_indices_up[Event][i,2].item()
            if(up):
                ax[0].scatter(self.Voxels_centers[ix_up,0,iz_up,0],
                              self.Voxels_centers[ix_up,0,iz_up,2],
                              s=1,
                              c='green',
                              label = 'Voxels hit up',alpha=.3)
                
        for j in range(len(self.Voxels_hit_indices_down[Event])):
            ix_down = self.Voxels_hit_indices_down[Event][j,0].item()
            iz_down = self.Voxels_hit_indices_down[Event][j,2].item()
            
            if(down):
                ax[0].scatter(self.Voxels_centers[ix_down,0,iz_down,0],
                              self.Voxels_centers[ix_down,0,iz_down,2],
                              s=1,
                              c='orange',
                              label = 'Voxels hit down',alpha=.3)
            if(i==0):
                ax[0].legend()
        
        
        ax[1].set_title('Reconstructed voxels hits YZ projection')
        ax[1].set_xlabel('Y [mm]')
        ax[1].set_ylabel('Z [mm]')
        
        #adjut axis limit
        ax[1].set_xlim(0,self.VOI_Ymax)
        ax[1].set_ylim(self.Z[0,Event]-4*self.Voxel_width,self.Z[5,Event]+4*self.Voxel_width)
        #plot detector layers
        
        ax[1].axhline(y=self.Z[2,Event], color='r', linestyle='-')
        ax[1].axhline(y=self.Z[3,Event], color='r', linestyle='-')
        ax[1].axhline(y=self.Z[1,Event], color='r', linestyle='-')
        ax[1].axhline(y=self.Z[4,Event], color='r', linestyle='-')
        ax[1].axhline(y=self.Z[0,Event], color='r', linestyle='-')
        ax[1].axhline(y=self.Z[5,Event], color='r', linestyle='-')
        #plot recorded hits

        ax[1].scatter(self.Y[0,Event],self.Z[0,Event],marker = 'x',s=100,color='blue')
        ax[1].scatter(self.Y[1,Event],self.Z[1,Event],marker = 'x',s=100,color='blue')
        ax[1].scatter(self.Y[2,Event],self.Z[2,Event], label='recorded hit',marker = 'x',s=100,color='blue')
        ax[1].scatter(self.Y[3,Event],self.Z[3,Event],marker = 'x',s=100,color='blue')
        ax[1].scatter(self.Y[4,Event],self.Z[4,Event],marker = 'x',s=100,color='blue')
        ax[1].scatter(self.Y[5,Event],self.Z[5,Event],marker = 'x',s=100,color='blue')
        
        
        #plot centers of hit voxels in XZ plane
        for i in range(len(self.Voxels_hit_indices_up[Event])):
            iy_up,iz_up = self.Voxels_hit_indices_up[Event][i,1].item(),self.Voxels_hit_indices_up[Event][i,2].item()
            if(up):

                ax[1].scatter(self.Voxels_centers[0,iy_up,iz_up,1],
                              self.Voxels_centers[0,iy_up,iz_up,2],
                              s=1,
                              c='green',
                              label = 'Voxel hit up',alpha=.3)
    
        for j in range(len(self.Voxels_hit_indices_down[Event])):
        
            iy_down,iz_down = self.Voxels_hit_indices_down[Event][j,1].item(),self.Voxels_hit_indices_down[Event][j,2].item()
            if(down):
                ax[1].scatter(self.Voxels_centers[0,iy_down,iz_down,1],
                              self.Voxels_centers[0,iy_down,iz_down,2],
                              s=1,
                              c='orange',
                              label = 'Voxel hit down',alpha=.3)
            if(i==0):
                ax[1].legend()
                
        if(filename is not None):
            fig.savefig('Figures/'+filename+'.png')





def plot_setup(Detector_,VOI_)->None:
        
    from matplotlib.patches import Rectangle


    fig, ax = plt.subplots(1, 2, figsize=(20, 7))
    ax[0].set_xlabel('X [mm]')
    ax[0].set_ylabel('Z [mm]')
    ax[0].set_title(' XZ view')
    #plot xz frame
    ax[0].set_xlim(Detector_.Xmin,Detector_.Xmax)
    ax[0].set_ylim(Detector_.Zmin-100,Detector_.Zmax+100)
    
    for i in range(len(Detector_.Z)):
        ax[0].axhline(y=Detector_.Z[i],color='r',linestyle ='-')


    ax[0].add_patch( Rectangle((VOI_.Xmin,VOI_.Zmin),
    abs(VOI_.Xmax-VOI_.Xmin), abs(VOI_.Zmax-VOI_.Zmin),
    fc ='none', 
    ec ='g',
    lw = 4,
    label = 'VOI'))

    ax[0].legend()


    ax[1].set_xlabel('Y [mm]')
    ax[1].set_ylabel('Z [mm]')
    ax[1].set_title('YZ view')
    #plot xz frame
    ax[1].set_xlim(Detector_.Ymin,Detector_.Ymax)
    ax[1].set_ylim(Detector_.Zmin-100,Detector_.Zmax+100)

    for i in range(len(Detector_.Z)):
        ax[1].axhline(y=Detector_.Z[i],color='r',linestyle ='-')

    ax[1].add_patch( Rectangle((VOI_.Ymin,VOI_.Zmin),
    abs(VOI_.Ymax-VOI_.Ymin), abs(VOI_.Zmax-VOI_.Zmin),
    fc ='none', 
    ec ='g',
    lw = 4,
    label = 'VOI'))

    ax[1].legend()
    plt.show()
    

    
    
    
    
    
    
    
    
    
    
    
    
    
def norm_tensor(tensor)->torch.tensor:
    return (tensor-torch.min(tensor))/(torch.max(tensor)-torch.min(tensor))

def plot_scores(scores, title:str=None):
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    scoresXZ = torch.flip(torch.transpose(norm_tensor(torch.sum(scores,dim=1)),0,1),[0,1])
    scoresXZ[:,:].tolist()

    scoresZY = torch.flip(torch.transpose(norm_tensor(torch.sum(scores,dim=0)),0,1),[0,1])
    scoresZY[:,:].tolist()

    scoresXY = norm_tensor(torch.sum(scores,dim=2))
    scoresXY[:,:].tolist()
    
    fig1,ax=plt.subplots(nrows=1, ncols=1,figsize=(6.6,5))
    ax.set_title('XZ view')
    ax.set_xlabel('x [a.u]')
    ax.set_ylabel('z [a.u]')
    im = ax.imshow(scoresXZ,cmap='viridis')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig1.colorbar(im, ax=ax,cax=cax)
    
    fig2,ax=plt.subplots(nrows=1, ncols=1,figsize=(6.6,5))
    ax.set_title('YZ view')
    ax.set_xlabel('y [a.u]')
    ax.set_ylabel('z [a.u]')
    im2 = ax.imshow(scoresZY,cmap='viridis')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig2.colorbar(im2, ax=ax,cax=cax)
    plt.show()
    

    fig3,ax=plt.subplots(nrows=1, ncols=1,figsize=(10,6.6))
    im = ax.imshow(scoresXY)
    fig3.colorbar(im, ax=ax)
    ax.set_title('XY view')
    ax.set_xlabel('y')
    ax.set_ylabel('x')

    if(title is not None):
        fig1.savefig(title+'XZ_view')
        fig2.savefig(title+'YZ_view')
        fig3.savefig(title+'XY_view')

    plt.show()


        
        
